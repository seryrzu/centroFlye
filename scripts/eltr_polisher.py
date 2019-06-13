#(c) 2019 by Authors
#This file is a part of centroFlye program.
#Released under the BSD license (see LICENSE file)

import argparse
import os
import subprocess
import math
import edlib
from collections import defaultdict
from utils.os_utils import smart_makedirs
from ncrf_parser import NCRF_Report
from utils.bio import write_bio_seqs, read_bio_seq, read_bio_seqs, compress_homopolymer

current_dir = os.path.dirname(os.path.realpath(__file__))

def read_reported_positions(read_positions_fn):
    pos = {}
    with open(read_positions_fn) as f:
        for line in f:
            line = line.strip().split(' ')
            r_id, p = line[0], line[1]
            if p == 'None':
                p = None
            else:
                p = int(p)
            pos[r_id] = p
    return pos


class ELTR_Polisher:
    def __init__(self, params):
        self.params = params
        if not os.path.isfile(params.unit):
            raise FileNotFoundError(f"File {params.unit} is not found")
        self.unit = read_bio_seq(params.unit)
        self.ncrf_report = NCRF_Report(params.ncrf)
        self.motif_alignments = self.ncrf_report.get_motif_alignments()
        smart_makedirs(params.outdir)
        self.read_placement = read_reported_positions(params.read_placement)
        self.max_pos = self.params.max_pos
        self.min_pos = self.params.min_pos

    def map_pos2read(self):
        pos2read = defaultdict(list)
        for r_id, pos in self.read_placement.items():
            if pos is None or pos > self.max_pos:
                continue
            record = self.ncrf_report.records[r_id]
            ma = self.motif_alignments[r_id]
            for i in range(len(ma)):
                if self.min_pos <= pos + i <= self.max_pos:
                    pos2read[pos + i].append((r_id, i))
        return pos2read

    def export_read_units(self, pos2read):
        filenames = {}
        for pos in pos2read:
            outdir = os.path.join(self.params.outdir, f'pos_{pos}')
            units_fn = os.path.join(outdir, 'read_units.fasta')
            longest_read_unit_fn = os.path.join(outdir, 'longest_read_unit.fasta')
            smart_makedirs(outdir)
            seqs = {}
            longest_read_unit, template_read = "", None
            for (r_id, p) in pos2read[pos]:
                r_al = self.motif_alignments[r_id][p].r_al
                r_al = r_al.upper().replace('-', '')
                seqs[f'gen_pos={pos}|r_id={r_id}|r_pos={p}'] = r_al
                if len(r_al) >= len(longest_read_unit):
                    longest_read_unit = r_al
                    template_read = r_id
            write_bio_seqs(units_fn, seqs)
            write_bio_seqs(longest_read_unit_fn, {template_read: longest_read_unit})
            filenames[pos] = (units_fn, longest_read_unit_fn)
        return filenames

    def run_polishing(self, read_unit_filenames):
        polishing_wrapper_fn = os.path.join(current_dir, 'polishing_wrapper.py')
        min_pos = min(read_unit_filenames.keys())
        max_pos = max(read_unit_filenames.keys())
        for pos in range(min_pos, max_pos + 1):
            print(pos, max_pos)
            units_fn, longest_read_unit_fn = read_unit_filenames[pos]
            pos_dir = os.path.dirname(units_fn)
            cmd = ['flye',
                   f'--{self.params.error_mode}-raw', units_fn,
                   '--polish-target', self.params.unit,
                   '-i', self.params.num_iters,
                   '-t', self.params.num_threads,
                   '-o', pos_dir]
            if self.params.output_progress:
                cmd.append('--output-progress')
            cmd = [str(x) for x in cmd]
            print(' '.join(cmd))
            out = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
            polished_seq_fn = os.path.join(pos_dir, f'polished_{self.params.num_iters}.fasta')
            # if can't align to unit or alignment is too short -- align to the longest unit in reads
            if out.decode("utf-8").find("No reads were aligned during polishing") != -1 or \
                    len(read_bio_seq(polished_seq_fn)) < 0.7 * len(self.unit):
                cmd = ['flye',
                    f'--{self.params.error_mode}-raw', units_fn,
                    '--polish-target', longest_read_unit_fn,
                    '-i', self.params.num_iters,
                    '-t', self.params.num_threads,
                    '-o', pos_dir]
                if self.params.output_progress:
                    cmd.append('--output-progress')
                cmd = [str(x) for x in cmd]
                print(' '.join(cmd))
                subprocess.check_call(cmd)

    def read_polishing(self, read_unit_filenames):
        min_pos = min(read_unit_filenames.keys())
        max_pos = max(read_unit_filenames.keys())
        polished_seqs = {}
        final_sequences = {}
        for i in range(1, self.params.num_iters + 1):
            for pos, (units_fn, longest_read_unit_fn) in read_unit_filenames.items():
                pos_dir = os.path.dirname(units_fn)
                polished_seq_fn = os.path.join(pos_dir, f'polished_{i}.fasta')
                polished_seq = read_bio_seq(polished_seq_fn)
                polished_seqs[pos] = polished_seq
            final_sequence = [polished_seqs[pos] for pos in range(min_pos, max_pos + 1)]
            final_sequence = ''.join(final_sequence)
            final_sequences[i] = final_sequence
        return final_sequences

    def compare_polished_sequences(self, final_sequences):
        report_fn = os.path.join(self.params.outdir, 'report.txt')
        with open(report_fn, 'w') as f:
            for i in range(1, self.params.num_iters):
                seq_i, seq_i1 = final_sequences[i], final_sequences[i+1]
                alignment = edlib.align(seq_i, seq_i1)
                print(f'Alignment polishing seq {i} vs {i+1}:', file=f)
                print(alignment, file=f)

                hpc_seq_i, hpc_seq_i1 = compress_homopolymer(final_sequences[i]), compress_homopolymer(final_sequences[i+1])
                alignment = edlib.align(hpc_seq_i, hpc_seq_i1)
                print(f'Alignment homopolymer compressed polishing seq {i} vs {i+1}:', file=f)
                print(alignment, file=f)


    def export_results(self, final_sequences):
        for i in range(1, self.params.num_iters + 1):
            final_sequence = final_sequences[i]
            final_sequence_hpc = compress_homopolymer(final_sequence)

            final_fn = os.path.join(self.params.outdir, f'final_sequence_{i}.fasta')
            write_bio_seqs(final_fn, {f'polished_repeat_{i}': final_sequence})

            final_hpc_fn = os.path.join(self.params.outdir, f'final_sequence_hpc_{i}.fasta')
            write_bio_seqs(final_hpc_fn, {f'polished_repeat_{i}': final_sequence_hpc})


    def run(self):
        pos2read = self.map_pos2read()
        read_unit_filenames = self.export_read_units(pos2read)
        self.run_polishing(read_unit_filenames)
        final_sequences = self.read_polishing(read_unit_filenames)
        self.compare_polished_sequences(final_sequences)
        self.export_results(final_sequences)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--read-placement", required=True)
    parser.add_argument("--unit", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--ncrf", required=True)
    parser.add_argument("--error-mode", default="pacbio")
    parser.add_argument("--num-iters", default=2, type=int)
    parser.add_argument("--num-threads", default=16, type=int)
    parser.add_argument("--output-progress", action='store_false')
    parser.add_argument("--min-pos", type=int, default=0)
    parser.add_argument("--max-pos", type=int, default=math.inf)
    params = parser.parse_args()

    ELTR_Polisher(params).run()


if __name__ == "__main__":
    main()
