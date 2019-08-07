# (c) 2019 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import argparse
import os
from subprocess import call, Popen, PIPE
from utils.bio import read_bio_seqs, read_bio_seq, write_bio_seqs
from utils.various import chunks2
from utils.os_utils import smart_makedirs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reads",
                        help="Path to centromeric reads in fasta format",
                        required=True)
    parser.add_argument("--repeat",
                        help="Path to the unit sequence",
                        required=True)
    parser.add_argument("-t",
                        "--threads",
                        help="Number of threads",
                        type=int,
                        default=30)
    parser.add_argument("-o",
                        "--outdir",
                        help="Output directory",
                        required=True)
    parser.add_argument("--ncrf-bin",
                        help="Path to binary of NCRF",
                        default='NCRF')
    params = parser.parse_args()
    smart_makedirs(params.outdir)

    repeat = read_bio_seq(params.repeat)

    reads = read_bio_seqs(params.reads)
    reads_split = chunks2(list(reads.keys()), params.threads)
    reads_chunks_fn = {}
    for i in range(len(reads_split)):
        reads_chunk = {k: reads[k] for k in reads_split[i]}
        outdir = os.path.join(params.outdir, 'split_reads')
        smart_makedirs(outdir)
        reads_fn = os.path.join(outdir, f'split_reads_{i}.fasta')
        reads_chunks_fn[i] = reads_fn
        write_bio_seqs(reads_fn, reads_chunk)

    ps = []
    ncrf_reports_fn = []
    for i, fn in reads_chunks_fn.items():
        outdir = os.path.join(params.outdir, 'ncrf_report')
        smart_makedirs(outdir)
        ncrf_report_fn = os.path.join(outdir, f'report_{i}.ncrf')
        with open(ncrf_report_fn, 'w') as f:
            p1 = Popen(['cat', fn], stdout=PIPE)
            p2 = Popen([params.ncrf_bin, f'unit:{repeat}'],
                       stdin=p1.stdout, stdout=f)
            ps.append(p2)
        ncrf_reports_fn.append(ncrf_report_fn)
    for p in ps:
        p.wait()

    final_report_fn = os.path.join(params.outdir, 'report.ncrf')
    with open(final_report_fn, 'w') as f:
        cmd1 = ['cat'] + ncrf_reports_fn
        p1 = Popen(cmd1, stdout=PIPE)
        cmd2 = f"grep -v -E end-of-file".split(' ')
        p2 = Popen(cmd2, stdin=p1.stdout, stdout=f)
        p2.wait()

    cmd = f'sed -i s/unit/{repeat}/g {final_report_fn}'
    call(cmd.split(' '))


if __name__ == "__main__":
    main()
