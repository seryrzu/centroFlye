import argparse
import os

__filedir__ = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.append(os.path.join(__filedir__, '..'))


from utils.bio import read_bio_seq, write_bio_seqs, compress_homopolymer
from utils.os_utils import smart_makedirs
from ncrf_parser import NCRF_Report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ncrf", help="Input NCRF", required=True)
    parser.add_argument("--seq", help="Input sequence", required=True)
    parser.add_argument("--buf",
                        help="Buffer on the sides to include",
                        type=int,
                        default=20)
    parser.add_argument("--outdir", help="Output dir", required=True)
    params = parser.parse_args()

    smart_makedirs(params.outdir)
    ncrf_report = NCRF_Report(params.ncrf)
    input_seq = read_bio_seq(params.seq)
    all_mas = ncrf_report.get_motif_alignments()
    for seq_id, mas in all_mas.items():
        record = ncrf_report.records[seq_id]
        units = {}
        coords = {}
        al_start = record.r_st
        alignment = record.r_al.replace('-', '')
        start = 0
        for ma in mas:
            ma_st = ma.start
            ma_en = ma.end
            seq_al = record.r_al[ma_st:ma_en]
            seq = seq_al.replace('-', '')
            end = start + len(seq)
            seq_st = input_seq[al_start+start-params.buf:al_start+start]
            seq_en = input_seq[al_start+end:end+al_start+params.buf]
            seq = seq_st + seq + seq_en
            ma_id = f'{seq_id}|st_{start + al_start}|en_{end - 1 + al_start}'
            units[ma_id] = seq
            coords[ma_id] = (start + al_start, end + al_start)
            # print(input_seq[start+al_start:end+al_start])
            # print(seq[params.buf:-params.buf])
            assert input_seq[start+al_start-len(seq_st):end+al_start+len(seq_en)] == seq
            start = end
        outfile = os.path.join(params.outdir, f'{seq_id}.fasta')
        write_bio_seqs(outfile, units)


if __name__ == "__main__":
    main()
