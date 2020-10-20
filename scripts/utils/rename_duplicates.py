import argparse
from collections import Counter

from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True)
    parser.add_argument("-o", required=True)
    params = parser.parse_args()

    cnt = Counter()
    with open(params.o, 'w') as o:
        for seq_record in SeqIO.parse(open(params.i), 'fasta'):
            cnt[seq_record.id] += 1
            seq_record.id += f'_{cnt[seq_record.id]}'
            print(seq_record.description)
            seq_record.description = ''
            SeqIO.write(seq_record, o, 'fasta')


if __name__ == "__main__":
    main()
