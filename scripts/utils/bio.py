# (c) 2019 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

from Bio import SeqIO
import numpy as np
from itertools import groupby


def read_bio_seq(filename):
    seqs = read_bio_seqs(filename)
    return str(list(seqs.values())[0])


def read_bio_seqs(filename):
    form = filename.split('.')[-1]
    if form == 'fa' or form == 'fna':
        form = 'fasta'
    elif form == 'fq':
        form = 'fastq'
    seqs = SeqIO.parse(filename, format=form)
    seqs = {seq.id: str(seq.seq) for seq in seqs}
    return seqs


trans = str.maketrans('ATGCatgc-', 'TACGtacg-')
def RC(s):
    return s.translate(trans)[::-1]


def write_bio_seqs(filename, seqs):
    with open(filename, 'w') as f:
        for seq_id, seq in seqs.items():
            print(f'>{seq_id}', file=f)
            print(seq, file=f)


def reverse_seq(seq):
    rev_seq = list(seq[::-1])
    for i in range(len(rev_seq)):
        if rev_seq[i] == 'R':
            continue
        sign = rev_seq[i][0]
        if sign == '+':
            rev_seq[i] = '-' + rev_seq[i][1:]
        elif sign == '-':
            rev_seq[i] = '+' + rev_seq[i][1:]
        else:
            raise ValueError
    return rev_seq


def gen_random_seq(length, bases=list("ACGT")):
    indexes = np.random.randint(len(bases), size=length)
    sequence = ''.join([bases[i] for i in indexes])
    return sequence


def compress_homopolymer(seq):
    return ''.join(x[0] for x in groupby(list(seq)))
