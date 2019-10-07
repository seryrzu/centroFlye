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


def hamming_distance(s1, s2, match_char=set()):
    #assert len(s1) == len(s2)
    nucl_ident = []
    for x, y in zip(s1, s2):
        if x in match_char or y in match_char:
            nucl_ident.append(0)
        else:
            nucl_ident.append(x != y)
    return sum(nucl_ident), len(nucl_ident)


def identity_shift(s1, s2, min_overlap, match_char=set()):
    best_identity, best_shift, best_hd, best_len = 0, None, None, None
    alt_shifts = []
    for shift in range(len(s1) - min_overlap):
        hd, cur_length = \
            hamming_distance(s1[shift:], s2, match_char=match_char)
        identity = 1 - hd / cur_length
        if identity == best_identity:
            alt_shifts.append(shift)
        if identity > best_identity:
            best_identity = identity
            best_shift = shift
            best_hd = hd
            best_len = cur_length
            alt_shifts = []
    return {'id': best_identity, 'shift': best_shift,
            'hd': best_hd, 'len': best_len,
            'alt_shifts': alt_shifts}


def OverlapAlignment(s1, s2, mismatch, sigma):
    n, m = len(s1) + 1, len(s2) + 1
    s1 = ' ' + s1
    s2 = ' ' + s2

    def cell(prev, char, score):
        return {'prev': prev, 'char': char, 'score': score}

    w = [[0] * m for i in range(n)]
    # for i in range(1, n):
    #     w[i][0] = w[i-1][0] - sigma
    for j in range(1, m):
        w[0][j] = w[0][j-1] - sigma

    for i in range(1, n):
        for j in range(1, m):
            a = w[i-1][j-1] + (1 if s1[i] == s2[j] else -mismatch)
            b, c = w[i-1][j] - sigma, w[i][j-1] - sigma
            w[i][j] = max(a, b, c)

    lrow_max = max(w[-1])
    jmax = [j for j in range(1, m) if w[-1][j] == lrow_max][0]
    a1, a2 = [], []
    i, j = n-1, jmax
    while i != 0 and j != 0:
        if w[i][j] == w[i-1][j-1] + (1 if s1[i] == s2[j] else -mismatch):
            a1.append(s1[i])
            a2.append(s2[j])
            i, j = i-1, j-1
        elif w[i][j] == w[i-1][j] - sigma:
            a1.append(s1[i])
            a2.append('-')
            i, j = i-1, j
        elif w[i][j] == w[i][j-1] - sigma:
            a1.append('-')
            a2.append(s2[j])
            i, j = i, j-1

    a1 = ''.join(a1[::-1])
    a2 = ''.join(a2[::-1])

    a1 = s1[1:(i+1)] + '|' + a1
    a2 = '-' * i + '|' + a2
    
    a1 += '|' + '-' * (m-jmax-1)
    a2 += '|' + s2[jmax+1:]
    
    assert len(a1) == len(a2)

    return w[n-1][jmax], a1, a2, i