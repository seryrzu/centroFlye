# (c) 2019 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

from Bio import SeqIO
import numpy as np
from itertools import groupby
import re


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


def parse_cigar(cigar, s1=None, s2=None):
    parsed_cigar = []
    st = 0
    cnt = dict.fromkeys(list("=XID"), 0)
    for mo in re.finditer(r'=|X|I|D', cigar):
        group = mo.group()
        pos = mo.start()
        region_len = int(cigar[st:pos])
        parsed_cigar.append((region_len, group))
        cnt[group] += region_len
        st = pos + 1
    if s1 is None or s2 is None:
        return parsed_cigar, cnt

    a1, a2 = [], []
    i1, i2 = 0, 0
    for region_len, group in parsed_cigar:
        if group in '=X':
            new_s1 = s1[i1:i1+region_len]
            new_s2 = s2[i2:i2+region_len]
            if group == '=':
                assert new_s1 == new_s2
            a1 += new_s1
            a2 += new_s2
            i1 += region_len
            i2 += region_len
        elif group == 'D':
            a1 += '-' * region_len
            a2 += s2[i2:i2+region_len]
            i2 += region_len
        elif group == 'I':
            a2 += '-' * region_len
            a1 += s1[i1:i1+region_len]
            i1 += region_len

    a1 = ''.join(a1)
    a2 = ''.join(a2)
    return parsed_cigar, cnt, a1, a2


assert parse_cigar('89=1X6=3X76=') == ([(89, '='), (1, 'X'), (6, '='), (3, 'X'), (76, '=')],
                                       {'=': 171, 'X': 4, 'I': 0, 'D': 0})
