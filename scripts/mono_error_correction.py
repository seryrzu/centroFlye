from collections import Counter
from copy import deepcopy

import numpy as np

from utils.bio import hamming_distance, min_cyclic_shift
from debruijn_graph import DeBruijnGraph, get_frequent_kmers
from sd_parser_83640e3 import get_stats


def get_ma(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)


def filter_lowercaserich_reads(monoreads, max_lowercase=0.1):
    filtered_strings = {}
    for r_id, string in monoreads.items():
        lowercase = [s.islower() for s in string]
        if np.mean(lowercase) <= max_lowercase:
            filtered_strings[r_id] = string
    return filtered_strings


def trim_read(monoread, max_gap, ma_window, gap_symb='?'):
    is_gap = [c == gap_symb for c in monoread]
    ma = get_ma(is_gap, N=ma_window)
    left = 0
    while left < len(ma) and ma[left] > max_gap:
        left += 1
    right = len(ma) - 1
    while right >= 0 and ma[right] > max_gap:
        right -= 1
    monoread.trim_read(left, right+ma_window+1)
    monoread.strip()
    return monoread


def trim_reads(monoreads, max_gap=0.2, ma_window=30):
    trimmed_reads = {}
    for r_id, monoread in monoreads.items():
        trimmed_read = \
            trim_read(monoread, max_gap=max_gap, ma_window=ma_window)
        trimmed_reads[r_id] = trimmed_read
    return trimmed_reads


def cut_gaprich_reads(monoreads, max_gap=0.05, min_length=100, gap_symb='?'):
    processed_reads = {}
    cut_cnt = 0
    total_parts_cnt = 0
    for r_id, monoread in monoreads.items():
        if len(monoread) == 0:
            processed_reads[r_id] = monoread
            continue
        gap_cnt = Counter(monoread)[gap_symb]
        gap_prop = gap_cnt / len(monoread)
        if gap_prop <= max_gap:
            processed_reads[(r_id, 0)] = monoread
        else:
            parts = monoread.split(c=gap_symb, min_length=min_length)
            part_cnt = sum(len(part.string) > 0 for part in parts.values())
            if part_cnt:
                cut_cnt += 1
                total_parts_cnt += part_cnt
            processed_reads.update(parts)
    return processed_reads, cut_cnt, total_parts_cnt


def correct_gaps(monostrings, max_gap=0.3, nhor=1,
                 k=3, min_mult=5000,
                 gap_symb='?'):
    frequent_kmers, frequent_kmers_read_pos = \
        get_frequent_kmers(monostrings, k=k, min_mult=min_mult)
    db = DeBruijnGraph(k=k)
    db.add_kmers(frequent_kmers, coverage=frequent_kmers)

    hors, _ = db.get_contigs()
    hors = [min_cyclic_shift(hor) for hor in hors]
    hors.sort()
    corrected_strings = {}
    for r_id, monostring in monostrings.items():
        for single_hor in hors:
            for i_nhor in range(1, nhor+1):
                hor = single_hor * i_nhor
                hor_len = len(hor)
                for i in range(len(monostring)-hor_len+1):
                    kmer = monostring[i:i+hor_len]
                    gap_cnt = Counter(kmer)[gap_symb]
                    if gap_cnt == 0 or gap_cnt / hor_len > max_gap:
                        continue
                    hd, _ = hamming_distance(kmer, hor,
                                             match_char=set(gap_symb))
                    if hd == 0:
                        monostring[i:i+hor_len] = list(hor)
        monostring.assert_validity()
    return monostrings


def error_correction(monoreads, inplace=True, verbose=False, hor_correction=True):
    if not inplace:
        monoreads = deepcopy(monoreads)

    if verbose:
        print('Stats for uncorrected reads')
        get_stats(monoreads, verbose=True, return_stats=False)
    filtered_reads = filter_lowercaserich_reads(monoreads)

    if verbose:
        print('\nStats for filtered reads')
        get_stats(filtered_reads, verbose=True, return_stats=False)

    trimmed_monoreads = trim_reads(filtered_reads)
    if verbose:
        print('\nStats for trimmed+filtered reads')
        get_stats(trimmed_monoreads, verbose=True, return_stats=False)

    cut_monoreads, cut_cnt, total_parts_cnt = \
        cut_gaprich_reads(trimmed_monoreads)
    if verbose:
        print('\nStats for cut_gaprich_reads+trimmed+filtered reads')
        print(f'# cut reads = {cut_cnt}')
        print(f'# total cut parts = {total_parts_cnt}')
        get_stats(cut_monoreads, verbose=True, return_stats=False)
    corrected_reads = cut_monoreads

    if hor_correction:
        hor_corrected = correct_gaps(cut_monoreads)
        corrected_reads = hor_corrected
        if verbose:
            print('\nStats for hor_corrected+cut_gaprich_reads+trimmed+filtered reads')
            get_stats(hor_corrected, verbose=True, return_stats=False)

    return corrected_reads
