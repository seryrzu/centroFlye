# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

from collections import Counter
import logging
import os

import edlib
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from utils.os_utils import smart_makedirs

logger = logging.getLogger("centroFlye.mapping.mapping")


def map_queries(queries, targets,
                max_nloc_target,
                max_ntarget_locs,
                neutral_symbs,
                max_dist,
                mode='HW'):

    def map_query(query, target):
        if neutral_symbs is None or len(neutral_symbs) == 0:
            add_matches = None
        else:
            add_matches = []
            alphabet = set(query)
            for neutral_symb in neutral_symbs:
                for c in alphabet:
                    add_matches.append((c, neutral_symb))
        align = edlib.align(query,
                            target,
                            mode=mode,
                            task='locations',
                            k=max_dist,
                            additionalEqualities=add_matches)
        locs = align['locations']
        for i, (s, e) in enumerate(locs):
            locs[i] = (s, e+1)
        return locs

    all_locations = {i: {} for i in range(len(targets))}
    for q_id, query in queries.items():
        q_locs = {}
        for i, target in enumerate(targets):
            locs = map_query(query, target)
            if len(locs) > 0:
                q_locs[i] = locs
        if len(q_locs) <= max_ntarget_locs:
            for i, locs in q_locs.items():
                assert len(locs) > 0
                if len(locs) <= max_nloc_target:
                    all_locations[i][q_id] = locs
    return all_locations


def get_coverage(locations, target_len=None, outdir=None, title=None,
                 plot_close=True):
    coverage = Counter()
    for read_locs in locations.values():
        for s, e in read_locs:
            e += 1
            coverage[s] += 1
            coverage[e] -= 1

    if target_len is None:
        target_len = max(coverage) + 1
    else:
        assert target_len > max(coverage) + 1

    coverage_list = [0] * target_len
    for k, cov in coverage.items():
        coverage_list[k] = coverage[k]
    coverage = np.cumsum(coverage_list)


    if outdir is not None:
        smart_makedirs(outdir)
        txt_fn = os.path.join(outdir, 'coverage.txt')
        with open(txt_fn, 'w') as f:
            for i in range(len(coverage)):
                print(i, coverage[i], file=f)

        uncovered_bases = np.where(coverage == 0)[0]
        uncovered_intervals = []
        uncovered_intervals.append(uncovered_bases[0])
        for j, k in zip(uncovered_bases[:-1], uncovered_bases[1:]):
            if k-j != 1:
                uncovered_intervals.append(j)
                uncovered_intervals.append(k)
        uncovered_intervals.append(uncovered_bases[-1])
        assert len(uncovered_intervals) % 2 == 0
        uncovered_intervals = zip(uncovered_intervals[::2],
                                uncovered_intervals[1::2])
        uncovered_intervals = list(uncovered_intervals)
        uncovered_intervals_fn = os.path.join(outdir,
                                              'uncovered_intervals.txt')
        with open(uncovered_intervals_fn, 'w') as f:
            for s, e in uncovered_intervals:
                print(s, e, file=f)

        for s, e in uncovered_intervals:
            plt.axvspan(s, e, facecolor='gray', alpha=0.3)
        plt.plot(coverage)
        if title is not None:
            plt.title(title)
        plt.grid(True)
        plt.xlabel('monomers')
        plt.ylabel('coverage')
        pdf_fn = os.path.join(outdir, 'coverage.pdf')
        plt.savefig(pdf_fn, format='pdf')
        if plot_close:
            plt.close()

    return coverage
