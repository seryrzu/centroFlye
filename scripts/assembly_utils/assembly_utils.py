# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

from collections import Counter
import os

import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from utils.os_utils import smart_makedirs


def get_coverage(locations, target_len=None, outdir=None, title=None):
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
        pdf_fn = os.path.join(outdir, 'coverage.pdf')
        txt_fn = os.path.join(outdir, 'coverage.txt')
        with open(txt_fn, 'w') as f:
            for i in range(len(coverage)):
                print(i, coverage[i], file=f)
        plt.plot(coverage)
        if title is not None:
            plt.title(title)
        plt.xlabel('monomers')
        plt.ylabel('coverage')
        plt.savefig(pdf_fn, format='pdf')

    return coverage


# def map_sequences()
