# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)


import argparse
from collections import defaultdict
import logging
import os
import sys

this_dirname = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(this_dirname, os.path.pardir))

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from scipy.ndimage.filters import maximum_filter1d, minimum_filter1d

from standard_logger import get_logger
from utils.bio import read_bio_seq
from utils.git import get_git_revision_short_hash
from utils.os_utils import smart_makedirs

logger = logging.getLogger("centroFlye.monomers.monomer_db")


def max_filter1d_valid(a, W):
    hW = (W-1)//2  # Half window size
    return maximum_filter1d(a, size=W)[hW:-hW]


def min_filter1d_valid(a, W):
    hW = (W-1)//2  # Half window size
    return minimum_filter1d(a, size=W)[hW:-hW]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reads-fai", required=True)
    parser.add_argument("--paf", required=True)
    parser.add_argument("--right", required=True, type=int)
    parser.add_argument("--left", required=True, type=int)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--ref", required=True)
    params = parser.parse_args()

    smart_makedirs(params.outdir)
    logfn = os.path.join(params.outdir, 'substr_surv_paf.log')
    global logger
    logger = get_logger(logfn,
                        logger_name='centroFlye: substr survival paf')

    logger.info(f'cmd: {sys.argv}')
    logger.info(f'git hash: {get_git_revision_short_hash()}')

    logger.info(f'Reading reference from {params.ref}')
    ref = read_bio_seq(params.ref)
    assert len(ref) == params.right - params.left + 1
    uncompr2compr = {params.left: params.left}
    for i in range(1, len(ref)):
        uncompr2compr[i+params.left] = uncompr2compr[i-1+params.left]
        if ref[i] != ref[i-1]:
            uncompr2compr[i+params.left] += 1
    logger.info(f'Finished reading reference')

    logger.info(f'Reading read ids from {params.reads_fai}')
    read_ids = set()
    with open(params.reads_fai) as f:
        for line in f:
            line = line.strip().split('\t')
            read_ids.add(line[0])
    logger.info(f'Finished reading read ids')

    logger.info(f'Reading coords from {params.paf}')
    coords = []
    with open(params.paf) as f:
        for line in f:
            line = line.strip().split('\t')
            r_id, st, en = line[0], int(line[7]), int(line[8])
            if r_id not in read_ids or params.left > st or en > params.right:
                continue
            st = uncompr2compr[st]
            en = uncompr2compr[en]
            coords.append((st, en))
    coords.sort(key=lambda x: (x[0], -x[1]))
    logger.info(f'Finished reading coords')

    surv = defaultdict(int)
    rightmost_end = 0
    for st, en in coords:
        rightmost_end = max(en, rightmost_end)
        surv[st] = rightmost_end - st

    substr_surv_tsv_fn = os.path.join(params.outdir, 'substr_surv.tsv')
    with open(substr_surv_tsv_fn, 'w') as f:
        for i in range(min(surv)+1, max(surv)):
            surv[i] = max(surv[i],
                          surv[i-1]-1,
                          0)
            print(f'{i}\t{surv[i]}', file=f)

    substr_surv_pdf_fn = os.path.join(params.outdir, 'substr_surv.pdf')
    surv = list(surv.values())

    W = 10000
    plt.plot(max_filter1d_valid(surv, W=W))
    plt.plot(min_filter1d_valid(surv, W=W))

    plt.xlabel('reference')
    plt.ylabel('max survived substring length')
    plt.legend([f'Max', f'Min'])
    plt.grid()
    plt.ylim((0, 30000))
    plt.savefig(substr_surv_pdf_fn, format='pdf')
    plt.close()


if __name__ == "__main__":
    main()
