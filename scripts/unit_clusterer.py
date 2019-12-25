# (c) 2019 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import argparse
import os
import statistics
import subprocess

from unit_extractor import get_period_info
from utils.bio import read_bio_seq, read_bio_seqs, write_bio_seqs
from utils.os_utils import smart_makedirs

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt


def get_units(input_dir):
    units = {}
    for f in os.scandir(input_dir):
        if f.is_dir():
            polished_fn = os.path.join(f.path, 'polished_2.fasta')
            unit = read_bio_seq(polished_fn)
            units[os.path.basename(f.path)] = unit
    return units


def select_median_seq(seqs):
    lens = [len(unit) for unit in seqs.values()]
    median_len = statistics.median(lens)
    for s_id in sorted(seqs.keys()):
        seq = seqs[s_id]
        if len(seq) == median_len:
            median_s_id = s_id
            median_seq = seqs[s_id]
            break
    return median_s_id, median_seq, median_len


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",
                        "--input",
                        help="Directory with read units",
                        required=True)
    # parser.add_argument("-r",
    #                     "--reads",
    #                     help="Input reads",
    #                     required=True)
    parser.add_argument("-o", "--outdir",
                        help="Output directory",
                        required=True)
    parser.add_argument("-b",
                        "--bin-size",
                        help="bin size",
                        type=int,
                        default=50)
    params = parser.parse_args()
    smart_makedirs(params.outdir)

    # reads = read_bio_seqs(params.reads)

    units = get_units(params.input)
    unit_lens = sorted(len(unit) for unit in units.values())
    periods, bin_convs, bin_left, bin_right = \
        get_period_info(unit_lens, bin_size=params.bin_size)

    # Currently support only one cluster
    filt_units = \
        {k: v for k, v in units.items() if bin_left <= len(v) <= bin_right}
    filt_units_fn = os.path.join(params.outdir, 'cluster_units.fasta')
    write_bio_seqs(filt_units_fn, filt_units)

    median_unit_id, median_unit, median_len = select_median_seq(filt_units)
    median_read_unit_fn = os.path.join(params.outdir, 'median_read_unit.fasta')
    write_bio_seqs(median_read_unit_fn,
                   {median_unit_id: median_unit})

    cmd = ['flye',
           f'--nano-raw', filt_units_fn,
           '--polish-target', median_read_unit_fn,
           '-i', 2,
           '-t', 50,
           '-o', params.outdir]
    cmd = [str(x) for x in cmd]
    subprocess.check_call(cmd)


if __name__ == "__main__":
    main()
