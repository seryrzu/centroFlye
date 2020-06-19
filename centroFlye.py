#!/usr/bin/env python3

# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import argparse
import os
import sys

this_dirname = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(this_dirname, 'scripts'))

from standard_logger import get_logger

from config.config import get_config
from sd_parser.sd_parser import SD_Report
from sequence_graph.idb_graph import get_idb_monostring_set
from scripts.utils.bio import read_bio_seqs
from utils.os_utils import smart_makedirs


def parse_args(config):
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--reads',
                        help='Path to centromeric reads',
                        required=True)
    parser.add_argument('-m', '--monomers',
                        help='Path to fasta with monomers used to call SD',
                        required=True)
    parser.add_argument('-o', '--outdir',
                        help='Output directory',
                        required=True)
    parser.add_argument('-s', '--sd-report',
                        help='Path to SD report',
                        required=True)
    parser.add_argument('-t', '--threads',
                        help='Number of threads',
                        type=int,
                        default=config['common']['threads'])
    parser.add_argument('--mode',
                        help='ont/hifi',
                        default='ont')
    params = parser.parse_args()
    if params.mode != 'ont':
        print('ERROR. Only the ont mode is supported so far')
        sys.exit(2)
    return params


class centroFlye:
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def run(self, params):
        sd_report = SD_Report(sd_report_fn=params.sd_report,
                              monomers_fn=params.monomers,
                              sequences_fn=params.reads,
                              mode=params.mode)
        idb_outdir = os.path.join(params.outdir, 'idb')
        dbs, all_frequent_kmers = \
            get_idb_monostring_set(string_set=sd_report.monostring_set,
                                   mink=self.config['idb']['mink'],
                                   maxk=self.config['idb']['maxk'],
                                   outdir=idb_outdir,
                                   mode=params.mode)


def main():
    config = get_config()
    params = parse_args(config)
    smart_makedirs(params.outdir)
    logfn = os.path.join(params.outdir, 'centroFlye.log')
    logger = get_logger(logfn,
                        logger_name='centroFlye')
    centroFlye(config, logger).run(params)


if __name__ == "__main__":
    main()
