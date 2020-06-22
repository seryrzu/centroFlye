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

from config.config import config, copy_config
from sd_parser.sd_parser import SD_Report
from sequence_graph.idb_graph import get_idb_monostring_set
from sequence_graph.db_graph_scaffolding import monoscaffolds2scaffolds
from utils.git import get_git_revision_short_hash
from utils.os_utils import smart_makedirs


def parse_args():
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
    def __init__(self, logger):
        self.logger = logger

    def run(self, params):
        sd_report = SD_Report(sd_report_fn=params.sd_report,
                              monomers_fn=params.monomers,
                              sequences_fn=params.reads,
                              mode=params.mode)
        idb_outdir = os.path.join(params.outdir, 'idb')
        dbs, all_frequent_kmers = \
            get_idb_monostring_set(string_set=sd_report.monostring_set,
                                   mink=config['idb']['mink'],
                                   maxk=config['idb']['maxk'],
                                   outdir=idb_outdir,
                                   mode=params.mode)
        db = dbs[config['idb']['maxk']]
        scaffolding_outdir = os.path.join(params.outdir, 'scaffolding')
        monoscaffolds2scaffolds(db, sd_report.monostring_set,
                                outdir=scaffolding_outdir)


def main():
    params = parse_args()
    smart_makedirs(params.outdir)
    copy_config(params.outdir)

    logfn = os.path.join(params.outdir, 'centroFlye.log')
    logger = get_logger(logfn,
                        logger_name='centroFlye')

    logger.info(f'cmd: {sys.argv}')
    logger.info(f'git hash: {get_git_revision_short_hash()}')

    centroFlye(logger).run(params)


if __name__ == "__main__":
    main()
