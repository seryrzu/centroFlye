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

from assembly_utils.assembly_utils import map_monoreads_to_monoassembly
from config.config import config, copy_config
from sd_parser.sd_parser import SD_Report, run_SD
from sequence_graph.db_graph_scaffolding import monoscaffolds2scaffolds
from sequence_graph.idb_graph import get_idb_monostring_set
from sequence_graph.path_graph import PathDeBruijnGraph
from utils.git import get_git_revision_short_hash
from utils.os_utils import expandpath, smart_makedirs


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
                        help='Path to SD report')
    parser.add_argument('-t', '--threads',
                        help='Number of threads',
                        type=int,
                        default=config['common']['threads'])
    parser.add_argument('--mode',
                        help='ont/hifi',
                        default='ont')
    parser.add_argument('-a', '--assembly',
                        help='Assembly for comparison')
    params = parser.parse_args()
    if params.mode != 'ont':
        print('ERROR. Only the ont mode is supported so far')
        sys.exit(2)

    params.reads = expandpath(params.reads)
    params.monomers = expandpath(params.monomers)
    params.outdir = expandpath(params.outdir)
    if params.sd_report is not None:
        params.sd_report = expandpath(params.sd_report)
    if params.assembly is not None:
        params.assembly = expandpath(params.assembly)

    # override config TODO make singleton
    # config['common']['threads'] = params.threads

    return params


class centroFlye:
    def __init__(self, logger):
        self.logger = logger

    def run(self, params):
        if params.sd_report is None:
            reads_sd_report_outdir = os.path.join(params.outdir,
                                                  'SD_reads')
            reads_sd_report_fn = run_SD(sequences_fn=params.reads,
                                        monomers_fn=params.monomers,
                                        outdir=reads_sd_report_outdir,
                                        n_threads=params.threads)
        else:
            reads_sd_report_fn = params.sd_report

        sd_report = SD_Report(sd_report_fn=reads_sd_report_fn,
                              monomers_fn=params.monomers,
                              sequences_fn=params.reads,
                              mode=params.mode)
        monoreads_set = sd_report.monostring_set

        if params.assembly is not None:
            assembly_sd_report_outdir = os.path.join(params.outdir,
                                                     'SD_assembly')
            assembly_sd_report_fn = run_SD(sequences_fn=params.assembly,
                                           monomers_fn=params.monomers,
                                           outdir=assembly_sd_report_outdir,
                                           n_threads=params.threads)
            assembly_sd_report = SD_Report(sd_report_fn=assembly_sd_report_fn,
                                           monomers_fn=params.monomers,
                                           sequences_fn=params.assembly,
                                           mode='assembly')
            assembly_stats_dir = os.path.join(params.outdir, 'assembly_stats')
            monoassembly = assembly_sd_report.monostring_set

            reads2assembly_dir = os.path.join(assembly_stats_dir,
                                              'monoreads2monoassembly')
            map_monoreads_to_monoassembly(monoreads_set=monoreads_set,
                                          monoassembly_set=monoassembly,
                                          outdir=reads2assembly_dir)

        else:
            monoassembly = None

        idb_outdir = os.path.join(params.outdir, 'idb')
        dbs, all_frequent_kmers = \
            get_idb_monostring_set(string_set=monoreads_set,
                                   mink=config['idb']['mink'],
                                   maxk=config['idb']['maxk'],
                                   assembly=monoassembly,
                                   outdir=idb_outdir,
                                   mode=params.mode)
        db = dbs[config['idb']['maxk']]
        path_db_outdir = os.path.join(params.outdir, 'path_db')
        path_db = PathDeBruijnGraph.from_mono_db(db=db,
                                                 monostring_set=monoreads_set,
                                                 assembly=monoassembly,
                                                 k=config['path_db']['k'],
                                                 outdir=path_db_outdir)
        scaffolding_outdir = os.path.join(params.outdir, 'scaffolding')
        monoscaffolds2scaffolds(path_db, sd_report.monostring_set,
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
    logger.info(f'config: {config}')

    centroFlye(logger).run(params)


if __name__ == "__main__":
    main()
