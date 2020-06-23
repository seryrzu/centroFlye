# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

from config.config import config
import logging
import os
from subprocess import check_call

import pandas as pd

from monomers.monomer_db import MonomerDB
from monomers.monostring_set import MonoStringSet

from utils.bio import read_bio_seqs
from utils.os_utils import smart_makedirs

logger = logging.getLogger("centroFlye.sd_parser.sd_parser")


class SD_Report:
    def __init__(self, sd_report_fn, monomers_fn, sequences_fn, mode,
                 hpc=True, cluster=True, correct_hybrids=None):
        logger.info('Reading SD Report')

        assert mode in ['ont', 'hifi', 'assembly']
        if correct_hybrids is None:
            correct_hybrids = False  # correction of hybrids is disabled
            # if mode == 'assembly':
            #     correct_hybrids = False
            # else:
            #     correct_hybrids = True

        logger.info(f'Mode is {mode}')

        logger.info(f'    sd_report_fn = {sd_report_fn}')
        logger.info(f'    monomers_fn  = {monomers_fn}')
        logger.info(f'    sequences_fn = {sequences_fn}')
        monomer_db = MonomerDB.from_fasta_file(monomers_fn, cluster=cluster)

        if hpc:
            names = ['s_id', 'monomer',
                     's_st', 's_en',
                     'identity',
                     'sec_monomer', 'sec_identity',
                     'homo_monomer', 'homo_identity',
                     'sec_homo_monomer', 'sec_homo_identity',
                     'reliability'
                     ]
        else:
            names = ['s_id', 'monomer',
                     's_st', 's_en',
                     'identity',
                     'sec_monomer', 'sec_identity',
                     'reliability'
                     ]

        logger.info('Reading SD Report from csv')
        report = pd.read_csv(sd_report_fn, sep='\t',
                             header=None,
                             names=names)
        logger.info('Finished reading SD Report from csv')

        logger.info('Reading sequences')
        sequences = read_bio_seqs(sequences_fn)
        logger.info('Finished reading sequences')
        monostring_set = \
            MonoStringSet.from_sd_report(report=report,
                                         sequences=sequences,
                                         monomer_db=monomer_db,
                                         mode=mode,
                                         correct_hybrids=correct_hybrids)

        logger.info('Creating monostrings dict')
        logger.info('Finished creating monostrings dict')

        self.monomer_db = monomer_db
        self.report = report
        self.monostring_set = monostring_set


def run_SD(sequences_fn, monomers_fn, outdir='.',
           outfn='final_decomposition.tsv',
           n_threads=config["common"]["threads"]):
    logger.info(f'Running SD on')
    logger.info(f'\tSequences = {sequences_fn}')
    logger.info(f'\tMonomers = {monomers_fn}')
    logger.info(f'\tOutdir = {outdir}')
    smart_makedirs(outdir)
    outfn = os.path.join(outdir, outfn)
    if os.path.isfile(outfn):
        logger.info(f'File {outfn} exists. Reusing')
        return outfn
    cmd = f'{config["binaries"]["SD"]} {sequences_fn} {monomers_fn} ' + \
          f'-t {n_threads} -o {outfn}'
    logger.info(cmd)
    cmd = cmd.split(' ')
    check_call(cmd)
    return outfn
