# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import logging

import pandas as pd

from monomers.monomer_db import MonomerDB
from monomers.monostring import MonoString

from utils.bio import read_bio_seqs
from utils.os_utils import expandpath

logger = logging.getLogger("centroFlye.sd_parser.sd_parser")

class SD_Report:
    def __init__(self, sd_report_fn, monomers_fn, sequences_fn, hpc=True):
        sd_report_fn = expandpath(sd_report_fn)
        monomers_fn = expandpath(monomers_fn)
        sequences_fn = expandpath(sequences_fn)

        logger.info('Reading SD Report')
        logger.info(f'    sd_report_fn = {sd_report_fn}')
        logger.info(f'    monomers_fn  = {monomers_fn}')
        logger.info(f'    sequences_fn = {sequences_fn}')
        monomer_db = MonomerDB.from_fasta_file(monomers_fn)

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


        logger.info('Creating monostrings dict')
        monostrings = {}
        for seq_id, seq_record in report.groupby('s_id'):
            seq_record = seq_record.sort_values(by=['s_st'])
            sequence = sequences[seq_id]
            monostring = MonoString.from_sd_record(seq_id=seq_id,
                                                   monomer_db=monomer_db,
                                                   sd_record=seq_record,
                                                   sequence=sequence)
            monostrings[seq_id] = monostring
        logger.info('Finished creating monostrings dict')

        self.monomer_db = monomer_db
        self.report = report
        self.monostrings = monostrings

