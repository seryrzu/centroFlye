# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import argparse
import os
import sys

this_dirname = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(this_dirname, os.path.pardir))

from sd_parser.sd_parser import SD_Report
from standard_logger import get_logger
from submonomers.submonomer_db import SubmonomerDB
from submonomers.submonostring_set import SubmonoStringSet
from utils.os_utils import smart_makedirs, expandpath


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sd-report", required=True)
    parser.add_argument("--monomers", required=True)
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--coverage", type=int, required=True)
    params = parser.parse_args()
    params.sd_report = expandpath(params.sd_report)
    params.monomers = expandpath(params.monomers)
    params.sequences = expandpath(params.sequences)
    params.outdir = expandpath(params.outdir)
    return params


def main():
    params = parse_args()
    smart_makedirs(params.outdir)
    logfile = os.path.join(params.outdir, 'partitioning.log')
    logger = get_logger(logfile,
                        logger_name="centroFlye: submonopartion")

    logger.info('Reading SD Report')
    logger.info(f'    sd_report_fn = {params.sd_report}')
    logger.info(f'    monomers_fn  = {params.monomers}')
    logger.info(f'    sequences_fn = {params.sequences}')
    sd_report = SD_Report(sd_report_fn=params.sd_report,
                          sequences_fn=params.sequences,
                          monomers_fn=params.monomers)
    logger.info('Finished reading SD report')

    monostring_set = sd_report.monostring_set

    logger.info(f'Extracting submonomer db with coverage = {params.coverage}')
    submonomer_db = SubmonomerDB.from_monostring_set(monostring_set,
                                                     coverage=params.coverage)
    submonomer_db_fn = os.path.join(params.outdir, 'submonomer_db.fasta')
    logger.info('Finished extracting')
    logger.info(f'Exporing database to {submonomer_db_fn}')
    submonomer_db.to_fasta(filename=submonomer_db_fn)

    logger.info(f'Extracting submonosequences')
    submonostring_set = \
        SubmonoStringSet.from_monostringset(monostring_set, submonomer_db)

    submonostring_set_fn = os.path.join(params.outdir, 'submonostring_set.tsv')
    logger.info(f'Exporting submonoreads to {submonostring_set_fn}')
    submonostring_set.to_tsv(submonostring_set_fn)


if __name__ == "__main__":
    main()
