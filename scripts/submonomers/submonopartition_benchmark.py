# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)


import argparse
import os
import sys

from joblib import Parallel, delayed
import numpy as np

import matplotlib.pyplot as plt

this_dirname = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(this_dirname, os.path.pardir))

from sd_parser.sd_parser import SD_Report
from standard_logger import get_logger
from submonomers.submonomer_db import SubmonomerDB
from submonomers.submonostring import SubmonoString
from submonomers.submonomers_extraction_benchmark import get_coverage
from submonomers.submonostring_set import SubmonoStringSet,\
                                          CorrectedSubmonoStringSet
from utils.os_utils import smart_makedirs



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sd-report-reads", required=True)
    parser.add_argument("--sd-report-assembly", required=True)
    parser.add_argument("--monomers", required=True)
    parser.add_argument("--reads", required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--outdir", required=True)
    params = parser.parse_args()
    return params


def get_submonostrings(params, logger):
    logger.info('Reading SD reports')
    sd_report_reads = SD_Report(sd_report_fn=params.sd_report_reads,
                                sequences_fn=params.reads,
                                monomers_fn=params.monomers)
    sd_report_assembly = SD_Report(sd_report_fn=params.sd_report_assembly,
                                   sequences_fn=params.assembly,
                                   monomers_fn=params.monomers)
    logger.info('Finished reading SD reports')
    monoassembly = \
        list(sd_report_assembly.monostring_set.monostrings.values())[0]
    monoreads_set = sd_report_reads.monostring_set

    coverage = get_coverage(monoassembly, monoreads_set)
    logger.info(f'Estimated coverage = {coverage:.5}x')

    logger.info(f'Extracting submonomers from the assembly')
    submonomers_assembly = \
        SubmonomerDB.from_monostring_set(sd_report_assembly.monostring_set,
                                         coverage=1)

    logger.info(f'Extracting submonoreads')
    submonoread_set = SubmonoStringSet.from_monostringset(monoreads_set,
                                                          submonomers_assembly)

    logger.info(f'Extracting submonoassembly')
    submonoassembly = SubmonoString.from_monostring(monoassembly,
                                                    submonomers_assembly)
    return submonoread_set, submonoassembly


def map_reads(submonoread_set, submonoassembly, n_jobs=30):
    def map_read(seq_id, raw_submonoread_str, raw_submonoassembly_str,
                 gap_symbs):
        read_len = len(raw_submonoread_str)
        best_hd = read_len + 1
        best_pos = []
        best_diff = []
        for i in range(len(raw_submonoassembly_str)-read_len+1):
            a_smsubstr = raw_submonoassembly_str[i:i+read_len]
            diff = [(r, a) for r, a in zip(raw_submonoread_str, a_smsubstr)
                    if r != a and r not in gap_symbs]
            hd = len(diff)
            if hd < best_hd:
                best_pos = [i]
                best_hd = hd
                best_diff = [diff]
            elif hd == best_hd:
                best_pos.append(i)
                best_diff.append(diff)
        return seq_id, best_hd, best_pos, best_diff

    results = Parallel(n_jobs=n_jobs, verbose=60)\
              (delayed(map_read)(submonoread.seq_id,
                                 submonoread.raw_submonostring,
                                 submonoassembly.raw_submonostring,
                                 submonoassembly.get_gap_symbols())
               for submonoread in submonoread_set.submonostrings.values())
    mappings = {}
    for seq_id, best_hd, best_pos, best_diff in results:
        if len(best_pos) == 1:
            pos = best_pos[0]
            mappings[seq_id] = pos
    return mappings


def estimate_rates_over_assembly(mappings,
                                 submonoassembly,
                                 submonoread_set,
                                 min_ident=0.3,
                                 min_cov=5):
    assembly_len = len(submonoassembly)
    errors = np.zeros(assembly_len)
    submonoequiv = np.zeros(assembly_len)
    coverage = np.zeros(assembly_len)

    gap_symbols = submonoassembly.get_gap_symbols()
    for r_id, p in mappings.items():
        submonoread = submonoread_set[r_id]
        read_len = len(submonoread)
        a_smsubstr = submonoassembly[p:p+read_len]
        diff = [(r, a) for r, a in zip(submonoread, a_smsubstr)
                if r != a and r not in gap_symbols]
        hd = len(diff)
        if hd / len(submonoread) > min_ident:
            continue
        for i, (cr, ca) in enumerate(zip(submonoread, a_smsubstr)):
            ia = i + p
            coverage[ia] += 1
            if cr != ca and cr not in gap_symbols:
                errors[ia] += 1
            if cr == submonoread.submonoequiv_symb:
                submonoequiv[ia] += 1

    error_rate = errors / coverage
    error_rate[coverage < min_cov] = 0

    submonoequiv_rate = submonoequiv / coverage
    submonoequiv_rate[coverage < min_cov] = 0

    return coverage, error_rate, submonoequiv_rate


def run_rate_stats(mappings, submonoassembly, submonoread_set,
                   cor_submonoread_set, outdir, logger):
    coverage, error_rate, submonoequiv_rate = \
        estimate_rates_over_assembly(mappings,
                                     submonoassembly,
                                     submonoread_set)
    plt.plot(coverage)
    plt.title('Coverage with uniquely mapped reads')
    plt.savefig(os.path.join(outdir, 'coverage.pdf'), format='pdf')
    plt.close()

    zero_cov = np.where(coverage == 0)[0]
    logger.info(f'Positions of zero coverage {zero_cov}')

    plt.plot(error_rate)
    plt.title('Error rate')
    plt.savefig(os.path.join(outdir, 'error_rate_assembly.pdf'), format='pdf')
    plt.close()

    high_error_rate = np.where(np.logical_and(error_rate > 0.9, coverage > 5))
    logger.info(f'Positions of high error rate {high_error_rate}')

    _, _, cor_submonoequiv_rate = \
        estimate_rates_over_assembly(mappings,
                                     submonoassembly,
                                     cor_submonoread_set)

    plt.plot(submonoequiv_rate, alpha=0.7)
    plt.plot(cor_submonoequiv_rate, alpha=0.7)
    plt.title('Submonoequivocal rate')
    plt.legend(['Before correction', 'After correction'])
    plt.savefig(os.path.join(outdir, 'submonoequiv_assembly.pdf'),
                format='pdf')
    plt.close()

    high_equiv_rate = np.where(np.logical_and(submonoequiv_rate > 0.9,
                                              coverage > 5))
    logger.info(f'Positions of high submonoequiv rate {high_equiv_rate}')
    high_cor_equiv_rate = np.where(np.logical_and(cor_submonoequiv_rate > 0.9,
                                                  coverage > 5))
    logger.info(f'Positions of high cor submonoequiv rate {high_cor_equiv_rate}')
    logger.info(f'Error rate            [50, 55, 60, 90, 95]%: {np.percentile(error_rate, [50, 55, 60, 90, 95])}')
    logger.info(f'Submonoequiv rate     [50, 55, 60, 90, 95]%: {np.percentile(submonoequiv_rate, [50, 55, 60, 90, 95])}')
    logger.info(f'Cor submonoequiv rate [50, 55, 60, 90, 95]%: {np.percentile(cor_submonoequiv_rate, [50, 55, 60, 90, 95])}')


def iscorrectsubmonomerclosest(mappings,
                               submonoread_set, cor_submonoread_set,
                               submonoassembly,
                               logger):
    is_among_closest = 0
    total_equiv = 0

    attempted_correct = 0
    valid_correct = 0
    impossible_correction = 0

    for r_id, map_pos in mappings.items():
        submonoread = submonoread_set[r_id]
        cor_submonoread = cor_submonoread_set[r_id]
        for i, c in enumerate(submonoread):
            if c == submonoread.submonoequiv_symb:
                sminst = submonoread.submonoinstances[i]
                dist2submonomers = sminst.dist2submonomers
                a_smindex = submonoassembly[map_pos + i]
                total_equiv += 1
                if a_smindex in dist2submonomers.keys():
                    is_among_closest += 1
                    min_dist = min(dist2submonomers.values())

                cor_sminst = cor_submonoread[i]
                if cor_sminst != submonoread.submonoequiv_symb:
                    attempted_correct += 1
                    if a_smindex not in dist2submonomers.keys():
                        impossible_correction += 1
                    if cor_sminst == a_smindex:
                        valid_correct += 1
    logger.info(f'Total equivocal positions {total_equiv}')
    logger.info(f'Assembly submonomer is among the closest {is_among_closest}')
    logger.info(f'Total attempted to correct {attempted_correct}')
    logger.info(f'Valid corrections {valid_correct}')
    logger.info(f'Impossible corrections (assembly submonomer not among closest) {impossible_correction}')
    logger.info(f'among closest / total {is_among_closest / total_equiv}')
    logger.info(f'(attempted_corr-impossible_corr) / among closest {(attempted_correct-impossible_correction) / is_among_closest}')
    logger.info(f'valid correction / attempt correction {valid_correct / attempted_correct}')


def main():
    params = parse_args()
    smart_makedirs(params.outdir)
    logfile = os.path.join(params.outdir, 'benchmark.log')
    logger = \
        get_logger(logfile,
                   logger_name="centroFlye: submonopartition_benchmark")
    submonoread_set, submonoassembly = get_submonostrings(params, logger)
    logger.info(f'Correcting submonostring set')
    cor_submonoread_set = \
        CorrectedSubmonoStringSet.from_submonostring_set(submonoread_set)
    logger.info(f'Mapping reads')
    mappings = map_reads(submonoread_set, submonoassembly)
    run_rate_stats(mappings, submonoassembly, submonoread_set,
                   cor_submonoread_set, params.outdir, logger)
    logger.info('Checking equivocal submonomer instances')
    iscorrectsubmonomerclosest(mappings=mappings,
                               submonoread_set=submonoread_set,
                               submonoassembly=submonoassembly,
                               logger=logger)


if __name__ == "__main__":
    main()
