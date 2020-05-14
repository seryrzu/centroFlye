# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)


import argparse
from collections import Counter
import os
import sys

from joblib import Parallel, delayed
import numpy as np
import networkx as nx

import matplotlib.pyplot as plt

this_dirname = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(this_dirname, os.path.pardir))

from debruijn_graph import iterative_graph, graph3col
from sd_parser.sd_parser import SD_Report
from standard_logger import get_logger
from submonomers.submonomer_db import SubmonomerDB
from submonomers.submonostring import SubmonoString
from submonomers.submonomers_extraction_benchmark import get_coverage
from submonomers.submonostring_set import SubmonoStringSet,\
                                          CorrectedSubmonoStringSet
from utils.os_utils import smart_makedirs, expandpath


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sd-report-reads", required=True)
    parser.add_argument("--sd-report-assembly", required=True)
    parser.add_argument("--monomers", required=True)
    parser.add_argument("--reads", required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--outdir", required=True)
    params = parser.parse_args()
    params.sd_report_reads = expandpath(params.sd_report_reads)
    params.sd_report_assembly = expandpath(params.sd_report_assembly)
    params.monomers = expandpath(params.monomers)
    params.reads = expandpath(params.reads)
    params.assembly = expandpath(params.assembly)
    params.outdir = expandpath(params.outdir)
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


def estimate_rates(mappings, submonoassembly, submonoread_set,
                   min_ident=0.3, min_cov=5):
    assembly_len = len(submonoassembly)
    errors_asm = np.zeros(assembly_len)
    submonoequiv_asm = np.zeros(assembly_len)
    coverage_asm = np.zeros(assembly_len)

    errors_subm = Counter()
    submonoequiv_subm = Counter()
    coverage_subm = Counter()

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
            coverage_asm[ia] += 1
            coverage_subm[ca] += 1
            if cr != ca and cr not in gap_symbols:
                errors_asm[ia] += 1
                errors_subm[ca] += 1
            if cr == submonoread.submonoequiv_symb:
                submonoequiv_asm[ia] += 1
                submonoequiv_subm[ca] += 1

    def norm_numpy(array, cov, min_cov=min_cov):
        rate = array / cov
        rate[cov < min_cov]=0
        return rate

    error_rate_asm = norm_numpy(errors_asm, coverage_asm)
    submonoequiv_rate_asm = norm_numpy(submonoequiv_asm, coverage_asm)

    def norm_dict(vals, cov, min_cov=min_cov):
        return {k: vals[k] / cov[k] for k in vals.keys()
                if cov[k] >= min_cov}

    error_rate_subm = norm_dict(errors_subm, coverage_subm)
    submonoequiv_rate_subm = norm_dict(submonoequiv_subm, coverage_subm)

    return coverage_asm,  error_rate_asm,  submonoequiv_rate_asm, \
           coverage_subm, error_rate_subm, submonoequiv_rate_subm


def run_rate_stats(mappings, submonoassembly, submonoread_set,
                   cor_submonoread_set, outdir, logger):
    coverage_asm, error_rate_asm, submonoequiv_rate_asm, \
        coverage_subm, error_rate_subm, submonoequiv_rate_subm = \
        estimate_rates(mappings,
                       submonoassembly,
                       submonoread_set)
    plt.plot(coverage_asm)
    plt.title('Coverage with uniquely mapped reads')
    plt.savefig(os.path.join(outdir, 'coverage_asm.pdf'), format='pdf')
    plt.close()

    _, cor_error_rate_asm, cor_submonoequiv_rate_asm, \
        cor_coverage_subm, cor_error_rate_subm, cor_submonoequiv_rate_subm = \
        estimate_rates(mappings,
                       submonoassembly,
                       cor_submonoread_set)

    plt.plot(error_rate_asm, alpha=0.7)
    plt.plot(cor_error_rate_asm, alpha=0.7)
    plt.title('Error rate')
    plt.legend(['Before correction', 'After correction'])
    plt.savefig(os.path.join(outdir, 'error_rate_assembly.pdf'), format='pdf')
    plt.close()

    plt.plot(submonoequiv_rate_asm, alpha=0.7)
    plt.plot(cor_submonoequiv_rate_asm, alpha=0.7)
    plt.title('Submonoequivocal rate')
    plt.legend(['Before correction', 'After correction'])
    plt.savefig(os.path.join(outdir, 'submonoequiv_assembly.pdf'),
                format='pdf')
    plt.close()

    plt.hist(error_rate_subm.values(), bins=100, alpha=0.7)
    plt.hist(cor_error_rate_subm.values(), bins=100, alpha=0.7)
    plt.yscale('log')
    plt.title('Histogram of Error Rate per monomer')
    plt.xlabel('Error rate')
    plt.ylabel('Count (log)')
    plt.legend(['Before correction', 'After correction'])
    plt.savefig(os.path.join(outdir, 'error_rate_per_submonomer.pdf'),
                format='pdf')
    plt.close()


    plt.hist(submonoequiv_rate_subm.values(), bins=100, alpha=0.7)
    plt.hist(cor_submonoequiv_rate_subm.values(), bins=100, alpha=0.7)
    plt.yscale('log')
    plt.title('Histogram of Submonoequivocal Rate per monomer')
    plt.xlabel('Submonoequivocal rate')
    plt.ylabel('Count (log)')
    plt.legend(['Before correction', 'After correction'])
    plt.savefig(os.path.join(outdir, 'submonoequiv_rate_per_submonomer.pdf'),
                format='pdf')
    plt.close()

    zero_cov = np.where(coverage_asm == 0)[0]
    logger.info(f'Positions of zero coverage {zero_cov}')

    high_error_rate = np.where(np.logical_and(error_rate_asm > 0.9, coverage_asm > 5))

    logger.info(f'Positions of high error rate {high_error_rate}')

    high_equiv_rate = np.where(np.logical_and(submonoequiv_rate_asm > 0.9,
                                              coverage_asm > 5))
    logger.info(f'Positions of high submonoequiv rate {high_equiv_rate}')
    high_cor_equiv_rate = np.where(np.logical_and(cor_submonoequiv_rate_asm > 0.9,
                                                  coverage_asm > 5))
    logger.info(f'Positions of high cor submonoequiv rate {high_cor_equiv_rate}')
    logger.info(f'Error rate            [50, 55, 60, 90, 95]%: {np.percentile(error_rate_asm, [50, 55, 60, 90, 95])}')
    logger.info(f'Cor error rate        [50, 55, 60, 90, 95]%: {np.percentile(cor_error_rate_asm, [50, 55, 60, 90, 95])}')
    logger.info(f'Submonoequiv rate     [50, 55, 60, 90, 95]%: {np.percentile(submonoequiv_rate_asm, [50, 55, 60, 90, 95])}')
    logger.info(f'Cor submonoequiv rate [50, 55, 60, 90, 95]%: {np.percentile(cor_submonoequiv_rate_asm, [50, 55, 60, 90, 95])}')


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


def three_col_graph_stats(submonoassembly, cor_submonoread_set, outdir,
                          mink=15, maxk=20, def_min_mult=2):
    # TODO refactor after the graph refactoring is complete
    _, dbs_assembly, uncompr_dbs_assembly, _, _ = \
        iterative_graph({'assembly': submonoassembly},
                        min_k=maxk, max_k=maxk,
                        outdir=os.path.join(outdir, 'assembly_graph'),
                        def_min_mult=1)

    _, dbs, uncompr_dbs, _, _ = \
        iterative_graph(cor_submonoread_set.cor_submonostrings,
                        min_k=mink, max_k=maxk,
                        outdir=os.path.join(outdir, 'cor_read_graph'),
                        def_min_mult=def_min_mult)
    gr3col = graph3col(uncompr_dbs_assembly[maxk], uncompr_dbs[maxk])
    gr3col_outdir = os.path.join(outdir, 'gr3col')
    smart_makedirs(gr3col_outdir)
    dot_file = os.path.join(gr3col_outdir, f'db_k{maxk}.dot')
    nx.drawing.nx_pydot.write_dot(gr3col, dot_file)


def main():
    params = parse_args()
    smart_makedirs(params.outdir)
    logfile = os.path.join(params.outdir, 'benchmark.log')
    logger = \
        get_logger(logfile,
                   logger_name="centroFlye: submonopartition_benchmark")
    submonoread_set, submonoassembly = get_submonostrings(params, logger)
    submonomer_db_fn = os.path.join(params.outdir,
                                    'submonomer_db.fasta')
    logger.info(f'Exporting submonomer_db to {submonomer_db_fn}')
    submonoread_set.submonomer_db.to_fasta(submonomer_db_fn)

    submonoassembly_fn = os.path.join(params.outdir, 'submonoassembly.tsv')
    logger.info(f'Exporting submonoassembly to {submonoassembly_fn}')
    submonoassembly.to_tsv(submonoassembly_fn)

    submonoread_set_fn = os.path.join(params.outdir, 'submonoread_set.tsv')
    logger.info(f'Exporting submonoreads to {submonoread_set_fn}')
    submonoread_set.to_tsv(submonoread_set_fn)

    logger.info(f'Correcting submonostring set')
    cor_submonoread_set = \
        CorrectedSubmonoStringSet.from_submonostring_set(submonoread_set)

    cor_submonoread_set_fn = os.path.join(params.outdir,
                                          'cor_submonoread_set.tsv')
    logger.info(f'Exporting corrected submonoreads to {cor_submonoread_set_fn}')
    cor_submonoread_set.to_tsv(cor_submonoread_set_fn)

    logger.info(f'Mapping reads')
    mappings = map_reads(submonoread_set, submonoassembly)
    run_rate_stats(mappings, submonoassembly, submonoread_set,
                   cor_submonoread_set, params.outdir, logger)
    logger.info('Checking equivocal submonomer instances')
    iscorrectsubmonomerclosest(mappings=mappings,
                               submonoread_set=submonoread_set,
                               cor_submonoread_set=cor_submonoread_set,
                               submonoassembly=submonoassembly,
                               logger=logger)

    logger.info('Building graphs')
    three_col_graph_stats(submonoassembly, cor_submonoread_set, params.outdir,
                          mink=15, maxk=20, def_min_mult=2)


if __name__ == "__main__":
    main()
