# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)


import argparse
import os
import sys

this_dirname = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(this_dirname, os.path.pardir))

import edlib

from sd_parser.sd_parser import SD_Report
from standard_logger import get_logger
from submonomers.submonomer_db import SubmonomerDB
from utils.os_utils import smart_makedirs
from utils.bio import perfect_overlap


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


def get_coverage(monoassembly, monoreads_set):
    assembly_len = len(monoassembly.raw_monostring)
    total_reads_len = sum(len(monoread.raw_monostring)
                          for monoread in monoreads_set.monostrings.values())
    return total_reads_len / assembly_len


def get_closest_submonomer(absent, submonomers, max_dist=5):
    best_dist = 1000
    best_submonomers = []
    for submonomer in submonomers:
        alignment = edlib.align(absent, submonomer, k=20)
        dist = alignment['editDistance']
        if dist == -1:
            continue
        if dist < best_dist:
            best_submonomers = [submonomer]
            best_dist = dist
        elif dist == best_dist:
            best_submonomers.append(submonomer)
    perfect_overlap_monomers = []
    for submonomer in best_submonomers:
        overlap = perfect_overlap(absent, submonomer)
        if len(overlap) >= min(len(absent), len(submonomer)) - max_dist:
            perfect_overlap_monomers.append(submonomer)
    return best_dist, perfect_overlap_monomers


def benchmark_read_submonomers_recruitment(assembly_submonomers, reads_submonomers):
    stats = {}
    inters = assembly_submonomers & reads_submonomers
    absent = assembly_submonomers - reads_submonomers
    stats['#subm_assembly'] = len(assembly_submonomers)
    stats['#subm_reads'] = len(reads_submonomers)
    stats['#overlap'] = len(inters)

    perfect_overlap_cnt = 0
    for monomer in absent:
        dist, perf_overlap_submonomers = \
            get_closest_submonomer(monomer, reads_submonomers)
        if len(perf_overlap_submonomers) == 1:
            perfect_overlap_cnt += 1
    stats['perfect overlap'] = perfect_overlap_cnt
    total_missing = len(assembly_submonomers) - len(inters) - perfect_overlap_cnt
    stats['total_missing'] = total_missing
    total_nonassembly = len(reads_submonomers) - len(inters) - perfect_overlap_cnt
    stats['total_nonassembly'] = total_nonassembly
    return stats


def print_stats(stats, logger):
    logger.info(f"# assembly_submonomers = {stats['#subm_assembly']}")
    logger.info(f"# reads_submonomers = {stats['#subm_reads']}")
    logger.info(f"# overlap = {stats['#overlap']}")
    logger.info(f"perfect_overlap_cnt = {stats['perfect overlap']}")
    logger.info(f"in assembly but missing in reads = {stats['total_missing']}")
    logger.info(f"in reads but missing in assembly = {stats['total_nonassembly']}")


def benchmark(submonomers_reads, submonomers_assembly, logger):
    n_monomers = submonomers_reads.monomer_db.get_size()
    assert n_monomers == submonomers_assembly.monomer_db.get_size()
    all_stats = []
    for monoindex in range(n_monomers):
        sm_reads_index = submonomers_reads.get_submonomers_by_mono_index(monoindex)
        sm_assembly_index = submonomers_assembly.get_submonomers_by_mono_index(monoindex)
        sm_reads_index_seqs = set(sm.seq for sm in sm_reads_index)
        sm_assembly_index_seqs = set(sm.seq for sm in sm_assembly_index)
        stats = benchmark_read_submonomers_recruitment(sm_assembly_index_seqs, sm_reads_index_seqs)
        logger.info(f'Statistics for monoindex {monoindex}')
        print_stats(stats, logger)
        all_stats.append(stats)
    total_stats = {k: sum(stats[k] for stats in all_stats)
                   for k in all_stats[0].keys()}
    logger.info(f'Total statistics')
    print_stats(total_stats, logger)


def main():
    params = parse_args()
    smart_makedirs(params.outdir)
    logfile = os.path.join(params.outdir, 'benchmark.log')
    logger = \
        get_logger(logfile,
                   logger_name="centroFlye: submonomers_extraction_benchmark")
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

    logger.info(f'Extracting submonomers from the reads')
    submonomers_reads = SubmonomerDB.from_monostring_set(monoreads_set,
                                                         coverage=coverage)
    submonomers_reads_fn = os.path.join(params.outdir,
                                        'submonomer_db_reads.fasta)
    submonomers_reads.to_fasta(submonomers_reads_fn)

    logger.info(f'Extracting submonomers from the assembly')
    submonomers_assembly = \
        SubmonomerDB.from_monostring_set(sd_report_assembly.monostring_set,
                                         coverage=1)
    submonomers_assembly_fn = os.path.join(params.outdir,
                                           'submonomer_db_reads.fasta)
    submonomers_assembly.to_fasta(submonomers_assembly_fn)

    benchmark(submonomers_reads=submonomers_reads,
              submonomers_assembly=submonomers_assembly,
              logger=logger)


if __name__ == "__main__":
    main()
