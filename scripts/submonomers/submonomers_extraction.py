# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)

import argparse
from collections import Counter
import logging
import os
import sys

this_dirname = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(this_dirname, os.path.pardir))

import edlib
import networkx as nx

from sd_parser.sd_parser import SD_Report
from standard_logger import get_logger
import submonomers.submonomer_db
from utils.bio import perfect_overlap, long_substr
from utils.os_utils import smart_makedirs, expandpath

logger = logging.getLogger("centroFlye.submonomers.submonomer_extraction")


def _reduce2overlaps(nucl_segms, max_len_diff=15, max_overlap_cut=15):
    overlap_graph = nx.Graph()
    overlaps = []
    for s1 in nucl_segms:
        for s2 in nucl_segms:
            if s2 >= s1:
                continue
            if abs(len(s1) - len(s2)) > max_len_diff:
                continue
            alignment = edlib.align(s1, s2, k=max_overlap_cut)
            if alignment['editDistance'] != -1:
                overlap = perfect_overlap(s1, s2)
                if len(overlap) >= max(len(s1), len(s2)) - max_overlap_cut:
                    overlap_graph.add_edge(s1, s2)
                    overlaps.append(len(overlap))

    reduced2overlaps = {s: s for s in nucl_segms}

    for cc in nx.connected_components(overlap_graph):
        cc = list(cc)
        overlap = long_substr(cc)
        for s in cc:
            reduced2overlaps[s] = overlap

    return reduced2overlaps


def _reduce2shortest(reduced2overlaps):
    seq2shortest = {s: s for s in reduced2overlaps.values()}
    for s1 in reduced2overlaps.values():
        for s2 in reduced2overlaps.values():
            if len(s2) >= len(seq2shortest[s1]):
                continue
            if s2 in s1:
                seq2shortest[s1] = s2
    reduced2shortest = {}
    for s, cs in reduced2overlaps.items():
        reduced2shortest[s] = seq2shortest[cs]
    return reduced2shortest


def _submonomer_seqs_for_monomer_extraction(nucl_segms, coverage,
                                            coverage_coef=0.25):
    cnt_segments = Counter(nucl_segms)
    freq_segments = set(seq for seq, cnt in cnt_segments.items()
                        if cnt >= coverage*coverage_coef)
    reduced2overlaps = _reduce2overlaps(freq_segments)
    reduced2shortest = _reduce2shortest(reduced2overlaps)
    submonomers_seqs = reduced2shortest.values()
    submonomers_seqs = list(set(submonomers_seqs))
    return submonomers_seqs


def _submonomer_seqs2submonomers(submonomer_seqs, st_submono_index, monomer):
    submonomer_list = []
    for i, seq in enumerate(submonomer_seqs, start=st_submono_index):
        submonomer = submonomers.submonomer_db.Submonomer(submono_index=i,
                                                          monomer=monomer,
                                                          seq=seq)
        submonomer_list.append(submonomer)
    return submonomer_list


def submonomer_db_extraction(monostring_set, coverage):
    monomer_db = monostring_set.monomer_db
    monoindexes = monomer_db.get_monoindexes()
    monomerinstances_dict = monostring_set.classify_monomerinstances()

    all_submonomers = []
    st_submono_index = 0
    for monoindex in monoindexes:
        logger.debug(f'Extracting submonomers for monoindex={monoindex}')
        mis = monomerinstances_dict[monoindex]
        nucl_segms = [mi.nucl_segment for mi in mis]
        submonomer_seqs = \
            _submonomer_seqs_for_monomer_extraction(nucl_segms, coverage)
        logger.debug(f'    Extracted {len(submonomer_seqs)} submonomers')

        monomer = monomer_db.monomers[monoindex]
        submonomers_list = _submonomer_seqs2submonomers(submonomer_seqs,
                                                        st_submono_index,
                                                        monomer)
        st_submono_index += len(submonomers_list)
        all_submonomers += submonomers_list
    submonomer_db = \
        submonomers.submonomer_db.SubmonomerDB(submonomers=all_submonomers,
                                               monomer_db=monomer_db)
    return submonomer_db


def _parse_args():
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
    params = _parse_args()
    smart_makedirs(params.outdir)

    logfile = os.path.join(params.outdir, 'extraction.log')
    logger = get_logger(logfile,
                        logger_name="centroFlye: submonomers_extraction")

    logger.info('Reading SD Report')
    logger.info(f'    sd_report_fn = {params.sd_report}')
    logger.info(f'    monomers_fn  = {params.monomers}')
    logger.info(f'    sequences_fn = {params.sequences}')
    sd_report = SD_Report(sd_report_fn=params.sd_report,
                          sequences_fn=params.sequences,
                          monomers_fn=params.monomers)
    logger.info('Finished reading SD report')

    logger.info(f'Extracting submonomer db with coverage = {params.coverage}')
    submonomer_db = \
        submonomers.submonomer_db.SubmonomerDB.from_monostring_set(
            sd_report.monostring_set,
            coverage=params.coverage)
    submonomer_db_fn = os.path.join(params.outdir, 'submonomer_db.fasta')
    logger.info('Finished extracting')
    logger.info(f'Exporing to {submonomer_db_fn}')
    submonomer_db.to_fasta(filename=submonomer_db_fn)


if __name__ == "__main__":
    main()
