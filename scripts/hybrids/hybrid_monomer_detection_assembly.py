# (c) 2020 by Authors
# This file is a part of centroFlye program.
# Released under the BSD license (see LICENSE file)


import argparse
from collections import defaultdict
import os
import sys

this_dirname = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(this_dirname, os.path.pardir))

from joblib import Parallel, delayed
import networkx as nx
from tqdm import tqdm

from sd_parser.sd_parser import SD_Report
from standard_logger import get_logger
from utils.bio import hybrid_alignment, write_bio_seqs, calc_identity
from utils.git import get_git_revision_short_hash
from utils.os_utils import smart_makedirs
from utils.various import list2str, fst_iterable, running_mean


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sd", required=True, help="SD report on assembly")
    parser.add_argument("--monomers", required=True, help="Monomers")
    parser.add_argument("--assembly", required=True, help="Assembly")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--max-ident-diff", type=int, default=.04)
    parser.add_argument("--min-sec-ident", type=float, default=.75)
    parser.add_argument("--threads", type=int, default=30)
    parser.add_argument("--score-mult", type=float, default=1.1)
    parser.add_argument("--max-sim", type=float, default=.95)
    parser.add_argument("--ident-mov-av-len", type=int, default=30)
    parser.add_argument("--min-moving-identity", type=int, default=.95)
    params = parser.parse_args()
    return params


def get_candidates(monoassembly,
                   min_sec_ident, max_ident_diff,
                   ident_mov_av_len,
                   min_moving_identity,
                   logger):
    candidates = []
    logger.info('Candidates:')
    identities = monoassembly.get_identities()
    identities_running_mean = \
        running_mean(identities, ident_mov_av_len)
    for pos, ident_running_mean in enumerate(identities_running_mean):
        if ident_running_mean < min_moving_identity:
            continue
        mi = monoassembly.monoinstances[pos + ident_mov_av_len//2]
        ident_diff = abs(mi.identity - mi.sec_identity)
        if mi.is_reliable() and \
                mi.sec_identity >= min_sec_ident and \
                ident_diff <= max_ident_diff:
            a = mi.get_monoid()
            b = mi.get_secmonoid()
            mono_ind1, mono_ind2 = mi.get_monoindex(), mi.get_secmonoindex()
            if mono_ind1 == mono_ind2:
                continue

            A = monoassembly.monomer_db.get_seqs_by_index(mi.get_monoindex())
            B = monoassembly.monomer_db.get_seqs_by_index(mi.get_secmonoindex())
            A, B = list(A), list(B)
            if len(A) > 1 or len(B) > 1:
                continue
            # only work with clusters of size 1
            A, B = A[0], B[0]

            candidates.append((pos, mi.nucl_segment, a, b, A, B))
            logger.info(f'pos={pos} ident={mi.identity}, '
                        f'sec_ident={mi.sec_identity} {a} {b}')
    return candidates


def filter_hybrid_alignments(hybrid_aligns, candidates,
                             score_mult, max_sim, logger):
    hybrids2info, info2hybrids = {}, {}
    hybrids2pos = defaultdict(list)
    pair2jumps = defaultdict(set)
    for i, (_, _, max_sc, sec_sc, jump, orient) in enumerate(hybrid_aligns):
        if max_sc >= sec_sc * score_mult:
            jump = jump[1:]
            pos, _, a, b, A, B = candidates[i]
            if orient == '<':
                a, b = b, a
                A, B = B, A
            hybrid = A[:jump[0]] + B[jump[1]:]
            ident_AB = calc_identity(A, B)
            ident_HA = calc_identity(hybrid, A)
            ident_HB = calc_identity(hybrid, B)
            max_ident = max(ident_AB, ident_HA, ident_HB)
            if max_ident <= max_sim:
                logger.info(f'pos={pos} a={a} b={b}')
                logger.info(f'pref={jump[0]} suff={jump[1]}')
                logger.info(f'ident(a, b)={ident_AB:0.2}')
                logger.info(f'ident(hybrid, a)={ident_HA:0.2}')
                logger.info(f'ident(hybrid, b)={ident_HB:0.2}')
                logger.info(f'maxscore={max_sc} secmaxscore = {sec_sc}')
                logger.info('---')
                info = (a, b, jump)
                hybrids2info[hybrid] = info
                info2hybrids[info] = hybrid
                pair2jumps[(a, b)].add(jump)
                hybrids2pos[hybrid].append(pos)

    main_hybrids2pos, main_hybrids2info = {}, {}
    for pair, jumps in pair2jumps.items():
        pair_hybrids_mult = {}
        all_pos = []
        for jump in jumps:
            info = (*pair, jump)
            hybrid = info2hybrids[info]
            positions = hybrids2pos[hybrid]
            pair_hybrids_mult[hybrid] = len(positions)
            all_pos += positions
        all_pos.sort()

        pop_hybrid = max(pair_hybrids_mult, key=pair_hybrids_mult.get)
        _, _, pop_jump = hybrids2info[pop_hybrid]
        logger.info(f'For pair {pair} selected jump: {pop_jump} '
                    f'all jumps: {jumps}')
        for jump in jumps:
            info = (*pair, jump)
            hybrid = info2hybrids[info]
            if jump == pop_jump:
                main_hybrids2pos[hybrid] = all_pos
                main_hybrids2info[hybrid] = info
    hybrids2pos, hybrids2info = main_hybrids2pos, main_hybrids2info

    logger.info('---')
    logger.info('Summary of extracted hybrids')
    for hybrid, info in hybrids2info.items():
        logger.info(f'info = {info}')
        logger.info(f'pos = {hybrids2pos[hybrid]}')
        logger.info(hybrid)
        logger.info('---')
    return hybrids2info, hybrids2pos


def cluster_hybrids(hybrids2info, hybrids2pos, max_sim, logger):
    logger.info('Constructing an identity graph')
    ident_graph = nx.Graph()
    for h1, info1 in hybrids2info.items():
        ident_graph.add_node(h1)
        for h2, info2 in hybrids2info.items():
            if h1 <= h2:
                continue
            ident = calc_identity(h1, h2)
            if ident >= max_sim:
                ident_graph.add_edge(h1, h2)
                logger.info(f'Ident {ident} b/w {info1} {info2}')

    logger.info('Finished constructing an identity graph')

    logger.info('Selecting of most popular hybrids in each component')
    filt_hybrids2info = {}
    filt_hybrids2pos = defaultdict(list)
    for i, cc in enumerate(nx.connected_components(ident_graph)):
        cc_size = len(cc)
        logger.info(f'Analysis of cc #{i} of size {cc_size}')
        if cc_size == 1:
            hybrid = fst_iterable(cc)
            filt_hybrids2pos[hybrid] = hybrids2pos[hybrid]
            filt_hybrids2info[hybrid] = hybrids2info[hybrid]
            logger.info(f'hybrid info = {hybrids2info[hybrid]}')
            logger.info(f'       pos  = {hybrids2pos[hybrid]}')
            continue

        hybrids_mult = {hybrid: len(hybrids2pos[hybrid]) for hybrid in cc}
        pop_hybrid = max(hybrids_mult, key=hybrids_mult.get)
        filt_hybrids2info[pop_hybrid] = hybrids2info[pop_hybrid]
        logger.info(f'Most popular hybrid info = {hybrids2info[pop_hybrid]}')
        logger.info(f'  sequence = {pop_hybrid}')
        logger.info('List of info and pos for all hybrids in cc')
        for hybrid in cc:
            logger.info(f'  info = {hybrids2info[hybrid]}')
            logger.info(f'  pos  = {hybrids2pos[hybrid]}')
            filt_hybrids2pos[pop_hybrid] += hybrids2pos[hybrid]
        filt_hybrids2pos[pop_hybrid].sort()
        logger.info(f'All positions for cc: {filt_hybrids2pos[pop_hybrid]}')
        logger.info('---')
    return filt_hybrids2pos, filt_hybrids2info


def extract_hybrids(monoassembly,
                    min_sec_ident, max_sim, max_ident_diff, score_mult,
                    ident_mov_av_len,
                    min_moving_identity,
                    threads, logger):
    candidates = get_candidates(monoassembly, min_sec_ident, max_ident_diff,
                                ident_mov_av_len,
                                min_moving_identity,
                                logger)
    hybrid_aligns = \
        Parallel(n_jobs=threads)(delayed(hybrid_alignment)(segm, A, B)
                                 for _, segm, _, _, A, B in tqdm(candidates))
    hybrids2info, hybrids2pos = \
        filter_hybrid_alignments(hybrid_aligns, candidates,
                                 score_mult, max_sim, logger)

    filt_hybrids2pos, filt_hybrids2info = \
        cluster_hybrids(hybrids2info, hybrids2pos, max_sim, logger)
    return filt_hybrids2pos, filt_hybrids2info


def export_hybrids(hybrids2pos, hybrids2info, outdir, monomer_db, sep='\t'):
    hybrids = {}
    tsv_outfile = os.path.join(outdir, 'hybrids.tsv')
    with open(tsv_outfile, 'w') as f:
        header = ['monomer1', 'monomer2', 'jump1', 'jump2', 'pos', 'seq']
        print(list2str(header), file=f)
        for hybrid, pos in hybrids2pos.items():
            a, b, (jump_fst, jump_snd) = hybrids2info[hybrid]
            outlist = [a, b, jump_fst, jump_snd, pos, hybrid]
            outline = list2str(outlist, sep=sep)
            print(outline, file=f)
            monoindex1 = monomer_db.id2index[a]
            monoindex2 = monomer_db.id2index[b]
            hybrid_id = f'{monoindex1}_{monoindex2}_{jump_fst}_{jump_snd}'
            hybrids[hybrid_id] = hybrid

    hybrids_outfile = os.path.join(outdir, 'hybrids.fasta')
    write_bio_seqs(hybrids_outfile, hybrids)

    combined = monomer_db.get_monomers_dict()
    combined.update(hybrids)
    combined_outfile = os.path.join(outdir, 'monomers_hybrids.fasta')
    write_bio_seqs(combined_outfile, combined)


def main():
    params = parse_args()
    smart_makedirs(params.outdir)
    logfn = os.path.join(params.outdir, 'hybrid_detection_assembly.log')
    logger = get_logger(logfn,
                        logger_name='centroFlye: hybrid_detection_assembly')

    logger.info(f'cmd: {sys.argv}')
    logger.info(f'git hash: {get_git_revision_short_hash()}')

    logger.info('Reading SD reports')
    sd_report = SD_Report(sd_report_fn=params.sd,
                          monomers_fn=params.monomers,
                          sequences_fn=params.assembly,
                          mode='assembly')
    logger.info('Finished reading SD reports')

    monomer_db = sd_report.monomer_db
    monoassembly = fst_iterable(sd_report.monostring_set.monostrings.values())
    hybrids2pos, hybrids2info = \
        extract_hybrids(monoassembly=monoassembly,
                        min_sec_ident=params.min_sec_ident,
                        max_sim=params.max_sim,
                        max_ident_diff=params.max_ident_diff,
                        score_mult=params.score_mult,
                        ident_mov_av_len=params.ident_mov_av_len,
                        min_moving_identity=params.min_moving_identity,
                        threads=params.threads,
                        logger=logger)
    export_hybrids(hybrids2pos, hybrids2info, params.outdir, monomer_db)


if __name__ == "__main__":
    main()
