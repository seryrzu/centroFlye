import argparse
from collections import Counter
from pprint import pprint

from joblib import Parallel, delayed
import numpy as np
from statsmodels.stats.proportion import proportions_ztest

from sd_parser import SD_Report
from utils.bio import read_bio_seqs, hybrid_alignment, RC, \
                      calc_identity, write_bio_seqs

# TODO: update this script with new API

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sd", required=True, help="SD report")
    parser.add_argument("--monomers", required=True, help="Monomers")
    parser.add_argument("--reads", required=True, help="Cen reads")
    parser.add_argument("--max-ident-diff", type=int, default=5)
    parser.add_argument("--min-sec-ident", type=int, default=80)
    parser.add_argument("--threads", type=int, default=30)
    parser.add_argument("--outfile", required=True)
    parser.add_argument("--genome-len", type=int, required=True)
    parser.add_argument("--isassembly", action='store_true')
    params = parser.parse_args()
    return params


def read_input(params):
    sd_report = SD_Report(SD_report_fn=params.sd,
                          monomers_fn=params.monomers)
    reads = read_bio_seqs(params.reads)
    return sd_report, reads


def get_putative_pairs(names):
    putative_pairs = []
    for m1 in names:
        for m2 in names:
            if m1 < m2:
                putative_pairs.append((m1, m2))
    return putative_pairs


def extract_hybrid(pair, sd_report, reads,
                   A, B,
                   max_ident_diff,
                   min_sec_ident,
                   threads,
                   coverage,
                   sign_lev=0.05,
                   isassembly=False):
    def extract_read_segms(sd_report, df, reads):
        segms = []
        for i, row in df.iterrows():
            r_id, st, en = row.r_id, row.r_st, row.r_en + 1
            monomer = row.monomer
            segm = reads[r_id][st:en]
            if sd_report.is_upper(monomer):
                segm = RC(segm)
            segms.append(segm)
        return segms

    def get_cuts(segms, A, B, t=threads):
        results = \
            Parallel(n_jobs=t)(delayed(hybrid_alignment)(segm, A, B)
                               for segm in segms)
        jumps = []
        for _, _, max_sc, sec_max_sc, jump, orient in results:
            print(max_sc, sec_max_sc, jump, orient)
            if max_sc >= sec_max_sc * 1.05:
                jump = jump[1:]
                jumps.append((jump, orient))
        # print(hybrids)
        return jumps

    def n_improved_pairs(segms, A, B, hybrid):
        n_improved = 0
        min_improved, max_improved = 1, 0
        for segm in segms:
            identA = calc_identity(segm, A)
            identB = calc_identity(segm, B)
            ident_hybrid = calc_identity(segm, hybrid)
            if ident_hybrid > max(identA, identB):
                n_improved += 1
                min_improved = min(min_improved, ident_hybrid)
                max_improved = max(max_improved, ident_hybrid)
        return n_improved, min_improved, max_improved

    results = {
        "status": False
    }
    df = sd_report.df
    bases = list(pair)
    bases += sd_report.get_rev_index(bases)
    AB = np.where(df.monomer.isin(bases) &
                  df.sec_monomer.isin(bases) &
                  (abs(df.identity - df.sec_identity) < max_ident_diff) &
                  (df.sec_identity > min_sec_ident) &
                  df.reliability.isin(['+']))[0]
    dfAB = df.iloc[AB]

    read_segms = extract_read_segms(sd_report, dfAB, reads)
    np.random.shuffle(read_segms)
    if len(read_segms) < coverage * 0.5:
        return results

    import pandas as pd
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(dfAB[['monomer', 'identity', 'sec_monomer', 'sec_identity']].head())

    if isassembly:
        assembly_segms = read_segms
        jumps = Counter(get_cuts(assembly_segms, A, B))
        if len(jumps) == 0:
            return results
        results["jumps"] = jumps
        mc_jump = jumps.most_common(1)[0][0]
        results["mc_jump"] = mc_jump
        mc_jump, orient = mc_jump

        if orient == '>':
            hybrid = A[:mc_jump[0]] + B[mc_jump[1]:]
        else:
            hybrid = B[:mc_jump[0]] + A[mc_jump[1]:]
        results["hybrid"] = hybrid
        results["status"] = True
        return results

    nreads = min(2000, len(read_segms)) // 2
    print(f'nreads = {nreads}')

    results["ntrain"] = nreads
    results["ntest"] = nreads


    train, test = read_segms[:nreads], read_segms[nreads:nreads*2]
    jumps = Counter(get_cuts(train, A, B))
    if len(jumps) == 0:
        return results
    print(jumps)
    results["jumps"] = jumps

    mc_jump = jumps.most_common(1)[0]

    if mc_jump[1] < 5:
        return results
    results["mc_jump"] = mc_jump[0]
    mc_jump, orient = mc_jump[0]

    if orient == '>':
        hybrid = A[:mc_jump[0]] + B[mc_jump[1]:]
    else:
        hybrid = B[:mc_jump[0]] + A[mc_jump[1]:]

    train_n_improved, min_train_ident, max_train_ident = \
        n_improved_pairs(train, A, B, hybrid)
    test_n_improved, min_test_ident, max_test_ident = \
        n_improved_pairs(test, A, B, hybrid)
    min_ident = min(min_train_ident, min_test_ident)
    max_ident = max(max_train_ident, max_test_ident)
    results["min_ident"] = min_ident
    results["max_ident"] = max_ident

    if train_n_improved == 0 or test_n_improved == 0:
        return results
    elif train_n_improved == len(train) or test_n_improved == len(test):
        results["train_n_improved"] = train_n_improved
        results["test_n_improved"] = test_n_improved
        results["hybrid"] = hybrid
        results["status"] = True
        return results

    stat, pval = proportions_ztest([train_n_improved, test_n_improved],
                                   [len(train), len(test)])
    results["stat"] = stat
    results["pval"] = pval
    results["train_n_improved"] = train_n_improved
    results["test_n_improved"] = test_n_improved
    results["hybrid"] = hybrid
    results["status"] = True
    return results


def main():
    params = parse_args()
    sd_report, reads = read_input(params)
    coverage = 171 * len(sd_report.df) / params.genome_len
    print(f'Estimated coverage {coverage}x')

    putative_pairs = get_putative_pairs(sd_report.get_monomer_encoded_names())
    # putative_pairs = [(1, 3)]
    # print(sd_report.monomer_names_map['S2C2H1L.4'], sd_report.monomer_names_map['S2C2H1L.2'])

    hybrids = {}
    for pair in putative_pairs:
        a, b = pair
        print(pair)
        print(sd_report.rev_monomer_names_map[a], sd_report.rev_monomer_names_map[b])

        A = sd_report.monomers[sd_report.rev_monomer_names_map[a]]
        B = sd_report.monomers[sd_report.rev_monomer_names_map[b]]
        extraction_res = extract_hybrid(pair, sd_report, reads,
                                        A, B,
                                        params.max_ident_diff,
                                        params.min_sec_ident,
                                        params.threads,
                                        coverage,
                                        isassembly=params.isassembly)
        if extraction_res["status"]:
            pprint(extraction_res, width=1)
            print("")
            hybrid, (jump, orient) = extraction_res["hybrid"], extraction_res["mc_jump"]
            hybrids[f'{a}_{b}|{jump[0]}_{jump[1]}|{orient}'] = hybrid
    write_bio_seqs(params.outfile, hybrids)


if __name__ == "__main__":
    main()
