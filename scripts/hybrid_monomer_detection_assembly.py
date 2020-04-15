import argparse
from joblib import Parallel, delayed
from tqdm import tqdm

from sd_parser import SD_Report
from utils.bio import read_bio_seq, hybrid_alignment, RC, write_bio_seqs

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sd", required=True, help="SD report on assembly")
    parser.add_argument("--monomers", required=True, help="Monomers")
    parser.add_argument("--assembly", required=True, help="Assembly")
    parser.add_argument("--outfile", required=True)
    parser.add_argument("--max-ident-diff", type=int, default=5)
    parser.add_argument("--min-ident", type=float, default=70)
    parser.add_argument("--threads", type=int, default=30)
    parser.add_argument("--score-mult", type=float, default=1.05)
    params = parser.parse_args()
    return params


def read_input(params):
    sd_report = SD_Report(SD_report_fn=params.sd,
                          monomers_fn=params.monomers)
    assembly = read_bio_seq(params.assembly)
    return sd_report, assembly


def get_candidates(sd_report, assembly, min_ident):
    candidates = []
    for p, row in sd_report.df.iterrows():
        sec_identity = row.sec_identity
        if sec_identity < min_ident:
            continue

        st, en = row.r_st, row.r_en + 1
        segm = assembly[st:en]

        m_id, sec_m_id = row.monomer, row.sec_monomer
        a = sd_report.rev_monomer_names_map[m_id]
        b = sd_report.rev_monomer_names_map[sec_m_id]
        A = sd_report.monomers[a]
        B = sd_report.monomers[b]

        if sd_report.is_upper(m_id):
            segm = RC(segm)
        candidates.append((p, segm, a, b, A, B))
    return candidates


def extract_hybrids(candidates, score_mult, threads):
    results = \
        Parallel(n_jobs=threads)(delayed(hybrid_alignment)(segm, A, B)
                                 for _, segm, _, _, A, B in tqdm(candidates))
    hybrids = {}
    for i, (_, _, max_sc, sec_max_sc, jump, orient) in enumerate(results):
        if max_sc >= sec_max_sc * score_mult:
            jump = jump[1:]
            p, _, a, b, A, B = candidates[i]
            if orient == '>':
                hybrid = A[:jump[0]] + B[jump[1]:]
            else:
                hybrid = B[:jump[0]] + A[jump[1]:]
            hybrid_name = f'position={p}|m1={a}|m2={b}|prefix={jump[0]}|suffix={jump[1]}|orient={orient}'
            hybrids[hybrid_name] = hybrid
            print(p, a, b, max_sc, sec_max_sc, jump, orient, hybrid)
    return hybrids


def main():
    params = parse_args()
    sd_report, assembly = read_input(params)
    candidates = get_candidates(sd_report, assembly, params.min_ident)
    hybrids = extract_hybrids(candidates, params.score_mult, params.threads)
    write_bio_seqs(params.outfile, hybrids)


if __name__ == "__main__":
    main()
