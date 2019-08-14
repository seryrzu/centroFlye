import argparse
import os
import json
import numpy as np
from collections import Counter, defaultdict
from utils.os_utils import smart_makedirs
from utils.bio import gen_random_seq, read_bio_seq, write_bio_seqs
from utils.json_utils import stringify_keys


def generate_mutations(unit, mult, div_rate, flank_len=200000):
    n_mut = np.random.binomial(n=len(unit)*mult, p=div_rate, size=1)[0]
    unit_ns = np.random.choice(range(mult), size=n_mut)
    unit_ns = Counter(unit_ns)

    units = defaultdict(lambda: unit)
    all_muts = defaultdict(list)
    for unit_n, n in unit_ns.items():
        muts_pos = np.random.choice(range(len(unit)), size=n, replace=False)
        new_unit = list(unit)
        for pos in muts_pos:
            bases = list("ACGT")
            bases.remove(new_unit[pos])
            new_unit[pos] = np.random.choice(bases)[0]
            all_muts[unit_n].append((int(pos), new_unit[pos]))
        units[unit_n] = ''.join(new_unit)
    tr = ''.join(units[i] for i in range(mult))
    left_flanked_tr = gen_random_seq(flank_len) + tr
    flanked_tr = left_flanked_tr + gen_random_seq(flank_len)
    return tr, left_flanked_tr, flanked_tr, all_muts


def output_results(tr, left_flanked_tr, flanked_tr, all_muts, output_dir):
    smart_makedirs(output_dir)
    write_bio_seqs(os.path.join(output_dir, 'tandem_repeat.fasta'),
                   {'sim_tr': tr})
    write_bio_seqs(os.path.join(output_dir,
                                'left_flanked_tandem_repeat.fasta'),
                   {'left_flanked_sim_tr': left_flanked_tr})
    write_bio_seqs(os.path.join(output_dir, 'flanked_tandem_repeat.fasta'),
                   {'flanked_sim_tr': flanked_tr})
    with open(os.path.join(output_dir, 'all_muts.json'), 'w') as f:
        all_muts = dict(all_muts)
        all_muts = stringify_keys(all_muts)
        print(json.dumps(all_muts), file=f)
    with open(os.path.join(output_dir, 'simulation.log'), 'w') as f:
        total_n_mut = sum(len(x) for x in all_muts.values())
        print(f'full_tr_len = {len(tr)}', file=f)
        print(f'total_n_mut = {total_n_mut}', file=f)
        for pos, muts in all_muts.items():
            print(f'{pos} : {len(muts)}', file=f)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--unit",
                        help="Unit of tandem repeat. Default: random")
    parser.add_argument("--unit-len",
                        help="Unit length used in case no unit is provided",
                        type=int,
                        default=200)
    parser.add_argument("--multiplicity",
                        help="Multiplicity of the repeat to generate",
                        required=True,
                        type=int)
    parser.add_argument("--div-rate",
                        help="Average divergence rate between blocks",
                        type=float,
                        required=True)
    parser.add_argument("-o", "--output",
                        help="Output directory",
                        required=True)
    parser.add_argument("--seed", help="Seed", type=int)
    params = parser.parse_args()

    if params.seed is not None:
        np.random.seed(params.seed)

    if params.unit is None:
        unit = gen_random_seq(length=params.unit_len)
    else:
        unit = read_bio_seq(params.unit)

    tr, left_flanked_tr, flanked_tr, all_muts = generate_mutations(unit, params.multiplicity, params.div_rate)
    output_results(tr, left_flanked_tr, flanked_tr, all_muts, params.output)


if __name__ == "__main__":
    main()
