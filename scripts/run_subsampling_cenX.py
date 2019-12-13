import argparse
import os
import random
import subprocess

from utils.bio import read_bio_seqs, write_bio_seqs
from utils.os_utils import smart_makedirs
from utils.trim_seqs import trim_seqs


seeds = [12398182, 3812983, 1723992]
preads = [0.95, 0.9, 0.8]
cut_freqs = [0, 0.1, 0.2]


SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
unit_fn = f'{SCRIPT_DIR}/../supplementary_data/DXZ1_rc.fasta'


def run(all_reads, seed, cut_freq, pread, outdir):
    print('New run')
    print(seed, cut_freq, pread)
    nreads = int(len(all_reads)*pread)
    read_ids = random.sample(all_reads.keys(), nreads)
    reads = {r_id: all_reads[r_id] for r_id in read_ids}
    reads = trim_seqs(reads, cut_freq)
    read_reduction = sum(len(v) for v in reads.values()) / \
                     sum(len(v) for v in all_reads.values())
    print('read_reduction', read_reduction)
    reads_fn = os.path.join(outdir, 'centromeric_reads.fasta')
    print(reads_fn)
    write_bio_seqs(reads_fn, reads)

    # Run NCRF
    ncrf_outdir = os.path.join(outdir, 'NCRF_rc')
    ncrf_cmd = ['python', '-u', f'{SCRIPT_DIR}/run_ncrf_parallel.py',
                '--reads', reads_fn,
                '-t', '50',
                '--outdir', ncrf_outdir,
                '--repeat', unit_fn]
    print(' '.join(ncrf_cmd))
    subprocess.call(ncrf_cmd)
    ncrf_fn = os.path.join(ncrf_outdir, 'report.ncrf')

    # Run k-mer recruitment
    recr_outdir = os.path.join(outdir, 'recruited_unique_kmers')
    recr_cmd = ['python', '-u', f'{SCRIPT_DIR}/distance_based_kmer_recruitment.py',
                '--ncrf', ncrf_fn,
                '--coverage', str(int(32 * pread * read_reduction)),
                '--min-coverage', '4',
                '--outdir', recr_outdir]
    print(' '.join(recr_cmd))
    subprocess.call(recr_cmd)
    kmers_fn = os.path.join(recr_outdir, 'unique_kmers_min_edge_cov_4.txt')

    # Run read placer
    placer_outdir = os.path.join(outdir, 'tr_resolution')
    placer_cmd = ['python', '-u', f'{SCRIPT_DIR}/read_placer.py',
                  '--ncrf', ncrf_fn,
                  '--genomic-kmers', kmers_fn,
                  '--outdir', placer_outdir]
    print(' '.join(placer_cmd))
    subprocess.call(placer_cmd)
    read_pos_fn = os.path.join(placer_outdir, 'read_positions.csv')

    # Run ELTR polisher
    polisher_outdir = os.path.join(outdir, 'polishing')
    polisher_cmd = ['python', '-u', f'{SCRIPT_DIR}/eltr_polisher.py',
                    '--read-placement', read_pos_fn,
                    '--outdir', polisher_outdir,
                    '--ncrf', ncrf_fn,
                    '--output-progress',
                    '--error-mode', 'nano',
                    '--num-iters', '4',
                    '--num-threads', '50',
                    '--unit', unit_fn]
    print(' '.join(polisher_cmd))
    subprocess.call(polisher_cmd)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--reads", required=True)
    parser.add_argument("-o", "--outdir", required=True)
    params = parser.parse_args()

    all_reads = read_bio_seqs(params.reads)
    smart_makedirs(params.outdir)

    for seed in seeds:
        for pread in preads:
            for cut_freq in cut_freqs:
                iter_outdir = \
                    os.path.join(params.outdir,
                                 f'pread_{pread}_cutfreq_{cut_freq}_seed_{seed}')
                smart_makedirs(iter_outdir)
                run(all_reads=all_reads,
                    seed=seed,
                    cut_freq=cut_freq,
                    pread=pread,
                    outdir=iter_outdir)


if __name__ == "__main__":
    main()
