import argparse
import sys
import math
# import os
# import random
# import subprocess
#
# from scripts.utils.bio import read_bio_seqs, write_bio_seqs
# from scripts.utils.os_utils import smart_makedirs
# from scripts.utils.trim_seqs import trim_seqs


# SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))


def parse_args():
    parser = argparse.ArgumentParser()

    # General
    gen_args = parser.add_argument_group("General")
    gen_args.add_argument("-i",
                          "--reads",
                          help="centromeric reads",
                          required=True)
    gen_args.add_argument("-o",
                          "--outdir",
                          help="output dir",
                          required=True)
    gen_args.add_argument("-u",
                          "--unit",
                          help="HOR (DXZ1 for rc)",
                          required=True)
    gen_args.add_argument("-c",
                          "--coverage",
                          help="Coverage with ultra-long (50kb+) reads",
                          type=int,
                          required=True)
    gen_args.add_argument("-t",
                          "--threads",
                          help="Number of threads",
                          type=int,
                          default=50)
    gen_args.add_argument("-k",
                          "--k-mer-len",
                          help="Length of k-mer",
                          type=int,
                          default=19)

    # Recruitment of unique k-mers
    recr_args = parser.add_argument_group("Unique k-mer Recruitment")
    recr_args.add_argument("--min-coverage",
                           help="minCov of an edge in Distance Graph",
                           type=int,
                           default=4)
    recr_args.add_argument('--min-nreads',
                           type=int,
                           default=0)
    recr_args.add_argument('--max-nreads',
                           type=int,
                           default=sys.maxsize)
    recr_args.add_argument('--min-distance',
                           help='Min distance (in units) to consider',
                           type=int,
                           default=1)
    recr_args.add_argument('--max-distance',
                           type=int,
                           help='Max distance (in units) to consider',
                           default=150)
    recr_args.add_argument('--bottom',
                           type=float,
                           help='Bottom parameter for rare k-mer recruitment',
                           default=0.9)
    recr_args.add_argument('--top',
                           type=float,
                           help='Top parameter for rare k-mer recruitment',
                           default=3.)
    recr_args.add_argument('--kmer-survival-rate',
                           type=float,
                           help='Estimated k-mer survival rate in reads',
                           default=0.34)
    recr_args.add_argument('--max-nonuniq',
                           help='Max #occurences of a kmer in a read',
                           type=int,
                           default=3)

    # Read placer
    placer_args = parser.add_argument_group("Read Placer args")
    placer_args.add_argument("--n-motif",
                             help="Number of motifs stuck together",
                             default=1,
                             type=int)
    placer_args.add_argument("--min-cloud-kmer-freq",
                             help="Minimal frequency of a kmer in the cloud",
                             default=2,
                             type=int)
    placer_args.add_argument("--min-kmer-mult",
                             help="Minimal frequency of a kmer in input",
                             default=2,
                             type=int)
    placer_args.add_argument("--min-unit",
                             help="Score[0]",
                             default=2,
                             type=int)
    placer_args.add_argument("--min-inters",
                             help="Score[1]",
                             default=10,
                             type=int)

    # Polisher
    polisher_args = parser.add_argument_group("Polisher")
    polisher_args.add_argument("--flye-bin",
                               help='Default path for Flye',
                               default='flye')
    polisher_args.add_argument("--error-mode",
                               help='Error mode: nano/pacbio',
                               default="pacbio")
    polisher_args.add_argument("--num-iters",
                               help='Number of iterations',
                               default=2,
                               type=int)
    polisher_args.add_argument("--output-progress",
                               help='Output progress of polisher',
                               action='store_false')
    polisher_args.add_argument("--min-pos",
                               help='Min position unit to polish',
                               type=int,
                               default=0)
    polisher_args.add_argument("--max-pos",
                               help='Max position unit to polish',
                               type=int,
                               default=math.inf)

    params = parser.parse_args()
    return params


# def run(all_reads, seed, cut_freq, pread, outdir):
#     print("New run")
#     print(seed, cut_freq, pread)
#     random.seed(seed)
#
#     nreads=int(len(all_reads)*pread)
#     read_ids=random.sample(all_reads.keys(), nreads)
#     reads={r_id: all_reads[r_id] for r_id in read_ids}
#     reads=trim_seqs(reads, cut_freq)
#     read_reduction=sum(len(v) for v in reads.values()) / \
#                      sum(len(all_reads[r_id]) for r_id in reads.keys())
#     print("read_reduction", read_reduction)
#     reads_fn=os.path.join(outdir, "centromeric_reads.fasta")
#     print(reads_fn)
#     write_bio_seqs(reads_fn, reads)
#
#     # Run NCRF
#     ncrf_outdir=os.path.join(outdir, "NCRF_rc")
#     ncrf_cmd=["python", "-u", f"{SCRIPT_DIR}/run_ncrf_parallel.py",
#                 "--reads", reads_fn,
#                 "-t", "50",
#                 "--outdir", ncrf_outdir,
#                 "--repeat", unit_fn]
#     print(" ".join(ncrf_cmd))
#     subprocess.call(ncrf_cmd)
#     ncrf_fn=os.path.join(ncrf_outdir, "report.ncrf")
#
#     # Run k-mer recruitment
#     recr_outdir=os.path.join(outdir, "recruited_unique_kmers")
#     recr_cmd=["python", "-u", f"{SCRIPT_DIR}/distance_based_kmer_recruitment.py",
#                 "--ncrf", ncrf_fn,
#                 "--coverage", str(int(32 * pread * read_reduction)),
#                 "--min-coverage", "4",
#                 "--outdir", recr_outdir]
#     print(" ".join(recr_cmd))
#     subprocess.call(recr_cmd)
#     kmers_fn=os.path.join(recr_outdir, "unique_kmers_min_edge_cov_4.txt")
#
#     # Run read placer
#     placer_outdir=os.path.join(outdir, "tr_resolution")
#     placer_cmd=["python", "-u", f"{SCRIPT_DIR}/read_placer.py",
#                   "--ncrf", ncrf_fn,
#                   "--genomic-kmers", kmers_fn,
#                   "--outdir", placer_outdir]
#     print(" ".join(placer_cmd))
#     subprocess.call(placer_cmd)
#     read_pos_fn=os.path.join(placer_outdir, "read_positions.csv")
#
#     # Run ELTR polisher
#     polisher_outdir=os.path.join(outdir, "polishing")
#     polisher_cmd=["python", "-u", f"{SCRIPT_DIR}/eltr_polisher.py",
#                     "--read-placement", read_pos_fn,
#                     "--outdir", polisher_outdir,
#                     "--ncrf", ncrf_fn,
#                     "--output-progress",
#                     "--error-mode", "nano",
#                     "--num-iters", "4",
#                     "--num-threads", "50",
#                     "--unit", unit_fn]
#     print(" ".join(polisher_cmd))
#     subprocess.call(polisher_cmd)


def main():
    params = parse_args()

    # all_reads = read_bio_seqs(params.reads)
    # smart_makedirs(params.outdir)

    # for seed in seeds:
    #     for pread in preads:
    #         for cut_freq in cut_freqs:
    #             if pread == 1 and cut_freq == 0:
    #                 continue
    #             iter_outdir = \
    #                 os.path.join(params.outdir,
    #                              f"pread_{pread}_cutfreq_{cut_freq}_seed_{seed}")
    #             smart_makedirs(iter_outdir)
    #             run(all_reads=all_reads,
    #                 seed=seed,
    #                 cut_freq=cut_freq,
    #                 pread=pread,
    #                 outdir=iter_outdir)


if __name__ == "__main__":
    main()
