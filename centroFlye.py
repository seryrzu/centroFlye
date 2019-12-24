import argparse
import math
import os

import shutil
import subprocess
import sys

from scripts.utils.os_utils import smart_makedirs
from scripts.utils.various import list2str, listEls2str


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

    # unit consensus
    unit_cons_args = parser.add_argument_group("Unit consensus reconstruction")
    unit_cons_args.add_argument("--cons-k-mer-len",
                                help="Length of frequent k-mer",
                                default=30,
                                type=int)
    # Polisher
    polisher_args = parser.add_argument_group("Polisher")
    polisher_args.add_argument("--flye-bin",
                               help='Default path for Flye',
                               default='flye')
    polisher_args.add_argument("--error-mode",
                               help='Error mode: nano/pacbio',
                               default="nano")
    polisher_args.add_argument("--num-polish-iters",
                               help='Number of iterations',
                               default=4,
                               type=int)
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


class CentroFlye:
    def __init__(self, params):
        self.params = params
        MAIN_DIR = os.path.dirname(os.path.realpath(__file__))
        SCRIPTS_DIR = os.path.join(MAIN_DIR, 'scripts')
        self.ncrf_script = os.path.join(SCRIPTS_DIR, "run_ncrf_parallel.py")
        self.recr_script = os.path.join(SCRIPTS_DIR,
                                        "distance_based_kmer_recruitment.py")
        self.placer_script = os.path.join(SCRIPTS_DIR, "read_placer.py")
        self.unit_reconstructor = \
            os.path.join(SCRIPTS_DIR,
                         "better_consensus_unit_reconstruction.py")
        self.polisher_script = os.path.join(SCRIPTS_DIR, "eltr_polisher.py")
        self.tandemPolisher = os.path.join(SCRIPTS_DIR,
                                           "ext",
                                           "tandemQUAST",
                                           "tandemquast.py")

    def run_NCRF(self):
        ncrf_outdir = os.path.join(self.params.outdir, "NCRF")
        ncrf_cmd = ["python", "-u", self.ncrf_script,
                    "--reads", self.params.reads,
                    "--threads", self.params.threads,
                    "--outdir", ncrf_outdir,
                    "--repeat", self.params.unit]
        ncrf_cmd = listEls2str(ncrf_cmd)
        print('Running NCRF:')
        print(list2str(ncrf_cmd))
        subprocess.call(ncrf_cmd)
        ncrf_fn = os.path.join(ncrf_outdir, "report.ncrf")
        return ncrf_fn

    def run_kmer_recr(self, ncrf_fn):
        recr_outdir = os.path.join(self.params.outdir,
                                   "recruited_unique_kmers")
        recr_cmd = ["python", "-u", self.recr_script,
                    "--ncrf", ncrf_fn,
                    "--coverage", self.params.coverage,
                    "--min-coverage", self.params.min_coverage,
                    "--outdir", recr_outdir,
                    "--min-nreads", self.params.min_nreads,
                    "--max-nreads", self.params.max_nreads,
                    "--min-distance", self.params.min_distance,
                    "--max-distance", self.params.max_distance,
                    "--bottom", self.params.bottom,
                    "--top", self.params.top,
                    "--kmer-survival-rate", self.params.kmer_survival_rate,
                    "--max-nonuniq", self.params.max_nonuniq]
        recr_cmd = listEls2str(recr_cmd)
        print("Running unique kmer recruitment")
        print(list2str(recr_cmd))
        subprocess.call(recr_cmd)
        kmers_fn = os.path.join(
            recr_outdir,
            f"unique_kmers_min_edge_cov_{self.params.min_coverage}.txt")
        return kmers_fn

    def run_read_placer(self, ncrf_fn, kmers_fn):
        placer_outdir = os.path.join(self.params.outdir, "tr_resolution")
        placer_cmd = ["python", "-u", self.placer_script,
                      "--ncrf", ncrf_fn,
                      "--genomic-kmers", kmers_fn,
                      "--outdir", placer_outdir,
                      "--n-motif", self.params.n_motif,
                      "--min-cloud-kmer-freq", self.params.min_cloud_kmer_freq,
                      "--min-kmer-mult", self.params.min_kmer_mult,
                      "--min-unit", self.params.min_unit,
                      "--min-inters", self.params.min_inters]
        placer_cmd = listEls2str(placer_cmd)
        print("Running tandem repeat resolution")
        print(list2str(placer_cmd))
        subprocess.call(placer_cmd)
        read_pos_fn = os.path.join(placer_outdir, "read_positions.csv")
        return read_pos_fn

    def run_unit_reconstructor(self, ncrf_fn):
        unit_star_fn = os.path.join(self.params.outdir,
                                    'cons_unit',
                                    'unit_star.fasta')
        unit_star_cmd = ["python", "-u", self.unit_reconstructor,
                         "--reads-ncrf", ncrf_fn,
                         "--unit", self.params.unit,
                         "-k", self.params.cons_k_mer_len,
                         "--output", unit_star_fn]
        unit_star_cmd = listEls2str(unit_star_cmd)
        print("Running unit star reconstruction")
        print(list2str(unit_star_cmd))
        subprocess.call(unit_star_cmd)
        return unit_star_fn

    def run_polisher(self, ncrf_fn, read_pos_fn, unit_star_fn):
        polisher_outdir = os.path.join(self.params.outdir, "polishing1")
        polisher_cmd = ["python", "-u", self.polisher_script,
                        "--read-placement", read_pos_fn,
                        "--outdir", polisher_outdir,
                        "--ncrf", ncrf_fn,
                        "--error-mode", self.params.error_mode,
                        "--num-iters", self.params.num_polish_iters,
                        "--num-threads", self.params.threads,
                        "--unit", unit_star_fn,
                        "--min-pos", self.params.min_pos]
        if self.params.max_pos != math.inf:
            polisher_cmd.append("--max-pos")
            polisher_cmd.append(self.params.max_pos)
        polisher_cmd = listEls2str(polisher_cmd)
        print("Running polishing -- first stage")
        print(list2str(polisher_cmd))
        subprocess.call(polisher_cmd)
        assembly_fn = \
            os.path.join(
                polisher_outdir,
                f'final_sequence_{self.params.num_polish_iters}.fasta'
            )
        return assembly_fn

    def run_tandemPolisher(self, assembly_fn):
        polisher_outdir = os.path.join(self.params.outdir, "polishing2")
        polisher_cmd = ["python", "-u", self.tandemPolisher,
                        "-r", self.params.reads,
                        "-o", polisher_outdir,
                        "--only-polish",
                        "-t", self.params.threads,
                        assembly_fn]
        polisher_cmd = listEls2str(polisher_cmd)
        print("Running polishing with tandemMapper")
        print(list2str(polisher_cmd))
        subprocess.call(polisher_cmd)
        polished_assembly_fn = \
            os.path.join(polisher_outdir, 'polished', 'polished_2.fasta')
        return polished_assembly_fn

    def copy_final_assembly(self, polished_assembly_fn):
        final_assembly_fn = os.path.join(self.params.outdir,
                                         'final_assembly.fasta')
        shutil.copyfile(polished_assembly_fn, final_assembly_fn)
        print(f"Final polished assembly is stored at {final_assembly_fn}")

    def run(self):
        smart_makedirs(self.params.outdir)
        ncrf_fn = self.run_NCRF()
        kmers_fn = self.run_kmer_recr(ncrf_fn)
        read_pos_fn = self.run_read_placer(ncrf_fn=ncrf_fn, kmers_fn=kmers_fn)
        unit_star_fn = self.run_unit_reconstructor(ncrf_fn=ncrf_fn)
        assembly_fn = self.run_polisher(ncrf_fn=ncrf_fn,
                                        read_pos_fn=read_pos_fn,
                                        unit_star_fn=unit_star_fn)
        polished_assembly_fn = self.run_tandemPolisher(assembly_fn=assembly_fn)
        self.copy_final_assembly(polished_assembly_fn)


def main():
    params = parse_args()
    CentroFlye(params).run()


if __name__ == "__main__":
    main()
