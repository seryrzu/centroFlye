import argparse
import os

from debruijn_graph import iterative_graph, scaffolding, read2scaffolds,\
                           cover_scaffolds_w_reads, extract_read_pseudounits, \
                           polish
from sd_parser_83640e3 import SD_Report
from mono_error_correction import error_correction

from utils.os_utils import smart_makedirs
from utils.bio import read_bio_seqs


SD_version = 83640e3


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sd-report',
                        help='Path to SD report (v.{SD_version})',
                        required=True)
    parser.add_argument('--monomers',
                        help='Path to fasta with monomers used to call SD',
                        required=True)
    parser.add_argument('--centromeric-reads',
                        help='Path to centromeric reads',
                        required=True)
    parser.add_argument('--outdir',
                        help='Output directory',
                        required=True)
    parser.add_argument('--min-k',
                        help='Min k in De DeBruijnGraph',
                        type=int,
                        default=100)
    parser.add_argument('--max-k',
                        help='Max k in De DeBruijnGraph',
                        type=int,
                        default=400)
    parser.add_argument('--min-mult',
                        help='Min mult for a frequent mono-k-mer',
                        type=int,
                        default=5)
    params = parser.parse_args()
    return params


def main():
    params = parse_args()
    smart_makedirs(params.outdir)

    print('Reading report')
    sd_report = SD_Report(SD_report_fn=params.sd_report,
                          monomers_fn=params.monomers)

    print('Error correcting monoreads')
    ec_monostrings = error_correction(sd_report.monostrings,
                                      verbose=True,
                                      inplace=False)

    print('Building the graph')
    contigs, dbs, all_frequent_kmers, all_frequent_kmers_read_pos = \
        iterative_graph(ec_monostrings,
                        min_k=params.min_k,
                        max_k=params.max_k,
                        outdir=os.path.join(params.outdir, 'idb'),
                        min_mult=params.min_mult)
    db = dbs[params.max_k]


    print('Mapping reads to the graph')
    mappings = db.map_reads(ec_monostrings, verbose=False)

    print('Scaffolding')
    scaffolds, edge_scaffolds = scaffolding(db, mappings)

    # Manual connection of two scaffolds for cen6
    # TODO
    cen6_scaffold = scaffolds[1] + scaffolds[0][db.k-1:]
    cen6_edge_scaffold = edge_scaffolds[1] + edge_scaffolds[0]

    print('Mapping reads to scaffolds')
    r2s = read2scaffolds(db, [cen6_edge_scaffold], mappings, ec_monostrings)

    print('Covering scaffolds with reads')
    scaf_read_coverage = cover_scaffolds_w_reads(r2s,
                                                 mappings,
                                                 [cen6_scaffold],
                                                 ec_monostrings, k=db.k)

    print('Extracting pseudounits and reads covering them')
    pseudounits, read_pseudounits = \
        extract_read_pseudounits(scaf_read_coverage,
                                 [cen6_scaffold],
                                 monostrings=ec_monostrings)

    print('Reading centromeric reads')
    centromeric_reads = read_bio_seqs(params.centromeric_reads)
    monomers = read_bio_seqs(params.monomers)

    print('Polishing')
    polish(scaffolds=[cen6_scaffold],
           pseudounits=pseudounits,
           read_pseudounits=read_pseudounits,
           reads=centromeric_reads,
           monomers=monomers,
           outdir='polishing')


if __name__ == "__main__":
    main()
