#!/usr/bin/env bash
set -e
set -o xtrace

#######
#USAGE:
path_to_CHM13=$1
OUTPUT=$2
THREADS=$3

${BASH_SOURCE%/*}/scripts/read_recruitment/run_read_recruitment.sh \
    ${path_to_CHM13}/rel3_bgzip.fastq.gz \
    ${OUTPUT}/centromeric_reads ${THREADS} 29000000

python ${BASH_SOURCE%/*}/scripts/run_ncrf_parallel.py \
            --reads ${OUTPUT}/centromeric_reads/centromeric_reads.fasta \
            -t ${THREADS} \
            --outdir ${OUTPUT}/NCRF_rc \
            --repeat supplementary_data/DXZ1_rc.fasta

python ${BASH_SOURCE%/*}/scripts/distance_based_kmer_recruitment.py \
        --ncrf ${OUTPUT}/NCRF_rc/report.ncrf \
        --coverage 32 \
        --min-coverage 4 \
        --outdir ${OUTPUT}/recruited_unique_kmers

python ${BASH_SOURCE%/*}/scripts/read_placer.py \
              --ncrf ${OUTPUT}/NCRF_rc/report.ncrf \
              --genomic-kmers ${OUTPUT}/recruited_unique_kmers/unique_kmers_min_edge_cov_4.txt \
              --outdir ${OUTPUT}/tr_resolution

python ${BASH_SOURCE%/*}/scripts/better_consensus_unit_reconstruction.py \
              --reads-ncrf ${OUTPUT}/NCRF_rc/report.ncrf \
              --unit supplementary_data/DXZ1_rc.fasta \
              --output ${OUTPUT}/DXZ1_star/DXZ1_rc_star.fasta

python ${BASH_SOURCE%/*}/scripts/eltr_polisher.py \
              --read-placement ${OUTPUT}/tr_resolution/read_positions.csv \
              --outdir ${OUTPUT}/polishing \
              --ncrf ${OUTPUT}/NCRF_rc/report.ncrf \
              --output-progress \
              --error-mode nano \
              --num-iters 4 \
              --num-threads ${THREADS} \
              --unit ${OUTPUT}/DXZ1_star/DXZ1_rc_star.fasta
