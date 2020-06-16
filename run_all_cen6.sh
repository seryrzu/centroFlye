#!/usr/bin/env bash
set -e
set -o xtrace

#######
#USAGE:
path_to_CHM13=$1
OUTPUT=$2
THREADS=$3

make -C ${BASH_SOURCE%/*}/scripts/read_recruitment

${BASH_SOURCE%/*}/scripts/read_recruitment/run_read_recruitment.sh \
    ${path_to_CHM13}/rel3.fastq.gz \
    ${OUTPUT}/centromeric_reads ${THREADS} 29000000 \
    ${BASH_SOURCE%/*}/supplementary_data/D6Z1.fasta \
    550

python ${BASH_SOURCE%/*}/scripts/ext/stringdecomposer/longreads_decomposer.py \
    -s ${OUTPUT}/centromeric_reads/centromeric_reads.fasta \
    -m ${BASH_SOURCE%/*}/supplementary_data/D6Z1_monomers.fasta \
    -t 50
mkdir ${OUTPUT}/string_decomposer_report
mv ${BASH_SOURCE%/*}/decomposition{.tsv,_alt.tsv} ${OUTPUT}/string_decomposer_report

python ${BASH_SOURCE%/*}/scripts/centroFlyeMono.py \
    --sd-report ${OUTPUT}/string_decomposer_report/decomposition.tsv \
    --monomers ${BASH_SOURCE%/*}/supplementary_data/D6Z1_monomers.fasta \
    --centromeric-reads ${OUTPUT}/centromeric_reads/centromeric_reads.fasta \
    --outdir ${OUTPUT}/centroFlyeMono_cen6
