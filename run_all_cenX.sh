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
    ${path_to_CHM13}/rel2.fastq.gz \
    ${OUTPUT}/centromeric_reads ${THREADS} 11100000

python ${BASH_SOURCE%/*}/centroFlye.py \
    --coverage 32 \
    --reads ${OUTPUT}/centromeric_reads/centromeric_reads.fasta \
    -t ${THREADS} \
    --outdir ${OUTPUT} \
    --unit supplementary_data/DXZ1_rc.fasta
