#!/bin/bash
set -e
set -o xtrace

#######
#USAGE:
NCRF=$1
OUTPUT=$2
#######

MIN_COVERAGE=5
COVERAGE=32

# taken from https://stackoverflow.com/a/246128
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

mkdir -p $OUTPUT

python $SCRIPT_DIR/../distance_based_kmer_recruitment.py --ncrf $NCRF --coverage $COVERAGE --min-coverages $MIN_COVERAGE --outdir $OUTPUT/first_1000 --max-nreads 1000
python $SCRIPT_DIR/../distance_based_kmer_recruitment.py --ncrf $NCRF --coverage $COVERAGE --min-coverages $MIN_COVERAGE --outdir $OUTPUT/second_1000 --min-nreads 1000 --max-nreads 2000

cat $OUTPUT/first_1000/unique_kmers_min_edge_cov_${MIN_COVERAGE}_RC.txt $OUTPUT/second_1000/unique_kmers_min_edge_cov_${MIN_COVERAGE}_RC.txt | sort | uniq > $OUTPUT/unique_kmers_min_edge_cov_${MIN_COVERAGE}_RC.txt
