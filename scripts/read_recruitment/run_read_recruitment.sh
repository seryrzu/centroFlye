#!/usr/bin/env bash
set -e
set -o xtrace

#######
#USAGE:
INPUT=$1
OUTPUT=$2
THREADS=$3
NREADS=$4
UNIT=$5


# TODO: update the "output" dir
# scripts/read_recruitment/run_read_recruitment.sh <path2CHM13>/rel2.fastq.gz results/centromeric_reads_test_2 50 11100000
#######

# taken from https://stackoverflow.com/a/246128
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

if [ -z "$UNIT" ] ; then
    UNIT="${SCRIPT_DIR}/../../supplementary_data/DXZ1_rc.fasta"
else
    UNIT="$(realpath $5)"
    echo $UNIT
fi

mkdir -p $OUTPUT

nreads_per_thread="$((NREADS / THREADS * 2))"
 
mkdir -p $OUTPUT/split_fa
zcat $INPUT | seqtk seq -A | awk -v x="$nreads_per_thread" -v y="$OUTPUT/split_fa" 'BEGIN {n_seq=0;} {if(n_seq%x==0){file=sprintf("%s/split_fasta_%d.fasta",y,n_seq);} print >> file; n_seq++; next;} { print >> file; }'

cd $OUTPUT/split_fa
mkdir -p ../split_rr
find . -name "*.fasta" | xargs -I {} -P $THREADS $SCRIPT_DIR/rr $UNIT {} ../split_rr/{}_cen.fasta 350
cd ../split_rr
cat *.fasta > ../centromeric_reads.fasta
