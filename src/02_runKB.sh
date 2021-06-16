#!/bin/bash

source activate single_cell_env

#kb --help

export INDEX=${1}
export GENO=${2}
export OUTDIR=${3}
export R1=${4}
export R2=${5}

echo 'Index file : '$INDEX
echo 'Transcript to gene file : '$GENO
echo $' '
echo 'Fastq 1 : '$R1
echo 'Fastq 2 : '$R2
echo $' '
echo 'Output directory : '$OUTDIR
echo $' '

# locale
export LANGUAGE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

cd $OUTDIR

kb count -i $INDEX --keep-tmp -g $GENO -o $OUTDIR -x 10xv3 -t 15 --verbose $R1 $R2 
#kb count -i $INDEX --keep-tmp -g $GENO -o $OUTDIR -x 10xv3 -t 10 --verbose --report $R1 $R2
