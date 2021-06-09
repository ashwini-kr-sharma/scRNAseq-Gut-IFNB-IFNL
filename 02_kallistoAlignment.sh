#!/bin/bash

ROOT="/icgc/dkfzlsdf/analysis/B080/sharma/boulantLab/sc_interferon"

#-------------------------------------------------------------------------------
# Create Conda virtual environment for Kallisto
#-------------------------------------------------------------------------------

# Make sure miniconda is installed - 
# see https://docs.conda.io/en/latest/miniconda.html

# conda create --name single_cell_env python=3.8 pip -y
# conda activate single_cell_env
# pip install kb-python
# pip install jupyter
#-------------------------------------------------------------------------------
# Download the reference transcriptome index and transcript to gene mapping
#-------------------------------------------------------------------------------

mkdir -p $ROOT$'/data/kallistoIndex'

if [[ -f $ROOT$'/data/kallistoIndex/index.idx' ]]; then

    echo "The index and transcript to gene files exists !! Doing nothing ..."
    
else
    echo "Downloading the index and transcript-to-gene files ..."
    
    cd $ROOT$'/data/kallistoIndex'

    kb ref -d human \
    --keep-tmp \
    --workflow 'standard' \
    -i $ROOT$'/data/kallistoIndex/index.idx' \
    -g $ROOT$'/data/kallistoIndex/t2g.txt' \
    --overwrite
fi

#-------------------------------------------------------------------------------
# Kalisto | bustools
#-------------------------------------------------------------------------------

conditions=('Int1' 'Int2' 'Int3' 'Int4' 'Int5' 'Int6' 'Int7' 'Int8')
techrep=('HJJCYBGXF' 'HJHG3BGXF')

SCRIPT=$ROOT'/src/merged_after_kb/02_runKB.sh'
INDEX=$ROOT'/data/kallistoIndex/index.idx'
GENO=$ROOT'/data/kallistoIndex/t2g.txt'

for TR in "${techrep[@]}"; do

  if [[ $TR == 'HJJCYBGXF' ]]; then
      VAL=1
  fi

  if [[ $TR == 'HJHG3BGXF' ]]; then
      VAL=2
  fi

  for TYPE in "${conditions[@]}"; do
  
    OUTDIR=$ROOT'/analysis/'$TR'_counts/kb/'$TYPE
    LOGDIR=$ROOT'/logs/alignment/'$TR'/'
    
    mkdir -p $OUTDIR
    mkdir -p $LOGDIR
    
    R1=$ROOT$'/data/sc_Interferon/Int_'$TR'/'$TR'_Pool_Interferon_20s001227-'$VAL'-1_Triana_lane1'$TYPE'_1_sequence.txt.gz'
    R2=$ROOT$'/data/sc_Interferon/Int_'$TR'/'$TR'_Pool_Interferon_20s001227-'$VAL'-1_Triana_lane1'$TYPE'_2_sequence.txt.gz'
    
    bsub -n 15 -R "rusage[mem=10G]" -W 10:00 -J $TYPE'-kb' \
    -o $LOGDIR$TYPE'_kb_python.out' \
    -e $LOGDIR$TYPE'_kb_python.err' \
    "$SCRIPT $INDEX $GENO $OUTDIR $R1 $R2"
  done
done
