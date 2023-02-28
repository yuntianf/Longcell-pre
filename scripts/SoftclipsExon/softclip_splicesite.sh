#!/bin/bash
#$ -N softclip_exon_extraction
#$ -o /home/stat/yuntianf/output
#$ -j y


#module load python/python-3.3
conda activate py3.7

data=FT_stim

script=$HOME/scripts/BarcodeMatch/softclips_exon/softclip_splicesite.py
bed=$HOME/data/BarcodeMatch/gene_bed_update/
toolkit='5'
bam=$HOME/LongReadsdata/LongReads_guppy5/${data}_sub/${data}.bam
out_folder=$HOME/interdata/LongReads/data/BarcodeMatch/${data}/${data}_softclips_exon/exon_reads/


if [[ -n $SGE_TASK_ID ]]; then
    F=$(ls -1 $bed/*.txt | sed -n ${SGE_TASK_ID}p)
else
    echo "need a file number ... perhaps run with 'qsub -t 1-N softclip_splicesite.sh'"
    exit 1
fi


python $script -b $bam -t $toolkit -g $F -o $out_folder
