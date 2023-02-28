#!/bin/bash
#$ -N umi_count
#$ -l m_mem_free=100G
#$ -o /home/stat/yuntianf/output
#$ -j y


#module load python/python-3.3
conda activate R

data=FT

script=$HOME/scripts/UMI_dedup/dedup/umi_count.R
input=$HOME/interdata/LongReads/data/BarcodeMatch/${data}/${data}_softclips_exon/sub_cell_exon/
thresh=5
out_folder=$HOME/interdata/LongReads/data/BarcodeMatch/${data}/${data}_softclips_exon/cell_gene_splice_count/


if [[ -n $SGE_TASK_ID ]]; then
    F=$(ls -1 $input/sub_cell_exon.* | sed -n ${SGE_TASK_ID}p)
else
    echo "need a file number ... perhaps run with 'qsub -t 1-N umi_count.sh'"
    exit 1
fi

prefix=sub_cell_gene_splice_count
file="${prefix}.${F#*.}"

Rscript $script $F $bc $thresh "${out_folder}/${file}"
