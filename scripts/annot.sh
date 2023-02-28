#!/bin/bash
# coding: utf-8

set -e

Rscript=Rscript
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
#                    Input parameters                     #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
gtf=gencode.v39.annotation.gtf
bed_folder=../gene_bed/

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
#                    Build Annotation                     #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+# 
script=./Auxiliary/gtf2bed.R
if [ -d "$bed_folder" ]; then
    echo "$bed_folder folder already exists"
else 
    mkdir $bed_folder
fi

if [ "$(ls -A $bed_folder)" ]; then
    echo "$bed_folder is not Empty, will skip the gtf to bed process!"
else
    $Rscript $script $gtf $bed_folder
fi

echo 'Build bed annotation finished!'
