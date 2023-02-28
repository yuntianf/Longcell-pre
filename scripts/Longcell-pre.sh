#!/bin/bash
# coding: utf-8


#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
#                    Default parameters                   #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
### environment
python=python  #path of python
Rscript=Rscript  #path of Rscript
parallel=parallel  #path of GNU parallel
cores=1

### softclips extraction
toolkit=5

### barcode match


### filter reads without barcode
num=4   # the number of split file for the input of UMI deduplication (for better parallization)

### UMI deduplication
UMI_len=10
thresh=5  # threshold for the similarity between UMI (quantified as needlman score, for longer UMI the thresh should be higher)

### output path

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
#                     Input parameters                    #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
showHelp() {
    cat << EOF  
Program: Longcell-pre (The tool to generate single cell isoform count from bam)
Version: 0.1

Usage: Longcell-pre [options]
Options:

    -h, -help,          --help                  Display help

    -b, -bam,           --bam                   Input bam files

    -d, -bed,           --bed                   Input folder for bed annotations

    -w, -whitelist,     --whitelist             Input barcode whitelist to do barcode match
    
    -c, -cores,         --cores                 Number of cores to use for preprocessing
    
    -t, -toolkit,       --toolkit               The 10X toolkit for library preparation (The location of the cell barcode and UMI in the read)
    
    -L, -blen,          --blen                  The length of the cell barcode
    
    -l, -ulen,          --ulen                  The length of the UMI
    
    -o, -outdir,        --outdir                The output directory

EOF
}

options=$(getopt -l "help,bam:,bed:,whitelist:,cores:,toolkit:,blen:,ulen:,outdir:" \
          -o "hb:d:w:c:t:L:l:o:" -a -- "$@")
[ $? -ne 0 ] && exit 1
eval set -- "$options"

while true
do
    case $1 in
        -h|--help) 
            showHelp
            exit 0
            ;;
        -b|--bam) 
            bam=$2
            echo "bam input bam is "$bam
            shift
            ;;
        -d|--bed)
            bed_folder=$2
            echo "bed input bed annotation is "$bed_folder
            shift
            ;;
        -w|--whitelist)
            barcodes=$2
            echo "barcode whitelist input is "$barcodes
            shift
            ;;
        -c|--cores)
            cores=$2
            echo "The number of cores is "$cores
            shift
            ;;
        -t|--toolkit)
            toolkit=$2
            shift
            ;;
        -L|--blen)
            blen=$2
            shift
            ;;
        -o|--outdir)
            outdir=$2
            echo "The output folder is "$outdir
            shift
            ;;
        --)
            shift
            break;;
    esac
    shift
done

if [ -z "$bam" ]; then
        echo 'The input bam file is missing' >&2
        exit 1
fi
if [ -z "$bed_folder" ]; then
        echo 'The input bed folder is missing' >&2
        exit 1
fi
if [ -z "$barcodes" ]; then
        echo 'The input barcode whitelist is missing' >&2
        exit 1
fi
if [ -z "$outdir" ]; then
        echo 'The output folder should be designated' >&2
        exit 1
fi

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
#                   make output folder                    #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
if [ -d "$outdir" ]; then
    echo "outdir folder already exists"
else
    mkdir $outdir
fi

mkdir -p $outdir/exon_reads
mkdir -p $outdir/softclips
mkdir -p $outdir/barcode_match
mkdir -p $outdir/sub_cell_exon
mkdir -p $outdir/cell_gene_splice_count
mkdir -p $outdir/cell_gene_exon_count

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
#       Extract softclips and splice sites                #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
script=./SoftclipsExon/softclip_splicesite.py
find $bed_folder -name EN* | $parallel -j $cores $python $script -b $bam -t $toolkit -g {} -o $outdir/exon_reads/
find $outdir/exon_reads/ -name "*" -type f -size 0c | xargs -n 1 rm -f
cat $outdir/exon_reads/EN* > $outdir/exon_reads/exon_reads.txt
cut -f 1,2 $outdir/exon_reads/exon_reads.txt > $outdir/softclips/softclips.txt

echo 'Extraction of softclips and splice sites finished!'

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
#                     Barcode Match                       #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
script=./BarcodeMatch/BarcodeMatch.py
$python $script -q $outdir/softclips/softclips.txt -c $barcodes -o "$outdir/barcode_match/bc.txt" -co $cores
echo 'Barcode Match finished!'

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
#             filter reads without barcode                #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
script=./BarcodeMatch/barcode_merge.R
bc=$outdir/barcode_match/bc.txt
exon=$outdir/exon_reads/exon_reads.txt
$Rscript $script $bc $exon $num $outdir/sub_cell_exon/

echo 'Reads filtering finished!'

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
#                   UMI deduplication                     #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
script=./UmiDedup/umi_count.R
input=$outdir/sub_cell_exon/
out=$outdir/cell_gene_splice_count/


find $input -name sub*.txt | $parallel -j $cores $Rscript $script -c {} -u $UMI_len -s $thresh -o $out/sub_cell_gene_splice_count.{#}.txt
awk 'FNR>1 || NR==1' $out/sub_cell_gene_splice_count.*.txt > $out/cell_gene_splice_count.txt

echo 'UMI deuplication finished!'


#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
#         exon bins to exon id and write to file          #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
script=./spliceob/createExonList.R
input=$outdir/cell_gene_splice_count/
bed=$bed_folder/gene_bed.rds
out=$outdir/cell_gene_exon_count/

find $input -name sub_cell_gene_splice_count.*.txt | $parallel -j $cores $Rscript $script {} $bed $out/sub_cell_gene_exon_count.{#}.txt
awk 'FNR>1 || NR==1' $out/sub_cell_gene_exon_count.*.txt > $out/cell_gene_exon_count.txt

echo 'exon transformation finished!'

script=./spliceob/saveExonList.R
input=$outdir/cell_gene_exon_count/cell_gene_exon_count.txt
out=$outdir/cell_gene_exon_count/

$Rscript $script $input $out
echo 'data saving finished!'

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
#                     output cleaning                     #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
rm -rf $outdir/exon_reads/
rm -rf $outdir/softclips/
rm -rf $outdir/sub_cell_exon/

rm -rf $outdir/cell_gene_splice_count/sub_cell_gene_splice_count.*.txt
rm -rf $outdir/cell_gene_exon_count/sub_cell_gene_exon_count.*.txt

echo 'All process finished!'
