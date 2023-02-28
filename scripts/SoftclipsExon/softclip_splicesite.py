#!/usr/bin/env python


import argparse, sys, os, re, pysam, csv, gzip, string,json
import numpy as np, pandas as pd

DNA_COMPLEMENT = str.maketrans("ACGT", "TGCA") 

SEARCH_OFFSET = 0
MAX_SEARCH_LEN_5 = 55-SEARCH_OFFSET  
MAX_SEARCH_LEN_3 = 55-SEARCH_OFFSET  
QUAL_THRESH = 30

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Define reusable modules for softclips extraction     #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+# 

def parse_commandline():
    parser=argparse.ArgumentParser()
    parser.add_argument('--bam', '-b', help='input bam file', type=str, required=True)
    parser.add_argument('--gene_bed', '-g', help='exon annotation for a gene', type=str, required=True)
    parser.add_argument('--toolkit', '-t', help='the location of barcodes and UMI, 5 or 3 prime', type=str, default='5')
    parser.add_argument('--out_path', '-o', help='output file path', type=str, default = "./")
    args=parser.parse_args()
    print(args, file=sys.stderr)
    return args

def extract_softclips(samrd,offset=0):
    soft_clip_fwd = 'NNNNNNNNN'; soft_clip_rev = 'NNNNNNNNN';
    sam_cigar = samrd.cigartuples
    if sam_cigar and sam_cigar[0][0] == 4:
        soft_clip_len = sam_cigar[0][1]
        soft_clip_fwd = samrd.query_sequence[:soft_clip_len+offset] 
    if sam_cigar and sam_cigar[-1][0] == 4:
        soft_clip_len = sam_cigar[-1][1]
        soft_clip_rev = samrd.query_sequence[(0-soft_clip_len-offset):]
    return {'fwd': soft_clip_fwd, 'rev': soft_clip_rev}

def reverse_complement(seq):
    return seq.translate(DNA_COMPLEMENT)[::-1] 

def polyA_rm(seq):
    loc = [m.start() for m in re.finditer('(?=AAAAAAAAAA)', seq)]
    if len(loc) == 0:
        return(seq)
    else:
        return(seq[max(loc):])
  
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Define reusable modules to identify splicesites for each reads    #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#  

def getBlocks(read,refer_start,refer_end,flank):
    start = read.reference_start+1
    end = start
    bins = read.cigartuples
    blocks = []
    
    for i in bins:
        if i[0]==0 or i[0]==2:
            end = end+i[1]    
        elif i[0] == 3:
            if start >= refer_start-flank and end <= refer_end+flank:
                blocks.append((start,end))
            start = end + i[1]
            end = start
    if start >= refer_start-flank and end <= refer_end+flank:
        blocks.append((start,end))
    return(blocks)

def annoSplicesite(exon):
    sites = []
    for i in exon:
        if len(sites) == 0 or sites[-1] + 1 != int(i[1]): 
            sites.append(int(i[1]))
        sites.append(int(i[2]))
    return(sites)

def splicesiteCheck(seq,anno,strand,inner_flank = 1):    
    if len(seq) == 2:
        return(-1)
    i = 1
    j = 0
    site = 0
    while i < len(seq)-1 and j < len(anno):
        if seq[i]  < anno[j]-inner_flank:
            i = i+1
        elif seq[i] >= anno[j]-inner_flank and seq[i] <= anno[j]+inner_flank:
            site = site+1
            i = i+1
            j = j+1
        else:
            j = j + 1
    return(site/(len(seq)-2))

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Define reusable modules to identify exons for each reads    #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#  

def createExonList(gene_input):
    exon_rds = list(csv.reader(gene_input, delimiter=" "))
    return exon_rds 

def printDebug(i, qname, chkrdnames):
    return True if (i <20 or qname in chkrdnames) else False; 

# function to detect polyA
def polyADetect(seq,bin = 20, count = 15):
    i = 0
    while (i + bin) < len(seq):
        if(seq.count("A",i,i+bin) >= count):
            return(True)
        i = i + 1
    return(False)

script_name = os.path.basename(__file__)
print("Running ", script_name)
  
CHKRD_NAMES = ['0c3ff184-b837-4671-a55a-24c7e5bc3df3', '57891cde-7423-48b7-a341-52ace43eba27']

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Check for valid arguments, and that files exist                             #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Files are assumed to be in current path (or fully qualified path is given)
args = parse_commandline()

sam_file = args.bam
if os.path.isfile(sam_file) and os.access(sam_file, os.R_OK):
    sam_input = pysam.Samfile(sam_file,'rb') if sam_file[-3:] == 'bam' else pysam.Samfile(sam_file,'r')
    sam_fname = os.path.basename(sam_file)
else:
    print("Unable to open sam/bam file for input:", sam_file)
    sys.exit(1)
  
gene_exon_bed = args.gene_bed

if os.path.isfile(gene_exon_bed) and os.access(gene_exon_bed, os.R_OK):
    gene_input = open(gene_exon_bed, 'r')
    gene_fname = os.path.basename(gene_exon_bed)
else:
    print("Unable to open gene_exon bed file for input:", gene_exon_bed)
    sys.exit(1)

  
out_detail = args.out_path + gene_fname[0:-4] + '.exon_reads.txt'
reads_out = open(out_detail, 'w')
reads_csv = csv.writer(reads_out, delimiter="\t", quoting=csv.QUOTE_MINIMAL)

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Read sam file and check if spans complete genomic coordinate range (from .bed)
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
i = 0; tot_rds = 0;
adapter_flank_bp = 6; last_qname = 'none';

gene_exons = createExonList(gene_input)
print(gene_exons)

gene_chr = gene_exons[0][0]
gene_start = int(gene_exons[0][1])
gene_end = int(gene_exons[-1][2])
gene_strand = 1 if gene_exons[0][4] == "+" else -1
gene_splicesite = annoSplicesite(gene_exons)

FLANK = 200

for samrd in sam_input.fetch(gene_chr,gene_start, gene_end): 
    tot_rds += 1
    if samrd.is_unmapped or samrd.is_secondary or samrd.mapping_quality < QUAL_THRESH:
        continue
    else:
        chrnm = sam_input.getrname(samrd.tid)
        if chrnm == gene_chr:
            if samrd.qname in CHKRD_NAMES or (samrd.reference_start <= gene_end and samrd.reference_end >= gene_start):
                i += 1
                if printDebug(i, samrd.qname, CHKRD_NAMES):
                    print("{2} Full transcript start: {0}, end: {1}.  Alignment per exon below".format(samrd.reference_start, samrd.reference_end, samrd.qname))

                ###splice site identification
                transcript_blocks = getBlocks(samrd,gene_start,gene_end,FLANK)

                if printDebug(i, samrd.qname, CHKRD_NAMES):
                    print(transcript_blocks)

                if len(transcript_blocks) == 0:
                    continue

                transcript_site = splicesiteCheck(list(sum(transcript_blocks,())), gene_splicesite,gene_strand) 

                ### softclips extraction    
                soft_clips = extract_softclips(samrd)

                if gene_strand == -1:
                    polyA_exist = polyADetect(reverse_complement(soft_clips['fwd']))
                else:
                    polyA_exist = polyADetect(soft_clips['rev'])

                if samrd.reference_end >= gene_end + FLANK:
                    polyA_exist = False

                if args.toolkit == '3':
                    soft_clips['fwd'] = polyA_rm(reverse_complement(soft_clips['fwd']))
                    soft_clips['rev'] = polyA_rm(soft_clips['rev'])
    
                    soft_clips['fwd'] = reverse_complement(soft_clips['fwd'])

                sc_5prime_len = len(soft_clips['fwd'])
                sc_3prime_len = len(soft_clips['rev'])
 
                i_start = max(sc_5prime_len-MAX_SEARCH_LEN_5-SEARCH_OFFSET, 0)
                i_end = min(MAX_SEARCH_LEN_3+SEARCH_OFFSET, sc_3prime_len)

                if sc_5prime_len > 16+SEARCH_OFFSET+1:
                    search_seq_5 = soft_clips['fwd'][i_start:-SEARCH_OFFSET] if SEARCH_OFFSET > 0 else soft_clips['fwd'][i_start:] 
                else:
                    search_seq_5 = ''
                if sc_3prime_len > 16+SEARCH_OFFSET+1:
                    search_seq_3 = soft_clips['rev'][SEARCH_OFFSET:i_end] 
                else:
                    search_seq_3 = ''
              
                #search_seq = [search_seq_5, reverse_complement(search_seq_3)]

                if args.toolkit == '5':
                    if gene_strand == -1:
                        search_seq = reverse_complement(search_seq_3)
                    else:
                        search_seq = search_seq_5
                else:
                    if gene_strand == -1:
                        search_seq = search_seq_5
                    else:
                        search_seq = reverse_complement(search_seq_3)               

                align_strand = "plus" if gene_strand == 1 else "minus"
    
                if len(search_seq)==0:
                    continue 

                transcript_blocks_str = [str(exon[0])+','+str(exon[1]) for exon in transcript_blocks] 
                out =[samrd.qname, search_seq, gene_fname[0:-4],'|'.join(transcript_blocks_str),transcript_site, polyA_exist,align_strand]
                reads_csv.writerow(out)
            


print("{0} total reads processed; {1} gene transcripts".format(tot_rds, i))


sam_input.close()
gene_input.close()
