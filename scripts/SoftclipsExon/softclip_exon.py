#!/usr/bin/env python


import argparse, sys, os, re, pysam, csv, gzip, string,json
import numpy as np, pandas as pd



DNA_COMPLEMENT = str.maketrans("ACGT", "TGCA") 

SEARCH_OFFSET = 0
MAX_SEARCH_LEN = 55-SEARCH_OFFSET  
polyA= "A"*10

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
  
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Define reusable modules to identify exons for each reads    #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#  

def createExonList(gene_input):
    exon_rds = list(csv.reader(gene_input, delimiter=" "))
    return exon_rds 

def getTranscriptSpan(gene_exons, min_overlap, exon_span):
    gene_strand = 1 if gene_exons[-1][3] > gene_exons[0][3] else -1
    if exon_span == 'all':
        max_start = int(gene_exons[-1][1]) - min_overlap
        min_end   = int(gene_exons[0][2]) + min_overlap
    else:
        last_exon = gene_exons[-1] if gene_strand == 1 else gene_exons[0]
        max_start = int(last_exon[2]) - min_overlap
        min_end   = int(last_exon[1]) + min_overlap
    return [gene_strand, max_start, min_end]
  
def isSpanningRead(exon_rds, gene_strand, exon_span):
    if exon_span == 'all':
        spanning_read = (exon_rds[0] > 0 and exon_rds[-1] > 0)
    elif gene_strand == 1:
        spanning_read = exon_rds[-1] > 0
    else:
        spanning_read = exon_rds[0] > 1
    return spanning_read
 
# function to find first index >= x 
def lowerIndex(arr, n, x): 
    l = 0
    h = n-1
    while (l <= h): 
        mid = int((l + h)/2) 
        if (arr[mid] >= int(x)): 
            h = mid - 1
        else: 
            l = mid + 1
    return l 
  
  
# function to find last index <= x 
def upperIndex(arr, n, x): 
    l = 0
    h = n-1
    while (l <= h): 
        mid = int((l + h)/2) 
        if (arr[mid] <= int(x)): 
            l = mid + 1
        else: 
            h = mid - 1
    return h 
  
# define paste0 function to write with the delimiters #
def paste0(string1, string2):
    texts = [string1 + string2]
    return texts 



# function to count elements within given range 
def countInRange(arr, n, x, y): 
  # initialize result 
    count = 0; 
    count = upperIndex(arr, n, y) - lowerIndex(arr, n, x) + 1; 
    return count

def printDebug(i, qname, chkrdnames):
    return True if (i <20 or qname in chkrdnames) else False; 



script_name = os.path.basename(__file__)
print("Running ", script_name)

MIN_OVERLAP = 12  
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
exon_span = "all"
transcript_info = getTranscriptSpan(gene_exons, MIN_OVERLAP, exon_span)
gene_strand = transcript_info[0]
transcript_max_start = transcript_info[1]
transcript_min_end   = transcript_info[2]

for samrd in sam_input.fetch(gene_chr,gene_start, gene_end): 
    tot_rds += 1
    if samrd.is_unmapped or samrd.is_secondary:
        continue
    else:
        chrnm = sam_input.getrname(samrd.tid)
        if chrnm == gene_chr:
            if samrd.qname in CHKRD_NAMES or (samrd.reference_start <= transcript_max_start and samrd.reference_end > transcript_min_end):
                i += 1
                if printDebug(i, samrd.qname, CHKRD_NAMES):
                    print("{2} Full transcript start: {0}, end: {1}.  Alignment per exon below".format(samrd.reference_start, samrd.reference_end, samrd.qname))
                ### softclips extraction    
                soft_clips = extract_softclips(samrd)
                 
                sc_5prime_len = len(soft_clips['fwd'])
                sc_3prime_len = len(soft_clips['rev'])
                if gene_strand == -1:
                    polyA_exist = polyA in reverse_complement(soft_clips['fwd'])
                else:
                    polyA_exist = polyA in soft_clips['rev']

                i_start = max(sc_5prime_len-MAX_SEARCH_LEN-SEARCH_OFFSET, 0)
                i_end = min(MAX_SEARCH_LEN+SEARCH_OFFSET, sc_3prime_len) 

                if sc_5prime_len > 16+SEARCH_OFFSET+1:   
                    if sc_3prime_len > 16+SEARCH_OFFSET+1:                            
                        search_seq_5 = soft_clips['fwd'][i_start:-SEARCH_OFFSET] if SEARCH_OFFSET > 0 else soft_clips['fwd'][i_start:] 
                        search_seq_3 = soft_clips['rev'][SEARCH_OFFSET:i_end] 
                    else:
                        search_seq_5 = soft_clips['fwd'][i_start:-SEARCH_OFFSET] if SEARCH_OFFSET > 0 else soft_clips['fwd'][i_start:] 
                        search_seq_3 = ''
                else:
                    search_seq_5 = ''
                    if sc_3prime_len > 16+SEARCH_OFFSET+1:  
                        search_seq_3 = soft_clips['rev'][SEARCH_OFFSET:i_end]
                    else:
                        search_seq_3 = ''
                if args.toolkit == '5':
                    if gene_strand == -1:
                        search_seq = reverse_complement(search_seq_3)
                        align_strand = "minus"
                    else:
                        search_seq = search_seq_5
                        align_strand = "plus"
                else:
                    if gene_strand == -1:
                        search_seq = reverse_complement(search_seq_5)
                        align_strand = "minus"
                    else:
                        search_seq = search_seq_3
                        align_strand = "plus"
                ###exon identification
                transcript_rds = [i + 1 for i in samrd.get_reference_positions()]
                exon_rds = np.zeros((len(gene_exons),), dtype=float)
                for j, exon in enumerate(gene_exons):
                    nr_exon_bases = countInRange(transcript_rds, len(transcript_rds), exon[1], exon[2])
                    exon_rds[j] = nr_exon_bases/(int(exon[2]) - int(exon[1]) + 1)
                    if printDebug(i, samrd.qname, CHKRD_NAMES):
                        print("Exon {0}: {1} aligned bases".format(exon[3], nr_exon_bases))

                if printDebug(i, samrd.qname, CHKRD_NAMES):
                    print(exon_rds)
                
                exon_rds_str = [str(exon)+'|' for exon in exon_rds]

                exon_rds_bin = np.where(exon_rds > 0, 1, 0)

                if sum(exon_rds_bin)!=0:
                    exon_rds_bin = np.where(exon_rds >= 0.9, 1, 0)
                    temp = np.where(exon_rds > 0)[0]
                    if exon_rds[min(temp)] >= 0.1:
                        exon_rds_bin[min(temp)] = 1
                    if exon_rds[max(temp)] >= 0.1:
                        exon_rds_bin[max(temp)] = 1
                    exon_skip_str = [str(exon[3])+'|' for exon in gene_exons]
                    for k, num in enumerate(exon_rds_bin):
                        if num == 0: exon_skip_str[k] = '-'
                    #reads_csv.writerow([samrd.qname, ''.join(exon_skip_str),gene_fname[0:-4],search_seq,align_strand])
                    out =[samrd.qname, ''.join(exon_skip_str),gene_fname[0:-4],search_seq,polyA_exist,
                          align_strand,''.join(exon_rds_str),samrd.reference_start,samrd.reference_end]
                    #out = [samrd.qname,search_seq_5,search_seq_3,polyA_exist,align_strand]
                    reads_csv.writerow(out)
            


print("{0} total reads processed; {1} gene transcripts".format(tot_rds, i))


sam_input.close()
gene_input.close()
