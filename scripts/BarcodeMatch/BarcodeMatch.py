#!/usr/bin/env python
# coding: utf-8

import os
import argparse
import numpy as np
import sys
import multiprocessing as ml
from functools import partial
import warnings
from pathos.multiprocessing import ProcessingPool as Pool

def parse_commandline():
    
    parser=argparse.ArgumentParser(description='do barcode matching between short and long reads')
    ### input files ###
    parser.add_argument('--project', '-p', help='project name', type=str, default = "Longcell_bc")
    parser.add_argument('--seq', '-q', help='input softclips', type=str, required=True)
    parser.add_argument('--barcodes', '-c', help='cellranger barcodes file', type=str, required=True)
    parser.add_argument('--cores', '-co', help='cores number for calculation', type=int, required=False, default=1)
    
    ### parameters (position) ###
    parser.add_argument('--mu', '-u', help='imputated start postions of barcodes', type=int, required=False, default=20)
    parser.add_argument('--sigma', '-s', help='start standard deviation of start postions', type=int, required=False, default=10)
    parser.add_argument('--sigma_start', '-ss', help='adjust parameter to avoid too small standard deviation during first several iterations', type=int, required=False, default=10)
    parser.add_argument('--kmer', '-k', help='k for kmer overlap between softclips and barcodes', type=int, required=False, default=8)
    parser.add_argument('--batch', '-b', help='the size of how many sequences manipulated at the same time', type=int, required=False, default=100)
    
    ### parameters (cosine scores) ###
    parser.add_argument('--top', '-t', help='barcodes with top n cos scores to use', type=int, required=False, default=8)
    parser.add_argument('--cos_thresh', '-ct', help='threshold for cos score to filter barcodes', type=int, required=False, default=0.25)
    
    ### output files ###
    parser.add_argument('--output', '-o', help='output filename', type=str, required=False,default="bc.txt")

    args=parser.parse_args()
    print(args)
    return args

def command(seq,barcodes,mu,sigma,sigma_start,kmer,batch,top,cos_thresh,output):
    os.system("\""+sys.path[0]+"/BarcodeMatch\" "+seq+" " + barcodes + " "+
    str(mu) + " " + str(sigma) + " " + str(sigma_start) + " " + str(kmer) + " " + 
    str(batch) + " " + str(top) + " " + str(cos_thresh) + " " + output)

def barcode_match(arg):
    reads_name = []
    softclips = []

    cores = arg.cores
    
    if cores == 1:
        command(arg.seq,arg.barcodes,arg.mu,arg.sigma,arg.sigma_start,
        arg.kmer,arg.batch,arg.top,arg.cos_thresh,arg.output)
            
            
    elif(cores > ml.cpu_count()):
        cores = ml.cpu_count()
        warnings.warn('only '+str(cores) + ' are available!', ResourceWarning)
        
    print("tasks will run on "+ str(cores) + " cores")
    
    if cores > 1:
        with open(arg.seq, "r") as f:
            for line in f.readlines():
                line = line.strip('\n')
                line = line.split('\t')
                reads_name.append(line[0])
                softclips.append(line[1])
        
        pool = Pool(processes=cores)
    
        read_num = len(softclips)
        batch = read_num//cores  + 1
    
        for i in range(0,cores):
            filename = arg.project + "_softclips_sub_"+str(i)+".txt"
            with open(filename,'w') as f: 
                for j in range((i)*batch,(i+1)*batch):
                    if(j >= read_num):
                        break
                    if(j == ((i+1)*batch - 1) or j == read_num - 1):
                        f.writelines(reads_name[j] + " " + softclips[j])
                    else:
                        f.writelines(reads_name[j] + " " + softclips[j] + "\n")
    
        seq = [arg.project + "_softclips_sub_"+ str(i) + ".txt" for i in range(0,cores)]
        out = [arg.project + "_output_sub_"+ str(i) +".txt" for i in range(0,cores)]
        
        pool.map(command, seq,[arg.barcodes]*cores,[arg.mu]*cores,[arg.sigma]*cores,
        [arg.sigma_start]*cores,[arg.kmer]*cores,[arg.batch]*cores,[arg.top]*cores,
        [arg.cos_thresh]*cores,out)
    
        for filename in seq:
            os.remove(filename)
        
        id = 0
        outfile=open(arg.output,'w')
        for filename in out:   
            for line in open(filename): 
                line = line.split(' ')
                line[0] = str(int(line[0]) + id*batch + 1)
                line = ' '.join(line)
                outfile.writelines(line)
                id = id + 1
            os.remove(filename)

args = parse_commandline()
barcode_match(args)





