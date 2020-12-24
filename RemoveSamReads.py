#!/usr/bin/env python
from Bio import SeqIO
#import itertools
import os,sys, argparse
import gzip

def ParseArg():
    p=argparse.ArgumentParser(description = 'Remove duplicated reads which have same sequences for both forward and reverse reads. Choose the one appears first.')
    p.add_argument('input1',type=str,metavar='reads1',help='left gz')
    p.add_argument('input2',type=str,metavar='reads2',help='right gz')
    p.add_argument('prefix',type=str,metavar='prefix')
    if len(sys.argv)==1:
        p.print_help()
        exit(0)
    return p.parse_args()

def Main():
    Unique_seqs=set()
    args=ParseArg()
    outfile1 = gzip.open(args.prefix+"_dup_R1.fq.gz","wt")
    outfile2 = gzip.open(args.prefix+"_dup_R2.fq.gz","wt")
    fastq_iter1 = SeqIO.parse(gzip.open(args.input1,"rt"),"fastq")
    fastq_iter2 = SeqIO.parse(gzip.open(args.input2,"rt"),"fastq")
    for rec1, rec2 in zip(fastq_iter1, fastq_iter2):
        if str((rec1+rec2).seq) not in Unique_seqs:
            SeqIO.write(rec1,outfile1,"fastq")
            SeqIO.write(rec2,outfile2,"fastq")
            Unique_seqs.add(str((rec1+rec2).seq))
    outfile1.close()
    outfile2.close()

if __name__ == '__main__':
    Main()
