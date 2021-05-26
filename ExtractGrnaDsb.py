#!/usr/bin/env python
import sys,re,collections,os
from Bio import SeqIO
def Getbdg(bdg):
	d=collections.defaultdict(dict)
	for x in bdg:
		x=x.rstrip()
		l=x.split("\t")
		for i in range(int(l[1]),int(l[2])):
			d[l[0]][i]=int(l[-1])
	return d
if __name__ == '__main__':
	prefix=sys.argv[1]
	Watson=prefix+"_Watson.bdg"
	Crick=prefix+"_Crick.bdg"
	extend=30
	Watson_bed=os.popen("bedtools slop -b %s -g chrom.size -i grna_genome.bed |bedtools intersect -a %s -b - -loj -nonamecheck|awk -F'\t' '$NF!=\".\"{print $1\"\t\"$2\"\t\"$3\"\t\"$4}'|sort -u"%(extend+10,Watson)).readlines()
	Crick_bed=os.popen("bedtools slop -b %s -g chrom.size -i grna_genome.bed |bedtools intersect -a %s -b - -loj -nonamecheck|awk -F'\t' '$NF!=\".\"{print $1\"\t\"$2\"\t\"$3\"\t\"$4}'|sort -u"%(extend+10,Crick)).readlines()
	dw=Getbdg(Watson_bed)
	dc=Getbdg(Crick_bed)
	amp={}
	for record in SeqIO.parse("GRCh38.primary_assembly.genome.fa","fasta"):
		amp[record.id]=str(record.seq)
	for x in open("grna_genome.bed"):
		x=x.rstrip()
		l=x.split("\t")
		start=int(l[1])-extend
		end=int(l[2])+extend
		if l[-1]=="-":
			z=int(l[1])-1
		else:
			z=int(l[2])
		Fr=open("%s_%s_%s_dsb.xls"%(prefix,l[0],l[3]),"w")
		Fr.write("chr\tposition\trelative_position\tletter\tWatson\tCrick\n")
		for i in range(start,end):
			relative_position=i-z
			fm=l[0]+"\t"+str(i)+"\t"+str(relative_position)+"\t"+str(amp[l[0]][i])+"\t"+str(dw[l[0]][i])+"\t"+str(dc[l[0]][i])
			Fr.write(fm+"\n")
		Fr.close()
