#!/usr/bin/env python
import pyBigWig,sys,os,argparse,pandas,collections
import numpy as np
from scipy.stats import poisson
import statsmodels.stats.multitest
def get_parser():
	parser = argparse.ArgumentParser(description="Hotspotcalling.py",formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--bw',dest="BigWig")
	parser.add_argument('--strand',dest="Strand",choices=["fwd","rev"],default="fwd")
	parser.add_argument('--qvalue',dest="Qvalue",type=float,default=0.01)
	parser.add_argument('--res',dest="Rest",default="no",help="Restriction sites bed file")
	parser.add_argument('--prefix',dest="Prefix")
	return parser
def is_overlap(dsb_li,cut_d):
	intervals=pandas.arrays.IntervalArray.from_tuples(cut_d[dsb_li[0]])
	if True in intervals.overlaps(pandas.Interval(dsb_li[1], dsb_li[2])):
		return True
	else:
		return False
if __name__ == '__main__':
	parser=get_parser()
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit()
	else:
		Args = parser.parse_args()
	Strand=Args.Strand
	Bw=Args.BigWig
	Qvalue=Args.Qvalue
	BedFile=Args.Rest
	Prefix=Args.Prefix
	cut_d=collections.defaultdict(list)
	if BedFile!="no":
		for x in open(BedFile):
			x=x.rstrip()
			l=x.split("\t")
			cut_d[l[0]].append((int(l[1]),int(l[2])))
	if Strand=="fwd":
		Strand="+"
	elif Strand=="rev":
		Strand="-"
	Fr_log=open(Prefix+".log","w")
	BwHandle=pyBigWig.open(Bw)
	BwValues=[]
	for x in BwHandle.chroms().items():
		Chr=x[0]
		flag=0
		if cut_d:
			inl=sorted(cut_d[Chr])
		else:
			inl=[]
		#Fr.write(x[0]+"\t"+str(x[1])+"\n")
		for n in BwHandle.intervals(Chr):
			#BdgFr.write(Chr+"\t"+str(n[0])+"\t"+str(n[1])+"\t"+str(n[2])+"\n")
			for i in range(n[0],n[1]):
				if n[2]>0:
					if inl and flag!=len(inl):
						start=inl[flag][0]
						end=inl[flag][1]
						if i<start:
							BwValues.append(n[2])
						elif i>=end:
							flag+=1
							if flag!=len(inl):
								start=inl[flag][0]
								if i<start:
									BwValues.append(n[2])
								else:
									Fr_log.write(Chr+"\t"+str(i)+"\t"+str(n[2])+"\n")		
							else:
								BwValues.append(n[2])
						else:
							Fr_log.write(Chr+"\t"+str(i)+"\t"+str(n[2])+"\n")

					else:
						BwValues.append(n[2])
	Fr_log.close()
	mu=np.mean(BwValues)
	print(mu)
	#sys.exit()
	Fr=open("%s_basepeaks.xls"%Prefix,"w")
	Fr.write("#mean %s\n"%mu)
	Fr.write("#Threshold qvalue %s\n"%Qvalue)
	Fr.write("#chr\tstart\tend\tname\tcount\tpvalue\tqvalue\n")
	Fr_bed=open("%s_temp_basepeaks.bed"%Prefix,"w")
	pl=[]
	fl=[]
	for x in BwHandle.chroms().items():
		Chr=x[0]
		for n in BwHandle.intervals(Chr):
			prob=poisson.cdf(float(n[2]),mu)
			pvalue=1-prob
			for i in range(n[0],n[1]):
				fm=Chr+"\t"+str(i)+"\t"+str(i+1)
				fl.append((fm,n[2],pvalue))
				pl.append(pvalue)
	qv=statsmodels.stats.multitest.multipletests(pl,method='fdr_bh') #Benjamini/Hochberg (non-negative)
	fd={}
	for fm,q in zip(fl,qv[1]):
		if q <= Qvalue:
			#print fm
			fd[fm[0]]=str(fm[1])+"\t"+str(fm[2])+"\t"+str(q)
			Fr_bed.write(fm[0]+"\t.\t%s\t%s"%(fm[1],Strand)+"\n")
	Fr_bed.close()
	BwHandle.close()
	os.system("bedtools sort -i %s_temp_basepeaks.bed|awk -F'\t' '{n+=1;print $1\"\t\"$2\"\t\"$3\"\t%s_basepeak\"n\"\t\"$5\"\t\"$6}' >%s_basepeaks.bed"%(Prefix,Prefix,Prefix))
	for x in open("%s_basepeaks.bed"%Prefix):
		x=x.rstrip()
		l=x.split("\t")
		p="\t".join(l[:3])
		Fr.write(p+"\t"+l[3]+"\t"+fd[p]+"\n")
	Fr.close()
	#os.system("bedtools merge -i %s_basepeaks.bed -c 5 -o mean|awk -F'\t' '{n+=1;print $1\"\t\"$2\"\t\"$3\"\t%s_merge_basepeak\"n\"\t\"$4}' >%s_merge_basepeaks.bed"%(Prefix,Prefix,Prefix))
	if BedFile!="no":
		os.system("bedtools intersect -v -a %s_basepeaks.bed -b %s >%s_basepeaks_nocut.bed"%(Prefix,BedFile,Prefix))
