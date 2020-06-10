#!/usr/bin/env python
import pyBigWig,sys,os
import numpy as np
from scipy.stats import poisson
import statsmodels.stats.multitest
Bw=sys.argv[1]
Strand=sys.argv[2]
Prefix=sys.argv[3]
Qvalue=float(sys.argv[4])
if Strand=="fwd":
	Strand="+"
elif Strand=="rev":
	Strand="-"
else:
	print("strand error")
	sys.exit()
BwHandle=pyBigWig.open(Bw)
BwValues=[]
for x in BwHandle.chroms().items():
	Chr=x[0]
	#Fr.write(x[0]+"\t"+str(x[1])+"\n")
	for n in BwHandle.intervals(Chr):
		#BdgFr.write(Chr+"\t"+str(n[0])+"\t"+str(n[1])+"\t"+str(n[2])+"\n")
		for i in range(n[0],n[1]):
			if n[2]>0:
				BwValues.append(n[2])
#Fr.close()
#BdgFr.close()
#BwHandle.close()
mu=np.mean(BwValues)
print(mu)
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
		fm=Chr+"\t"+str(n[0])+"\t"+str(n[1])
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
os.system("bedtools merge -i %s_basepeaks.bed -c 5 -o mean|awk -F'\t' '{n+=1;print $1\"\t\"$2\"\t\"$3\"\t%s_merge_basepeak\"n\"\t\"$4}' >%s_merge_basepeaks.bed"%(Prefix,Prefix,Prefix))
