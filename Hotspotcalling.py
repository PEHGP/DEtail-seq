#!/usr/bin/env python
import pyBigWig,sys,os,argparse,pandas,collections
import numpy as np
from scipy.stats import poisson
import statsmodels.stats.multitest
def get_parser():
	parser = argparse.ArgumentParser(description="Hotspotcalling.py",formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--bw',dest="BigWig",required=True)
	parser.add_argument('--strand',dest="Strand",choices=["fwd","rev"],default="fwd")
	parser.add_argument('--qvalue',dest="Qvalue",type=float,default=0.01)
	parser.add_argument("--filter",dest="Filter",help="A list of comma-delimited chromosome names,containing those chromosomes that should be excluded for computing the hotspotcalling and normalization.",default=None)
	parser.add_argument('--res',dest="Rest",default=None,help="Restriction sites bed file")
	parser.add_argument("--norm",dest="Norm",action="store_true",help="normalizing bw file")
	parser.add_argument("--debug",dest="Debug",action="store_true")
	parser.add_argument('--prefix',dest="Prefix",required=True)
	return parser
def DealBw(Bw,Prefix,CutDict,FilterChrom,Norm=False,Rest=False,Debug=False,Scale=0):#Divide by scale
	Fr_log=open(Prefix+".log","w")
	BwHandle=pyBigWig.open(Bw)
	Header=[]
	ChrDic=BwHandle.chroms()
	for Chr in ChrDic:
		ChrLength=ChrDic[Chr]
		if Chr in FilterChrom:
			continue
		Header.append((Chr,ChrLength))
	if Norm and Rest:
		BwNorm=pyBigWig.open("%s_norm_nocut.bw"%Prefix, "w")
		if Debug:
			BwNormDebug=pyBigWig.open("%s_norm.bw"%Prefix, "w") #check (max cut site)==1?
			BwNormDebug.addHeader(Header)
	if Norm and not Rest:
		BwNorm=pyBigWig.open("%s_norm.bw"%Prefix, "w")
	if Norm:
		BwNorm.addHeader(Header)
		Fr_log.write("Scale:%s # Divide by scale\n"%Scale)
	BwValues=[]
	CutValues=[]
	for x in BwHandle.chroms().items():
		Chr=x[0]
		if Chr in FilterChrom:
			continue
		flag=0
		if CutDict:
			inl=sorted(CutDict[Chr])
		else:
			inl=[]
		#Fr.write(x[0]+"\t"+str(x[1])+"\n")
		for n in BwHandle.intervals(Chr):
			#BdgFr.write(Chr+"\t"+str(n[0])+"\t"+str(n[1])+"\t"+str(n[2])+"\n")
			if n[2]>0:
				for i in range(n[0],n[1]):
					if inl and flag!=len(inl):
						start=inl[flag][0]
						end=inl[flag][1]
						if i<start:
							BwValues.append(n[2])
							if Norm:
								BwNorm.addEntries([Chr],[i],[i+1],[n[2]/Scale])
							if Norm and Rest and Debug:
								BwNormDebug.addEntries([Chr],[i],[i+1],[n[2]/Scale])
						elif i>=end:
							flag+=1
							if flag!=len(inl):
								start=inl[flag][0]
								if i<start:
									BwValues.append(n[2])
									if Norm:
										BwNorm.addEntries([Chr],[i],[i+1],[n[2]/Scale])
									if Norm and Rest and Debug:
										BwNormDebug.addEntries([Chr],[i],[i+1],[n[2]/Scale])
								else:
									Fr_log.write(Chr+"\t"+str(i)+"\t"+str(n[2])+"\n")
									CutValues.append(n[2])
									if Norm:
										BwNorm.addEntries([Chr],[i],[i+1],[0.0])
									if Norm and Rest and Debug:
										BwNormDebug.addEntries([Chr],[i],[i+1],[n[2]/Scale])		
							else:
								BwValues.append(n[2])
								if Norm:
									BwNorm.addEntries([Chr],[i],[i+1],[n[2]/Scale])
								if Norm and Rest and Debug:
									BwNormDebug.addEntries([Chr],[i],[i+1],[n[2]/Scale])
						else:
							Fr_log.write(Chr+"\t"+str(i)+"\t"+str(n[2])+"\n")
							CutValues.append(n[2])
							if Norm:
								BwNorm.addEntries([Chr],[i],[i+1],[0.0])
							if Norm and Rest and Debug:
								BwNormDebug.addEntries([Chr],[i],[i+1],[n[2]/Scale])
					else:
						BwValues.append(n[2])
						if Norm:
							BwNorm.addEntries([Chr],[i],[i+1],[n[2]/Scale])
						if Norm and Rest and Debug:
							BwNormDebug.addEntries([Chr],[i],[i+1],[n[2]/Scale])
			else:
				if Norm:
					BwNorm.addEntries([Chr],[n[0]],[n[1]],[n[2]])
				if Norm and Rest and Debug:
					BwNormDebug.addEntries([Chr],[n[0]],[n[1]],[n[2]])

	Fr_log.close()
	BwHandle.close()
	if not Norm:
		Mean=np.mean(BwValues)
		MaxCut=max(CutValues)
		return Mean,MaxCut
	else:
		BwNorm.close()
		if Debug:
			BwNormDebug.close()
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
	RestBedFile=Args.Rest
	Norm=Args.Norm
	Debug=Args.Debug
	Prefix=Args.Prefix
	Filter=Args.Filter
	CutDict=collections.defaultdict(list)
	if RestBedFile:
		for x in open(RestBedFile):
			x=x.rstrip()
			l=x.split("\t")
			CutDict[l[0]].append((int(l[1]),int(l[2])))
		Rest=True
	if Strand=="fwd":
		Strand="+"
	elif Strand=="rev":
		Strand="-"
	FilterChrom=[]
	if Filter:
		FilterChrom=Filter.split(",")
	Mean,MaxCut=DealBw(Bw,Prefix,CutDict,FilterChrom,Norm=False,Rest=False,Debug=False,Scale=0)
	print(Mean)
	#sys.exit()
	Fr=open("%s_basepeaks.xls"%Prefix,"w")
	Fr.write("#mean %s\n"%Mean)
	Fr.write("#Threshold qvalue %s\n"%Qvalue)
	Fr.write("#chr\tstart\tend\tname\tcount\tpvalue\tqvalue\n")
	Fr_bed=open("%s_temp_basepeaks.bed"%Prefix,"w")
	pl=[]
	fl=[]
	BwHandle=pyBigWig.open(Bw)
	for x in BwHandle.chroms().items():
		Chr=x[0]
		if Chr in FilterChrom:
			continue
		for n in BwHandle.intervals(Chr):
			if n[2]>0:
				prob=poisson.cdf(float(n[2]),Mean)
				pvalue=1-prob
				for i in range(n[0],n[1]):
					fm=Chr+"\t"+str(i)+"\t"+str(i+1)
					fl.append((fm,n[2],pvalue))
					pl.append(pvalue)
		#print(Chr)
		#os.system("free -g")
	BwHandle.close()
	qv=statsmodels.stats.multitest.multipletests(pl,method='fdr_bh') #Benjamini/Hochberg (non-negative)
	fd={}
	for fm,q in zip(fl,qv[1]):
		if q <= Qvalue:
			#print fm
			fd[fm[0]]=str(fm[1])+"\t"+str(fm[2])+"\t"+str(q)
			Fr_bed.write(fm[0]+"\t.\t%s\t%s"%(fm[1],Strand)+"\n")
	Fr_bed.close()
	os.system("bedtools sort -i %s_temp_basepeaks.bed|awk -F'\t' '{n+=1;print $1\"\t\"$2\"\t\"$3\"\t%s_basepeak\"n\"\t\"$5\"\t\"$6}' >%s_basepeaks.bed"%(Prefix,Prefix,Prefix))
	for x in open("%s_basepeaks.bed"%Prefix):
		x=x.rstrip()
		l=x.split("\t")
		p="\t".join(l[:3])
		Fr.write(p+"\t"+l[3]+"\t"+fd[p]+"\n")
	Fr.close()
	#os.system("bedtools merge -i %s_basepeaks.bed -c 5 -o mean|awk -F'\t' '{n+=1;print $1\"\t\"$2\"\t\"$3\"\t%s_merge_basepeak\"n\"\t\"$4}' >%s_merge_basepeaks.bed"%(Prefix,Prefix,Prefix))
	if RestBedFile:
		os.system("bedtools intersect -nonamecheck -v -a %s_basepeaks.bed -b %s >%s_basepeaks_nocut.bed"%(Prefix,RestBedFile,Prefix))
	if Norm:
		if Rest:
			Scale=MaxCut/10000 #10000 may be need change
		else:
			Scale=Mean
		DealBw(Bw,Prefix,CutDict,FilterChrom,Norm=Norm,Rest=Rest,Debug=Debug,Scale=Scale)