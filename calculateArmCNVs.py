import sys
import seaborn as sns
import pandas as pd
import numpy as np
import scipy.stats
from collections import defaultdict
from collections import OrderedDict
from matplotlib import pyplot as plt

#Script for comparing Aneufinder with Copy-scAT
#Function for creating a dictionary from the Aneufinder input
def createDictionaryFromBed(bedfile):
    bed_dict=defaultdict(list)
    lines=bedfile.readlines()
    l=[]
    for line in lines:
        if line.strip()[0] != "t":
            chrom = "chr"+str(int(line.strip().split("\t")[0]))
            start = int(line.strip().split("\t")[1])
            end = int(line.strip().split("\t")[2])
            somy = int(line.strip().split("\t")[3])
            l.append([chrom, start, end, somy])

    for chrom, start, end,somy in l:
        bed_dict[chrom,start,end].append(somy)
    return(bed_dict)

#Calculate sudobulk CNV per chromosome arm from scWGS
def calculateArm(chromName, startArm, endArm, chromArm, bed_dict):
    arm_dict = defaultdict(list)
    l=[]
    for key in bed_dict:
        if chromName == key[0] and startArm < key[1] and endArm >= key[2]:
            avg_somy=sum(bed_dict[key]) / len(bed_dict[key])
            l.append([chromName, chromArm, avg_somy])
    for chromName, chromArm, avg_somy in l:
        arm_dict[chromName+chromArm].append(avg_somy)
    return(arm_dict)

def calculateCopyscAT(csv):
    averageCNV=csv.mean(axis=0)
    #averageCNV = averageCNV.to_frame().T
    averageCNV = averageCNV.to_dict()
    #print(averageCNV)
    return(averageCNV)

if __name__ =="__main__":
    fin=open("/home/katia/Helmholz/epiAneufinder/SNU_WGS/window_1e5/binsize_1e+05_stepsize_1e+05_CNV.removeCluster2.bed_CNV.converted.bed")
    bed_dict = createDictionaryFromBed(fin)
    armCoordinates=open("/home/katia/Helmholz/epiAneufinder/revisions/Copy-scAT/hg38_arms_positions.tsv")
    armHeader=armCoordinates.readline()
    armLines=armCoordinates.readlines()
    copyscat=pd.read_csv("/home/katia/Helmholz/epiAneufinder/revisions/Copy-scAT/SNU601clean_cnv_cnv_scores.csv")
    copyscatCNV=calculateCopyscAT(copyscat)
    WGSdict=defaultdict(list)
    for line in armLines:
        chromName=line.strip().split("\t")[0]
        chromArm=line.strip().split("\t")[1]
        startArm=int(line.strip().split("\t")[2])
        endArm=int(line.strip().split("\t")[3])
        arm_dict=calculateArm(chromName, startArm, endArm, chromArm,bed_dict)
        for key in arm_dict:
            WGSdict[key].append(sum(arm_dict[key])/len(arm_dict[key]))
            #print(key, sum(arm_dict[key])/len(arm_dict[key]))
    WGSdict=OrderedDict(sorted(WGSdict.items(), key=lambda t: t[0]))
    copyscatCNV=OrderedDict(sorted(copyscatCNV.items(), key=lambda t: t[0]))
    dd = defaultdict(list)
    for d in (WGSdict, copyscatCNV):  # you can list as many input dicts as you want here
        for key, value in d.items():
            dd[key].append(value)
    l0 = np.array([dd[i][0] for i in sorted(dd.keys())])
    l1=np.array([dd [i][1] for i in sorted(dd.keys())])
    print("Pearson Correlation : ", scipy.stats.pearsonr(np.hstack([dd [i][0] for i in sorted(dd.keys())]), [dd [i][1] for i in sorted(dd.keys())]))