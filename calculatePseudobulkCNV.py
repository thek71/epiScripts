import fargv
import sys
import numpy as np
import statistics

#Script for assigning the CNV status in discrete values, using the mean value per window and +2*standard deviation or higher for gain/-2*standard deviation for loss.

def calculateCNV(lines):
    chr=[]
    start=[]
    end=[]
    cnv=[]
    bulkValues=[]
    newLine=[]
    for line in lines:
        chr.append(line.strip().split("\t")[0])
        start.append(line.strip().split("\t")[1])
        end.append(line.strip().split("\t")[2])
        bulkValues.append(float(line.strip().split("\t")[3]))
    meanBulk=statistics.mean(bulkValues)
    print(meanBulk)
    sdBulk=statistics.stdev(bulkValues)
    print(sdBulk)
    for i in bulkValues:
        if i<=meanBulk-2*sdBulk:
            cnv.append("1")
        elif i>=meanBulk+2*sdBulk:
            cnv.append("3")
        else:
            cnv.append("2")
    for (a,b,c,d) in zip(chr,start,end,cnv):
        newLine.append(a+"\t"+b+"\t"+c+"\t"+d)
    return(newLine)


if __name__ == "__main__":
    p={"inputBed":"/home/katia/Helmholz/epiAneufinder/revisions/GSM4861367_COLO320HSR_rep1_atac/epiAneufinder_results/colo320HSP_rep1_results_table.bed", 
       "outputBed":"/home/katia/Helmholz/epiAneufinder/revisions/GSM4861367_COLO320HSR_rep1_atac/epiAneufinder_results/colo320HSP_rep1_results_table.CNV.bed"}
    p, _ = fargv.fargv(p)
    fin=open(p.inputBed,"r")
    lines=fin.readlines()
    newLines=calculateCNV(lines)
    with open(p.outputBed, "w") as outfile:
        outfile.write("\n".join(newLines))
