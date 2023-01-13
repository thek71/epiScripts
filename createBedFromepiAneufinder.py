import sys
import seaborn as sns
import pandas as pd
import fargv
import numpy as np
import scipy.stats
from collections import defaultdict
from matplotlib import pyplot as plt
from sklearn.metrics import r2_score, mean_absolute_error
#plt.style.use('seaborn-whitegrid')
#sns.set_theme()

#Script for converting the table output of epiAnreufinder in pseudobulk bed format.

def createDictionaryFromTable(table):
    snu_dict=table.set_index(['seq', 'start', 'end']).T.to_dict('list')
    return(snu_dict)

def calculatePopulationSomies(atac_dict):
    gain_atac = []
    loss_atac = []
    base_atac = []
    counts=0
    for k in atac_dict.keys():
        if k[0]!=0:  # selecting for all chromosomes
            counts=counts+1
            #Calculating pseudobulk representation for the scATAC. 0 is loss, 1 is disomic and 2 is gain
            #If the user changes notation it should be changed here as well
            loss_atac.append(atac_dict[k].count(0) / len(atac_dict[k]))
            base_atac.append(atac_dict[k].count(1) / len(atac_dict[k]))
            gain_atac.append(atac_dict[k].count(2) / len(atac_dict[k]))
    print("Count Bins:",counts)
    return(atac_dict.keys(),loss_atac, base_atac, gain_atac)

def createLinePlot(loss_atac, base_atac, gain_atac):
    new_base_atac = [x * 2 for x in base_atac]
    new_gain_atac = [x * 3 for x in gain_atac]
    atac_plot = [sum(x) for x in zip(new_gain_atac, new_base_atac, loss_atac)]
    return(atac_plot)
    #print(atac_plot)

if __name__ =="__main__":
    p = {"epiTable": "/tmp/1.csv", "outBed": "/tmp/2.csv"}
    p, _ = fargv.fargv(p)
    snu_full=pd.read_csv(p.epiTable, sep=" ")
    snu_dict=createDictionaryFromTable(snu_full)
    keys,loss_atac, base_atac, gain_atac = calculatePopulationSomies(snu_dict)
    #print(loss_atac)
    atac_plot=createLinePlot(loss_atac, base_atac, gain_atac)
    bed=[x for x in zip(keys,atac_plot)]
    f=open(p.outBed,"w")
    for j in bed:
        line='\t'.join(str(s) for s in j)
        f.write(line+"\n")
    f.close()