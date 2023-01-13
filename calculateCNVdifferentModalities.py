import sys
import seaborn as sns
import pandas as pd
import numpy as np
import scipy.stats
import fargv
from collections import defaultdict
from matplotlib import pyplot as plt
from sklearn.metrics import r2_score, mean_absolute_error


# plt.style.use('seaborn-whitegrid')
# sns.set_theme()


# Script for comparing the scATAC from multiome and sigleome
def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

#Function for creating a dictionary from the epiAneufinder data
def createDictionaryFromTable(table):
    snu_dict=table.set_index(['seq', 'start', 'end']).T.to_dict('list')
    return(snu_dict)

#Function for creating pseudo-bulk for both scATAC datasets
def calculatePopulationSomies(atac_dict1,atac_dict2):
    gain_atac1=[]
    loss_atac1 = []
    base_atac1 = []
    gain_atac2 = []
    loss_atac2 = []
    base_atac2 = []
    common_keys = set(atac_dict1).intersection(atac_dict2) #filtering for the common CNV locations between the two datasets
    sort_common_keys=sorted(common_keys)
    #print(sort_common_keys)
    common_size=len(common_keys)
    counts=0
    for k in sort_common_keys:
        #if k[0]!=0: #selecting for all chromosomes
        if k[0]!=0:  # selecting for all chromosomes
            counts=counts+1
            #Calculating pseudobulk representation for atac_dict1.  0 is loss, 1 is disomic and 2 is gain
            loss_atac1.append(atac_dict1[k].count(0)/len(atac_dict1[k]))
            base_atac1.append(atac_dict1[k].count(1) / len(atac_dict1[k]))
            gain_atac1.append(atac_dict1[k].count(2) / len(atac_dict1[k]))
            #Calculating pseudobulk representation for the scATAC. 0 is loss, 1 is disomic and 2 is gain
            #If the user changes notation it should be changed here as well
            loss_atac2.append(atac_dict2[k].count(0) / len(atac_dict2[k]))
            base_atac2.append(atac_dict2[k].count(1) / len(atac_dict2[k]))
            gain_atac2.append(atac_dict2[k].count(2) / len(atac_dict2[k]))
    print("Count Bins:",counts)
    return(loss_atac1,base_atac1, gain_atac1, loss_atac2, base_atac2, gain_atac2, common_size)

#Function for calculating different metrics between the two datasets and creating a line plot of the pseudoibulk data
def createLinePlot(loss_atac1,base_atac1, gain_atac1, loss_atac2, base_atac2, gain_atac2):
    new_base_atac1 = [x * 2 for x in base_atac1]
    new_base_atac2 = [x * 2 for x in base_atac2]
    new_gain_atac1 = [x * 3 for x in gain_atac1]
    new_gain_atac2 = [x * 3 for x in gain_atac2]
    atac1_plot=[sum(x) for x in zip(new_gain_atac1,new_base_atac1,loss_atac1)]
    atac2_plot = [sum(x) for x in zip(new_gain_atac2, new_base_atac2, loss_atac2)]
    atac2_array=np.array(atac2_plot)
    atac1_array=np.array(atac1_plot)
    #outf=open("genome.csv","w")
    #both = np.concatenate([atac_array[:, None], wgs_array[:, None]], axis=1)
    #np.savetxt(outf, both, delimiter=",")
    #outf.close()
    #print(np.corrcoef(atac_array,wgs_array))
    print("Pearson Correlation : ",scipy.stats.pearsonr(atac2_array, atac1_array))
    print("Spearman Correlation : ", scipy.stats.spearmanr(atac2_array, atac1_array)[0])
    print("Kendall Correlation : ", scipy.stats.kendalltau(atac2_array, atac1_array)[0])

if __name__ =="__main__":
    p = {"multiome_cnv": "/tmp/1.csv", "singleome_cnv": "/tmp/2.csv"}
    p, _ = fargv.fargv(p)
    atac1_set=pd.read_csv(p.multiome_cnv, sep=" ")
    atac2_set=pd.read_csv(p.singleome_cnv, sep=" ")
    atac_dict1=createDictionaryFromTable(atac1_set)
    atac_dict2 = createDictionaryFromTable(atac2_set)
    loss_atac1,base_atac1, gain_atac1, loss_atac2, base_atac2, gain_atac2, common_size = calculatePopulationSomies(atac_dict1,atac_dict2)
    print("Common bins:", str(common_size))
    createLinePlot(loss_atac1,base_atac1, gain_atac1, loss_atac2, base_atac2, gain_atac2)