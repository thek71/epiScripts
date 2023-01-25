import sys
import seaborn as sns
import pandas as pd
import numpy as np
import copy
import scipy.stats
from collections import defaultdict
from matplotlib import pyplot as plt
from sklearn.metrics import r2_score, mean_absolute_error
#plt.style.use('seaborn-whitegrid')
#sns.set_theme()


#Script for calculating the correlation between a normal sample and an ideal diploid sample.
#Input is the table with the results from epiAneufinder.

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

#Function for creating a dictionary from the epiAneufinder data
def createDictionaryFromTable(table):
    snu_dict=table.set_index(['seq', 'start', 'end']).T.to_dict('list')
    #print(snu_dict)
    return(snu_dict)

#Function for creating the ideal normal sample, by copying the dictionary created from the epiAneufinder data
def createNormalBase(dictionary):
    base_dict=copy.deepcopy(dictionary)
    for i, j in base_dict.items():  # use iteritems in py2k
        if j == 2:
            base_dict[i] =-2
    return(base_dict)

def calculatePopulationSomies(atac_dict,base_dict):
    gain_wgs=[]
    loss_wgs = []
    base_wgs = []
    gain_atac = []
    loss_atac = []
    base_atac = []
    common_keys = set(base_dict).intersection(atac_dict) #filtering for the common CNV locations between the two datasets
    sort_common_keys=sorted(common_keys)
    #print(sort_common_keys)
    counts=0
    for k in sort_common_keys:
        if k[0]!=0: #selecting for all chromosomes
        #if k[0]!=0:  # selecting for all chromosomes
            counts=counts+1
            #Calculating pseudobulk representation for the scWGS. 1 is loss, 2 is disomic and 3 is gain
            loss_wgs.append((base_dict[k].count(1)+base_dict[k].count(0))/len(base_dict[k]))
            base_wgs.append(base_dict[k].count(2) / len(base_dict[k]))
            gain_wgs.append(base_dict[k].count(3) / len(base_dict[k]))
            #Calculating pseudobulk representation for the scATAC. 0 is loss, 1 is disomic and 2 is gain
            #If the user changes notation it should be changed here as well
            loss_atac.append(atac_dict[k].count(0) / len(atac_dict[k]))
            base_atac.append(atac_dict[k].count(1) / len(atac_dict[k]))
            gain_atac.append(atac_dict[k].count(2) / len(atac_dict[k]))
    print("Count Bins:",counts)
    return(loss_wgs,base_wgs, gain_wgs, loss_atac, base_atac, gain_atac)

def createLinePlot(loss_wgs, base_wgs, gain_wgs, loss_atac, base_atac, gain_atac):
    new_base_wgs = [x * 2 for x in base_wgs]
    new_base_atac = [x * 2 for x in base_atac]
    new_gain_wgs = [x * 3 for x in gain_wgs]
    new_gain_atac = [x * 3 for x in gain_atac]
    wgs_plot=[sum(x) for x in zip(new_gain_wgs,new_base_wgs,loss_wgs)]
    atac_plot = [sum(x) for x in zip(new_gain_atac, new_base_atac, loss_atac)]
    atac_array=np.array(atac_plot)
    wgs_array=np.array(wgs_plot)
    #outf=open("genome.csv","w")
    #both = np.concatenate([atac_array[:, None], wgs_array[:, None]], axis=1)
    #np.savetxt(outf, both, delimiter=",")
    #outf.close()
    #print(np.corrcoef(atac_array,wgs_array))
    print("Pearson Correlation : ",scipy.stats.pearsonr(atac_array, wgs_array))
    print("Spearman Correlation : ", scipy.stats.spearmanr(atac_array, wgs_array)[0])
    print("Kendall Correlation : ", scipy.stats.kendalltau(atac_array, wgs_array)[0])

    difference_array = np.subtract(atac_array, wgs_array)
    squared_array = np.square(difference_array)
    mse = squared_array.mean()
    print("Mean Square Error: ",squared_array.mean())

    print("R square: ",r2_score(atac_array, wgs_array))
    print("Mean Absolute Error: ",mean_absolute_error(atac_array,wgs_array))
    x = list(range(len(wgs_plot)))
    borders_hct=[0,2232,4601,6553,8420,10174,11815,13381,14796,15979,17283,18597,19894,20841,21688,22475,23251,24036,24811,25349,25960,26307,26654]
    borders_snu=[0, 2231, 4593, 6542, 8406, 10161, 11844, 13386, 14802, 15922, 17224, 18529, 19830, 20783, 21658, 22445, 23211, 23985, 24725, 25280, 25881, 26218, 26557]
    #borders_1e6=[0,218,450,643,823,996,1159,1311,1449,1558,1684,1809,1935,2027,2113,2189,2264,2339,2410,2463,2520,2552,2583]
    #borders_101=[0,2284,4672,6645,8528,10362,0.86837989469864670.868379894698646712021,13599,15037,16229,17545,18877,20197,21167,22063,22897,23705,24521,25304,25883,26495,26864,27233]
    #borders_10000_85=[0,2230,4592,6541,8406,10161,11844,13386,14802,15921,17221,18525,19826,20779,21654,22441,23207,23980,24719,25274,25875,26212,26551]
    plt.plot(x, wgs_plot, color='#98d1d1', label="Base")
    plt.plot(x, atac_plot, color='#df979e', label="ATAC")
    for border in borders_snu:
        plt.axvline(border, color='gray')
    plt.title("SNU601 scATAC br15 minsizeCNV=0 compared to scDNA")
    plt.xlim((0, len(atac_plot)))
    plt.legend()
    plt.show()

if __name__ =="__main__":
    snu_full=pd.read_csv("/home/katia/Helmholz/epiAneufinder/revisions/Satpathy_BoneMarrow_msCNV0/epiAneufinder_results/results_table.tsv", sep=" ")
    snu_dict = createDictionaryFromTable(snu_full)
    base_dict = createNormalBase(snu_dict)
    print(base_dict)
    loss_wgs, base_wgs, gain_wgs, loss_atac, base_atac, gain_atac = calculatePopulationSomies(snu_dict,base_dict)
    createLinePlot(loss_wgs, base_wgs, gain_wgs, loss_atac, base_atac, gain_atac)