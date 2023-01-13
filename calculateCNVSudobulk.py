import sys
import seaborn as sns
import pandas as pd
import numpy as np
import scipy.stats
from collections import defaultdict
from matplotlib import pyplot as plt
from sklearn.metrics import r2_score, mean_absolute_error
#plt.style.use('seaborn-whitegrid')
#sns.set_theme()


#Script for comparing the scATAC with the scWGW dataset
#Inputs are the bed file containing the results from Aneufinder software for calling CNVs from scWGS data and the table with the results from epiAneufinder.
#The bed file from Aneufinder has been pre-processed so that the chr-star-end columns correspond in size with the epiAneufinder results.

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

#Function for creating a dictionary from the Aneufinder input
def createDictionaryFromBed(bedfile):
    bed_dict=defaultdict(list)
    lines=bedfile.readlines()
    l=[]
    for line in lines:
        if line.strip()[0] != "t":
            #chrom = "chr"+str(int(line.strip().split("\t")[0]))
            chrom=int(line.strip().split("\t")[0])
            start = int(line.strip().split("\t")[1])
            end = int(line.strip().split("\t")[2])
            somy = int(line.strip().split("\t")[3])
            l.append([chrom, start, end, somy])
            #d['chr1', 100000, 200000][0:10].count(1)

    for chrom, start, end,somy in l:
        bed_dict[chrom,start,end].append(somy)
        #print(bed_dict)
    return(bed_dict)

#Function for filtering the Aneufinder dictionary in the same way as the epiAneufinder results.
#Bins having more than 85% of the cells with zero status are removed
def filterDictionary(dict):
    new_dict = {k: v for k, v in bed_dict.items() if (bed_dict[k].count(0)/len(bed_dict[k]))<0.85}
    return(new_dict)

#Function for creating a dictionary from the epiAneufinder data
def createDictionaryFromTable(table):
    snu_dict=table.set_index(['seq', 'start', 'end']).T.to_dict('list')
    #print(snu_dict)
    return(snu_dict)

#Function for creating pseudo-bulk for both scATAC and scWGS
def calculatePopulationSomies(atac_dict,wgs_dict):
    gain_wgs=[]
    loss_wgs = []
    base_wgs = []
    gain_atac = []
    loss_atac = []
    base_atac = []
    common_keys = set(wgs_dict).intersection(atac_dict) #filtering for the common CNV locations between the two datasets
    sort_common_keys=sorted(common_keys)
    print(sort_common_keys)
    counts=0
    for k in sort_common_keys:
        #if k[0]!=0: #selecting for all chromosomes
        if k[0]!=0:  # selecting for all chromosomes
            counts=counts+1
            #Calculating pseudobulk representation for the scWGS. 1 is loss, 2 is disomic and 3 is gain
            loss_wgs.append((wgs_dict[k].count(1)+wgs_dict[k].count(0))/len(wgs_dict[k]))
            base_wgs.append(wgs_dict[k].count(2) / len(wgs_dict[k]))
            gain_wgs.append(wgs_dict[k].count(3) / len(wgs_dict[k]))
            #Calculating pseudobulk representation for the scATAC. 0 is loss, 1 is disomic and 2 is gain
            #If the user changes notation it should be changed here as well
            loss_atac.append(atac_dict[k].count(0) / len(atac_dict[k]))
            base_atac.append(atac_dict[k].count(1) / len(atac_dict[k]))
            gain_atac.append(atac_dict[k].count(2) / len(atac_dict[k]))
    print("Count Bins:",counts)
    return(loss_wgs,base_wgs, gain_wgs, loss_atac, base_atac, gain_atac)

#Function for calculating different metrics between the two datasets and creating a line plot of the pseudoibulk data
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
    #borders_1e6=[0,218,450,643,823,996,1159,1311,1449,1558,1684,1809,1935,2027,2113,2189,2264,2339,2410,2463,2520,2552,2583]
    #borders_101=[0,2284,4672,6645,8528,10362,12021,13599,15037,16229,17545,18877,20197,21167,22063,22897,23705,24521,25304,25883,26495,26864,27233]
    #borders_10000_85=[0,2230,4592,6541,8406,10161,11844,13386,14802,15921,17221,18525,19826,20779,21654,22441,23207,23980,24719,25274,25875,26212,26551]
    plt.plot(x, wgs_plot, color='#98d1d1', label="GS")
    plt.plot(x, atac_plot, color='#df979e', label="ATAC")
    for border in borders_hct:
        plt.axvline(border, color='gray')
    plt.title("HCT116 scATAC br15 minsizeCNV=0 compared to scDNA")
    plt.xlim((0, len(atac_plot)))
    plt.legend()
    plt.show()

if __name__ =="__main__":
    fin=open("/home/katia/Helmholz/epiAneufinder/HCT116/WGS/HCT116_WT_T0_2N/BROWSERFILES/method-edivisive/binsize_1e+05_stepsize_1e+05_CNV.converted.bed")
    snu_full=pd.read_csv("/home/katia/Helmholz/epiAneufinder/revisions/HCT_br15_msCNV0/epiAneufinder_results/results_table.tsv", sep=" ")
    #sample1=pd.read_csv("/home/katia/Helmholz/epiAneufinder/revisions/GSM4861381_COLO320HSR_rep8_atac/epiAneufinder_results/colo320HSP_rep8_results_table.tsv", sep=" ")
    #sample1_dict=createDictionaryFromTable(sample1)
    #sample2 = pd.read_csv("/home/katia/Helmholz/epiAneufinder/revisions/GSM4861379_COLO320HSR_rep7_atac/epiAneufinder_results/colo320HSP_rep7_results_table.tsv", sep=" ")
    snu_dict = createDictionaryFromTable(snu_full)
    bed_dict=createDictionaryFromBed(fin)
    filtered_dict=filterDictionary(bed_dict)
    loss_wgs, base_wgs, gain_wgs, loss_atac, base_atac, gain_atac = calculatePopulationSomies(snu_dict,bed_dict)
    createLinePlot(loss_wgs, base_wgs, gain_wgs, loss_atac, base_atac, gain_atac)
    #print(gain_atac)
    #print(loss_wgs,loss_atac)
