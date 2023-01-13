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

#Function for creating a dictionary from the epiAneufinder data
def createDictionaryFromTable(table):
    snu_dict=table.set_index(['seq', 'start', 'end']).T.to_dict('list')
    return(snu_dict)

def calculatePopulationSomies(atac_dict, density_dict):
    gain_atac = []
    loss_atac = []
    base_atac = []
    common_keys = set(density_dict).intersection(atac_dict) #filtering for the common CNV locations between the two datasets
    sort_common_keys=sorted(common_keys)
    filtered_density_dict = {k: v for k, v in density_dict.items() if k in sort_common_keys}
    #print(sort_common_keys)
    counts=0
    for k in sort_common_keys:
        #if k[0]!=0: #selecting for all chromosomes
        if k[0]!=0:  # selecting for all chromosomes
            counts=counts+1
            #Calculating pseudobulk representation for the scATAC. 0 is loss, 1 is disomic and 2 is gain
            #If the user changes notation it should be changed here as well
            loss_atac.append(atac_dict[k].count(0) / len(atac_dict[k]))
            base_atac.append(atac_dict[k].count(1) / len(atac_dict[k]))
            gain_atac.append(atac_dict[k].count(2) / len(atac_dict[k]))
    print("Count Bins:",counts)
    return(loss_atac, base_atac, gain_atac, filtered_density_dict)

#Function for calculating different metrics between the two datasets and creating a line plot of the pseudoibulk data
def createLinePlot(density_dict, loss_atac, base_atac, gain_atac):
    new_base_atac = [x * 2 for x in base_atac]
    new_gain_atac = [x * 3 for x in gain_atac]
    atac_plot = [sum(x) for x in zip(new_gain_atac, new_base_atac, loss_atac)]
    atac_array=np.array(atac_plot)
    density_array=[x for x in density_dict.values()]
    x = list(range(len(atac_plot)))
    plt.plot(x,density_array)
    plt.plot(x, atac_plot, color='orange', label="ATAC")
    plt.show()
    #print(density_array)
    print("Pearson Correlation : ",scipy.stats.pearsonr(atac_array, density_array))
    print("Spearman Correlation : ", scipy.stats.spearmanr(atac_array, density_array))
    print("Kendall Correlation : ", scipy.stats.kendalltau(atac_array, density_array))

if __name__ =="__main__":
    density_table=pd.read_csv("/home/katia/Helmholz/epiAneufinder/Hg38_geneDensity.csv", sep="\t")
    snu_full=pd.read_csv("/home/katia/Helmholz/epiAneufinder/revisions/SNU601_br15/epiAneufinder_results/results_table.tsv", sep=" ")
    snu_dict=createDictionaryFromTable(snu_full)
    density_dict=createDictionaryFromTable(density_table)
    loss_atac, base_atac, gain_atac , filtered_density_dict= calculatePopulationSomies(snu_dict,density_dict)
    #print(filtered_density_dict)
    createLinePlot(filtered_density_dict, loss_atac, base_atac, gain_atac)