import sys
import seaborn as sns
import pandas as pd
import numpy as np
import scipy.stats
from collections import defaultdict
from matplotlib import pyplot as plt
from sklearn.metrics import r2_score, mean_absolute_error

#Function for creating a dictionary from the epiAneufinder data
def createDictionaryFromTable(table):
    snu_dict=table.set_index(['seq', 'start', 'end']).T.to_dict('list')
    #print(snu_dict)
    return(snu_dict)

if __name__ =="__main__":
    atac=pd.read_csv("/home/katia/Helmholz/epiAneufinder/revisions/GSM4861367_COLO320HSR_rep1_atac_v3blacklist/epiAneufinder_results/results_table.tsv", sep=" ")
    atac_dict=createDictionaryFromTable(atac)
    for k in atac_dict:
        disomic=atac_dict[k].count(1)
        gain=atac_dict[k].count(2)
        loss=atac_dict[k].count(0)
        total=disomic+gain+loss
        percent_base=100*(disomic/total)
        percent_loss = 100 * (loss / total)
        percent_gain = 100 * (gain / total)
        print(k,percent_loss,percent_base,percent_gain)