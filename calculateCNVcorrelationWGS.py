import sys
import seaborn as sns
import pandas as pd
import numpy as np
import scipy.stats
import fargv
from collections import defaultdict
from matplotlib import pyplot as plt
from sklearn.metrics import r2_score, mean_absolute_error, mutual_info_score, normalized_mutual_info_score
#plt.style.use('seaborn-whitegrid')
#sns.set_theme()


#Script for comparing the scATAC with the WGS dataset
#Inputs are the bed file containing the results from cnvkit software for calling CNVs from WGS data and the table with the results from epiAneufinder.
#The input file from CNVkit is the intersect between the resulting calls and the epiAneufinder windows.

standarize = lambda x: ((x - x.mean()) / (x.std() + .0000000000000001)) #The standarization formula for bringing the two datasets in the same scale
standarize2 = lambda x: ((x / (x.std() + .0000000000000001)))
normalize = lambda x: ((2*(x-x.min())/(x.max()-x.min()))+1)


def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

#Function for creating a dictionary from the WGS input
def createDictionaryFromBedCNVkit(bedfile):
    bed_dict=defaultdict(list)
    lines=bedfile.readlines()
    l=[]
    for line in lines:
        if line.strip()[0] != "t":
            #chrom = "chr"+str(int(line.strip().split("\t")[0]))
            if line.strip().split("\t")[6]!=".":
                chrom=int((line.strip().split("\t")[6]))
                start = int(line.strip().split("\t")[7])+1
                end = int(line.strip().split("\t")[8])
                somy = int(line.strip().split("\t")[5])
            if somy<2:
                new_somy=1
            elif somy==2:
                new_somy=2
            else:
                new_somy=3
            l.append([chrom, start, end, new_somy])
            #d['chr1', 100000, 200000][0:10].count(1)

    for chrom, start, end,somy in l:
        bed_dict[chrom,start,end].append(somy)
        #print(bed_dict)
    #print(len(bed_dict))
    return(bed_dict)

#create dictionary from aneufinder bed file (created using convertAneufinderToepiBed.py)
def createDictionaryFromBed(bedfile):
    bed_dict=defaultdict(list)
    lines=bedfile.readlines()
    l=[]
    for line in lines:
        #print(line)
        if line.strip()[0] != "t":
            #chrom = "chr"+str(int(line.strip().split("\t")[0]))
            chrom=int(line.strip().split("\t")[0])
            start = int(line.strip().split("\t")[1])
            end = int(line.strip().split("\t")[2])
            somy = float(line.strip().split("\t")[3])
            if somy >=20:
                new_somy=20
            else:
                new_somy=somy

            l.append([chrom, start, end, new_somy])
            #d['chr1', 100000, 200000][0:10].count(1)

    for chrom, start, end,somy in l:
        bed_dict[chrom,start,end].append(somy)
        #print(bed_dict)
    #print((bed_dict))
    return(bed_dict)



#Function for creating a dictionary from the epiAneufinder data
def createDictionaryFromTable(table):
    snu_dict=table.set_index(['seq', 'start', 'end']).T.to_dict('list')
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
    #print(sort_common_keys)
    counts=0
    for k in sort_common_keys:
        if k[0]==9: #selecting chromosome
        #if k[0]!=0:  # selecting for all chromosomes
            counts=counts+1
            #print(wgs_dict[k])
            #Calculating pseudobulk representation for the WGS. 1 is loss, 2 is disomic and 3 is gain
            loss_wgs.append(wgs_dict[k].count(1)/len(wgs_dict[k]))
            base_wgs.append(wgs_dict[k].count(2) / len(wgs_dict[k]))
            gain_wgs.append(wgs_dict[k].count(3) / len(wgs_dict[k]))
            #Calculating pseudobulk representation for the scATAC. 0 is loss, 1 is disomic and 2 is gain
            #If the user changes notation it should be changed here as well
            loss_atac.append(atac_dict[k].count(0) / len(atac_dict[k]))
            base_atac.append(atac_dict[k].count(1) / len(atac_dict[k]))
            gain_atac.append(atac_dict[k].count(2) / len(atac_dict[k]))
    print("Count Bins:",counts)
    #print(len(loss_atac), len(base_atac), len(gain_atac))
    return(loss_wgs,base_wgs, gain_wgs, loss_atac, base_atac, gain_atac)

def calculatePopulationSomiesWGS(atac_dict,wgs_dict):
    #new_atac_dict = {k: [sum(atac_dict[k])/len(atac_dict[k])] for k in atac_dict.keys() if k[0]==4}
    new_atac_dict = {k: [sum(atac_dict[k]) / len(atac_dict[k])] for k in atac_dict.keys()}
    common_keys = set(wgs_dict).intersection(new_atac_dict) #filtering for the common CNV locations between the two datasets
    sort_common_keys=sorted(common_keys)
    #print(sort_common_keys)
    common_wgs_dict = {key: wgs_dict[key] for key in sort_common_keys}
    common_atac_dict = {key: new_atac_dict[key] for key in sort_common_keys}
    return(common_wgs_dict, common_atac_dict)

#Function for calculating different metrics between the two datasets and creating a line plot of the pseudoibulk data
def createLinePlotAneufinder(common_wgs_dict, common_atac_dict,gaussian_sigma=0,filter_edges="constant"):
    wgs_list=([i for i in common_wgs_dict.values()])
    atac_list= [i for i in common_atac_dict.values()]
    wgs_array = np.hstack(wgs_list)
    smoothed_wgs_array = scipy.ndimage.gaussian_filter1d(wgs_array, gaussian_sigma, axis=- 1, order=0, output=None,
                                                   mode=filter_edges)
    atac_array = np.hstack(atac_list)
    #atac_array=np.asarray(tuple(common_atac_dict.values()))
    #wgs_array = np.asarray(tuple(common_wgs_dict.values()))
    print("MIS:", normalized_mutual_info_score(atac_array,wgs_array, average_method='min'))
    print("Pearson Correlation : ", scipy.stats.pearsonr(standarize(atac_array),standarize(smoothed_wgs_array)))
    print("Spearman Correlation : ", scipy.stats.spearmanr(atac_array, smoothed_wgs_array)[0])
    print("Kendall Correlation : ", scipy.stats.kendalltau(atac_array, smoothed_wgs_array)[0])

    difference_array = np.subtract(atac_array, smoothed_wgs_array)
    squared_array = np.square(difference_array)
    mse = squared_array.mean()
    print("Mean Square Error: ",squared_array.mean())

    print("R square: ",r2_score(atac_array, smoothed_wgs_array))
    print("Mean Absolute Error: ",mean_absolute_error(atac_array,smoothed_wgs_array))
    x = list(range(len(smoothed_wgs_array)))
    plt.plot(x,(standarize(atac_array)), color='#df979e', label="ATAC")
    plt.plot(x,(standarize(smoothed_wgs_array)+0.4), color='#98d1d1', label="WGS")
    borders_colorep1 = [0, 2109, 4392, 6241, 8068,9661, 11085, 12484, 13785, 14854, 16070, 17287, 18527, 19236, 20058, 20758, 21504, 22212, 22836, 23312, 23867, 24172, 24482]
    #plt.plot(x,((atac_array)), color='#df979e', label="ATAC")
    #plt.plot(x,((wgs_array)), color='#98d1d1', label="WGS")
    plt.title("COLO320 scATAC br15 minsizeCNV=0 compared to WGS")
    for border in borders_colorep1:
        plt.axvline(border, color='gray')
    plt.xlim((0, len(atac_array)))
    plt.legend()
    plt.show()

    plt.show()
    return(scipy.stats.pearsonr(standarize(atac_array),standarize(smoothed_wgs_array)))
def createLinePlot(loss_wgs, base_wgs, gain_wgs, loss_atac, base_atac, gain_atac):
    new_base_wgs = [x * 2 for x in base_wgs]
    new_base_atac = [x * 2 for x in base_atac]
    new_gain_wgs = [x * 3 for x in gain_wgs]
    new_gain_atac = [x * 3 for x in gain_atac]
    wgs_plot=[sum(x) for x in zip(new_gain_wgs,new_base_wgs,loss_wgs)]
    atac_plot = [sum(x) for x in zip(new_gain_atac, new_base_atac, loss_atac)]
    atac_array=np.array(atac_plot)
    wgs_array=np.array(wgs_plot)
    print("Pearson Correlation : ",scipy.stats.pearsonr(standarize(atac_array), standarize(wgs_array)))
    print("Spearman Correlation : ", scipy.stats.spearmanr(standarize(atac_array), standarize(wgs_array)))
    print("Kendall Correlation : ", scipy.stats.kendalltau(standarize(atac_array), standarize(wgs_array)))

    difference_array = np.subtract(atac_array, wgs_array)
    squared_array = np.square(difference_array)
    mse = squared_array.mean()
    print("Mean Square Error: ",squared_array.mean())

    print("R square: ",r2_score(atac_array, wgs_array))
    print("Mean Absolute Error: ",mean_absolute_error(atac_array,wgs_array))
    x = list(range(len(wgs_plot)))
    borders_colorep1 = [0, 2109, 4392, 6241, 8068, 9661, 11085, 12484, 13785, 14854, 16070, 17287, 18527, 19236, 20058,
                        20758, 21504, 22212, 22836, 23312, 23867, 24172, 24482]
    borders=[0,2231,4594,6543,8408,10163,11846,13388,14804,15924,17226,18531,19832,20785,21660,22447,23214,23988,24728,25238,25884,26221,26560]
    #borders_1e6=[0,218,450,643,823,996,1159,1311,1449,1558,1684,1809,1935,2027,2113,2189,2264,2339,2410,2463,2520,2552,2583]
    #borders_101=[0,2284,4672,6645,8528,10362,12021,13599,15037,16229,17545,18877,20197,21167,22063,22897,23705,24521,25304,25883,26495,26864,27233]
    #borders_10000_85=[0,2230,4592,6541,8406,10161,11844,13386,14802,15921,17221,18525,19826,20779,21654,22441,23207,23980,24719,25274,25875,26212,26551]
    #plt.plot(x, wgs_plot, color='blue', label="GS")
    #plt.plot(x, atac_plot, color='orange', label="ATAC")
    #plt.plot(x, standarize(wgs_array), color='blue', label="GS")
    #plt.plot(x,standarize(atac_array),color='orange', label="ATAC")
    #both = np.concatenate([standarize(wgs_array)[:, None], standarize(atac_array)[:, None]], axis=1)
    #x = list(range(len(w../../Colo320HSP/aneufinder/COLO320HSR_WGS_perBin_noChr.bedgs_array)))
    plt.plot(x, (normalize(atac_array)), color='#df979e', label="ATAC")
    plt.plot(x, (normalize(wgs_array)), color='#98d1d1', label="WGS")
    #plt.plot(both)
    for border in borders_colorep1:
        plt.axvline(border, color='gray')
    plt.title("COLO320 scATAC br15 minsizeCNV=0 compared to WGS (CNVkit)")
    plt.xlim((0, len(atac_plot)))
    plt.legend()
    plt.show()


if __name__ =="__main__":
    p = {"filter_edges": "nearest",
         "atac_input": "/home/katia/Helmholz/epiAneufinder/revisions/GSM4861367_COLO320HSR_rep1_atac/epiAneufinder_results/colo320HSP_rep1_results_table_nochr.tsv",
         "wgs_reads": "/home/katia/Helmholz/epiAneufinder/Colo320HSP/aneufinder/COLO320HSR_WGS_perBin_noChr.bed",
         "smooth_sigma": .1,
         "title": "WGS vs scATAC"}
    p, _ = fargv.fargv(p)
    harmonics = {}
    #fin=open("/home/katia/Helmholz/epiAneufinder/Colo320HSP/aneufinder/COLO320_HSR_WGS.CNV_nochr.bed")
    fin = open(p.wgs_reads)
    snu_full=pd.read_csv(p.atac_input, sep=" ")
    #l = ["seq", "start", "end"]
    #cellCluster = open("/home/katia/Helmholz/epiAneufinder/revisions/GSM4861371_COLO320HSR_rep3_atac/epiAneufinder_results/cluster2_cells.txt", "r")
    #lines = cellCluster.readlines()
    #for line in lines:
    #    l.append(line.strip())
    #snu_select=snu_full.loc[:,l]
    snu_dict=createDictionaryFromTable(snu_full)
    bed_dict=createDictionaryFromBed(fin)
    common_wgs_dict, common_atac_dict = calculatePopulationSomiesWGS(snu_dict,bed_dict)
    #for sigma in range(1, 10, 1):
    #    p.smooth_sigma = sigma
    for _ in (0,):
        correlation=createLinePlotAneufinder(common_wgs_dict, common_atac_dict,gaussian_sigma=p.smooth_sigma, filter_edges=p.filter_edges)
        harmonics[p.smooth_sigma] = correlation
        print(harmonics)

    #loss_wgs, base_wgs, gain_wgs, loss_atac, base_atac, gain_atac = calculatePopulationSomies(snu_dict, bed_dict)
    #createLinePlotAneufinder(loss_wgs, base_wgs, gain_wgs, loss_atac, base_atac, gain_atac)
    #print(gain_atac)
    #print(loss_wgs,loss_atac)
