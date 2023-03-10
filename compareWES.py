import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import fargv
import scipy
import scipy.ndimage

#Helping script for the comparison of WES data with the output of epineufinder.
#Calculated the pseudobulk dataset from the epiAneufinder results and compares them to the normalized WES dataset.
#Inputs are the csv output file from epiAneufinder, the WES count file of the reads for tumor and control (one file) for the genomic bins (calculated with featureCounts)
# and the count file with the number of exon bases per bin (calculated with bedtools).
#As output is the plot of the two signal modalities and the corresponding correlation coefficient.
#The user can decide which type of edge filtering to use and if all chromosomes or only one should be calculated.

standarize = lambda x: ((x - x.mean()) / (x.std() + .0000000000000001)) #The standarization formula for bringing the two datasets in the same scale

def countNplicitiesFromTable(table, chr_filter="All"):
    #Function that takes as input a csv file/output from epiAneufinder and calculates the pseudobulk dataset
    atac_dict=table.set_index(['seq', 'start', 'end']).T.to_dict('list')
    loss_atac=[]
    base_atac=[]
    gain_atac=[]
    for k in atac_dict:
        (chr, start, end) = k
        if chr == chr_filter or chr_filter=="All":
            loss_atac.append(atac_dict[k].count(0) / len(atac_dict[k]))
            base_atac.append(atac_dict[k].count(1) / len(atac_dict[k]))
            gain_atac.append(atac_dict[k].count(2) / len(atac_dict[k]))
    return(loss_atac,base_atac,gain_atac)

# def calculateWESTumorControl(wes_reads):
#     tumor_reads=[]
#     control_reads=[]
#     for l in wes_reads:
#         control_reads.append(int(l.strip().split("\t")[-1]))
#         tumor_reads.append(int(l.strip().split("\t")[-2]))
#     control_array = np.array(control_reads)
#     tumor_array = np.array(tumor_reads)
#     both = np.concatenate([standarize(control_array)[:, None], standarize(tumor_array)[:, None]], axis=1)
#     #plt.plot(np.convolve(both[:, 0], np.ones(100) / 100, mode='valid'), label="control")
#     #plt.plot(np.convolve(both[:, 1], np.ones(100) / 100, mode='valid'), label="tumor")
#     #plt.legend()
#     #plt.show()
#     return(both)


def preprocessATAC(loss_wgs, base_wgs, gain_wgs):
    #Scoring function of the CNVs. After calculating per bin the average CNVs per bin, gains/losses/disomies are scored with 3/1/2 accordingly.
    new_base_wgs = [x * 2 for x in base_wgs]
    new_gain_wgs = [x * 3 for x in gain_wgs]
    wgs_plot = np.array([sum(x) for x in zip(new_gain_wgs, new_base_wgs, loss_wgs)])
    #x = list(range(len(wgs_plot)))
    return(wgs_plot)

def preprocessWES(exon_coverage, wes_reads, gaussian_sigma=0, chr_filter="All", filter_edges="constant"):
    #Function that calculates the WES signal. Inputis the count file with the reads per bin and the exon bases file for the respective bins.
    fraction=[]
    coverage=[]
    for tsv_line in wes_reads:
        line = tsv_line.strip().split("\t")
        if line[1] == chr_filter or chr_filter == "All":
            if int(line[7]) != 0:
                fraction.append(int(line[6]) / int(line[7]))  # Normalization of the WES tumor signal by the control.
            else:
                fraction.append(1)
    for tsv_line in exon_coverage:
        line = tsv_line.strip().split("\t")
        if line[0] == chr_filter or chr_filter == "All":
            coverage.append(float(line[-1])) #Normalization of the WES tumor signal by the exon bases.
    normalised_wes = np.array([x * y for x, y in zip(fraction, coverage)])
    # normalised_wes = np.array([math.log((x * y),2) for x, y in zip(fraction, coverage)])
    # normalised_wes = np.array([math.log(x,2) for x in fraction])
    # normalised_wes = np.array([x for x in fraction])
    smoothed_wes = scipy.ndimage.gaussian_filter1d(normalised_wes, gaussian_sigma, axis=- 1, order=0, output=None, mode=filter_edges) #Gaussian filter on the normalized WES signal
    return(smoothed_wes)

if __name__ =="__main__":
    "  "
    p = {"chr": "All",
         "filter_edges": "nearest",
        "atac_input":"../SU008_Tumor_Post/SU008_Tumor_Post_finalTable.csv",
       "wes_reads": "../../WES/SU008/SU008_Tumor_Post_CNV.txt",
       "exon_coverage":"../../WES/SU008/SU008_Tumor_Post_CNVpositions.exonCoverage.txt",
       "smooth_sigma":.1,
         "title":"WES vs scATAC"}
    p,_ = fargv.fargv(p)
    harmonics = {}
    #By uncommenting the next three lines, commenting the for-loop and changing the range different sigmas can be calculated iteratively.
    #for sigma in range(4,40 , 2):
    #    p.smooth_sigma = sigma
        #atac_dict=createDictionaryFromTable(atac_input)
    for _ in (0,):
        atac_loss, atac_base, atac_gain = countNplicitiesFromTable(pd.read_csv(p.atac_input, sep=" "), chr_filter=p.chr)
        atac_array=preprocessATAC(atac_loss, atac_base, atac_gain)

        wes_reads=list(open(p.wes_reads,'r').readlines())[2:]



        lines_wes_coverage=open(p.exon_coverage,'r').readlines()
        bulk_array=preprocessWES(lines_wes_coverage, wes_reads, gaussian_sigma=p.smooth_sigma, chr_filter=p.chr, filter_edges=p.filter_edges)

        #bulk_array=np.array(normalised_wes)
        both = np.concatenate([standarize(bulk_array)[:, None], standarize(atac_array)[:, None]], axis=1)
        plt.plot(both)
        labels=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
        #Borders of chomosomes by number of bins for the protting function.
        su006_pre_chr_borders=[0,2148,4310,6247,8064,9794,11426,12934,14338,15360,16629,17916,19198,20141,21010,21779,22545,23310,24042,24588,25169,25501,25828]
        su008_pre=[0,2164,4420,6339,8076,9723,11335,12814,14209,15282,16526,17773,19040,19888,20705,21469,22188,22930,23654,24194,24777,25091,25410]
        su006_msCNV0=[0, 2118, 4261, 6141, 7940, 9623, 11105, 12589, 13967, 14979,16245, 17514, 18762, 19472, 20318, 21063, 21824, 22560, 23239, 23774, 24349, 24677,24990]
        su008_msCNV0=[0, 2135, 4372, 6235, 7956, 9558, 11025, 12484, 13855, 14918, 16160, 17392, 18627, 19267, 20062, 20802, 21518, 22234, 22905, 23434, 24011, 24325, 24631]
        for border in su008_msCNV0:
           plt.axvline(border, color='gray')
        plt.title(f"{p.title} for Chromosome {p.chr} \nCorrelation coefficient: {np.corrcoef(both.T)[1,0]:.5}")
        print(p.chr,{np.corrcoef(both.T)[1,0]:.5})
        difference_array = np.subtract(atac_array, bulk_array)
        squared_array = np.square(difference_array)
        mse = squared_array.mean()
        print("Mean Square Error: ", squared_array.mean())
        plt.legend([f"WES sigma: {p.smooth_sigma}", "scATAC"])
        plt.xlabel("epiAneufinder bins")
        plt.ylabel("Standardized values")
        plt.xlim((0,len(bulk_array)))

        plt.show()
        harmonics[p.smooth_sigma]=np.corrcoef(both.T)[1,0]
        print(harmonics)
        #python3 compareWES.py -atac_input=../../revisions/SU006pre_br15_msCNV0/epiAneufinder_results/SU006pre_br15_msCNV0_results_table.tsv -wes_reads=../../WES/SU006_pre_msCNV0/SU006pre_br15_msCNV0_CNV.txt -exon_coverage=../../WES/SU006_pre_msCNV0/SU006pre_br15_msCNV0_CNVpositions_exonCovergare.txt -smooth_sigma=12