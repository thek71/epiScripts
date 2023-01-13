import sys
import seaborn as sns
import pandas as pd
import numpy as np
import fargv
import scipy.stats
from collections import defaultdict
from matplotlib import pyplot as plt
#plt.style.use('seaborn-whitegrid')
#sns.set_theme()

def createDictionaryFromTable(table):
    snu_dict=table.set_index(['seqnames', 'start', 'end']).T.to_dict('list')
    return(snu_dict)

def countNplicitiesFromTable(table, chr_filter="All"):
    atac_dict=table.set_index(['seqnames', 'start', 'end']).T.to_dict('list')
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

def preprocessATAC(loss_atac, base_atac, gain_atac):
    new_base_atac = [x * 2 for x in base_atac]
    new_gain_atac = [x * 3 for x in gain_atac]
    atac_plot = np.array([sum(x) for x in zip(new_gain_atac, new_base_atac, loss_atac)])
    #x = list(range(len(wgs_plot)))
    return(atac_plot)

if __name__ =="__main__":
    "  "
    p = {"chr": "All",
         "filter_edges": "nearest",
        "atac_input":"/home/katia/Helmholz/epiAneufinder/results/SNU601/somiesdt_with_coord_info.tsv",
         "title":"WES vs scATAC",
         "label1":"label1"}
    p,_ = fargv.fargv(p)
    atac_loss, atac_base, atac_gain = countNplicitiesFromTable(pd.read_csv(p.atac_input, delimiter=","), chr_filter=p.chr)
    atac_array=preprocessATAC(atac_loss, atac_base, atac_gain)
    borders_6t = [0, 1556, 2971,4077,4903,5834,6857,7731,8466,9113,9906,10700,11507,11955,12462,12999,13485,14089,14469,14933,15329,15490,15744]
    borders_8t=[0,1329,2505,3445,4084,4837,5713,6429,7039,7577,8275,8934,9645,10014,10497,10962,11388,11946,12234,12668,13003,13146,13386]
    borders_8i=[0,2216,4572,6517,8376,10130,11799,13342,14747,15837,17110,18416,19713,20659,21533,22321,23094,23857,24595,25145,25730,26057,26391]
    borders_6i=[0,2224,4584,6530,8389,10144,11813,13359,14768,15864,17139,18445,19742,20688,21563,22357,23133,23899,24637,25188,25113,26102,26439]
    su008_chr_borders = [0, 2174, 4386, 6265, 7913, 9591, 11106, 12563, 13924, 14966, 16164, 17371, 18611, 19431, 20210,
                         20973, 21678, 22416, 23121, 23649, 24231, 24526, 24854]
    su006_pre_chr_borders = [0, 2178, 4353, 6297, 8110, 9847, 11491, 13011, 14655, 16175, 17577, 18584, 19840, 21135,
                             22418, 23360, 24229, 25008, 25776, 26529, 27261, 27803, 28386, 28699, 29022]
    borders_8i_post=[]
    borders_8t_post=[0,1491,2839,3939,4757,5723,6759,7614,8334,8990,9780,10552,11384,11834,12362,12894,13377,13981,14342,14804,15182,15349,15598]

    borders_6t_6 = [0,217,449,642,822,995,1158,1310,1448,1555,1681,1806,1932,2024,2110,2186,2261,2336,2408,2461,2517,2549,2580]
    borders_8t_6 = [0,217,449,642,822,995,1158,1310,1447,1554,1680,1805,1931,2023,2109,2185,2259,2334,2406,2459,2515,2547,2578]
    borders_8i_6 = [0,217,449,642,822,995,1158,1310,1448,1557,1683,1808,1934,2026,2112,2188,2263,2338,2410,2463,2519,2551,2582]
    borders_6i_6 = [0,217,449,642,822,995,1158,1310,1448,1558,1684,1809,1935,2027,2113,2189,2264,2339,2411,2464,2520,2552,2583]
    #su008_chr_borders_6 = [0]
    #su006_pre_chr_borders_6 = [0]
    borders_8i_post_6 = []
    borders_8t_post_6 = []

    plt.plot(atac_array, label=p.label1)
    for border in borders_6t_6:
        plt.axvline(border, color='gray')
    plt.xlim((0, len(atac_array)))
    plt.legend()
    plt.title(p.title)
    plt.show()
