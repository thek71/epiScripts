import matplotlib.pyplot as plt
import numpy as np


labels = ['90%', '80%', '70%', '60%','50%','40%','30%','20%','10%']
full_labels = ['100','90%', '80%', '70%', '60%','50%','40%','30%','20%','10%']
full_labels_pos=[i for i,_ in enumerate(full_labels)]
labels_pos = [i for i, _ in enumerate(labels)]
print(labels_pos)
bins_cnv = [26566,26554,26538,26523,26507,26491,26453,26331,26147]
cells_cnv = [3123,3089,3043,2977,2883,2748,2489,1911,381]
percent_cells_cnv=[100*(i/3157) for i in cells_cnv]
percent_bins_cnv=[100*(i/26571) for i in bins_cnv]
plt.plot(labels_pos, percent_cells_cnv, color='red', label="%cells")
plt.plot(labels_pos, percent_bins_cnv, color='black',label="%bin")
plt.xlabel("Percentage of sampling")
plt.xticks(labels_pos,labels)
plt.legend()
plt.title("Percentage of cells and bins per sub-sample SNU601")
#plt.show()
plt.savefig('/home/katia/Helmholz/epiAneufinder/results/subsampling/Subsampling_percent_cell_bins_per_subsampling.pdf')
plt.close()
#plt.xlabel("Percentage of sampling")
#plt.ylabel("Percentage of cells")
#plt.title("Percentage of cells in the sub-sample SNU601")
#plt.xticks(labels_pos,labels)
#plt.savefig('/home/katia/Helmholz/epiAneufinder/results/subsampling/Subsampling_percentage_cells_in_calls.pdf')
#plt.close()

fragments=[70117,64738,58944,52679,45843,38420,30298,21323,11447]
cells_sample=[3547,3555,3541,3534,3527,3528,3502,3462,3326]
fragments_percent=[100*(i/75013) for i in fragments]
cells_sample_percents=[100*(i/3575) for i in cells_sample]
#print(fragments_percent,"\n",cells_sample_percents)
#print(fragments_percent)
#labels_pos = [i for i, _ in enumerate(fragments)]
plt.plot(labels_pos, cells_sample_percents, color='green', label="%cells")
plt.plot(labels_pos, fragments_percent, color='blue',label="%fragments per cell")
plt.xlabel("Percentage of sampling")
#plt.ylabel("Percentage of fragments")
plt.title("Percentage of average cells and fragments per cell in the sub-sample SNU601")
plt.xticks(labels_pos,labels)
plt.legend()
#plt.show()
plt.savefig('/home/katia/Helmholz/epiAneufinder/results/subsampling/Subsampling_average_fragmment_and_cells_in_subsample.pdf')
plt.close()

reads=[147358.63,132921.31,118548.79,104138.51,89439.96,74679.96,59725.61,45124.65,30430.35,15838.18]
plt.plot(full_labels_pos,reads,marker='o', color='b',label="read pairs")
plt.xlabel("Percentage of sampling")
plt.xticks(full_labels_pos,full_labels)
plt.title("Average number of read pairs in sample")
plt.legend()
plt.savefig("/home/katia/Helmholz/epiAneufinder/results/subsampling/Subsampling_average_read_pairs.pdf")
plt.close()

consistent=[7961813,7656636,7327474,7027427,6680176,6205541,5406589,3863784,607115]
inconsistent=[1604297,2131863,2540924,2850808,3126338,3362239,3475098,3126932,733891]
percent_inconsistent=[100*(i/(i+j)) for i,j in zip(inconsistent,consistent)]
plt.plot(labels_pos, percent_inconsistent, color='blue', label="%inconsistent gains/losses")
plt.xlabel("Percentage of sampling")
plt.xticks(labels_pos,labels)
plt.legend()
plt.savefig("/home/katia/Helmholz/epiAneufinder/results/subsampling/Inconsistent_gains_and_losses.pdf")
plt.show()
plt.close()
#print(percent_inconsistent)




