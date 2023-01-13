import sys

fin=open(sys.argv[1],"r")

lines=fin.readlines()
countGain=0
countDisomy=0
countLoss=0
countGainAtac=0
countDisomyAtac=0
countLossAtac=0
countGainWGS=0
countDisomyWGS=0
countLossWGS=0
for line in lines:
    atacCNV=int(line.strip().split("\t")[3])
    wgsCNV=int(line.strip().split("\t")[8])
    bases=int(line.strip().split("\t")[-1])
    if atacCNV==2 and wgsCNV==2:
        countDisomy=countDisomy+1
    if atacCNV==3 and wgsCNV>=3:
        countGain=countGain+1
        #print(line)
    if atacCNV==1 and wgsCNV<=1:
        countLoss=countLoss+1
    if atacCNV==2 and wgsCNV!=2:
        if wgsCNV<=1:
            countLossWGS=countLossWGS+1
            #print(line.strip())
        else:
            countGainWGS=countGainWGS+1
            #print(line.strip())
    if atacCNV!=2 and wgsCNV==2:
        if atacCNV==1:
            countLossAtac=countLossAtac+1
            #print(line.strip())
        else:
            countGainAtac=countGainAtac+1
print("Common Gains:"+str(countGain),"("+str((countGain/len(lines))*100)+"% of total ATAC windows), ("+str((countGain/(countGainWGS+countGainAtac+countGain))*100)+"% of all gains)")
print("Common Losses:"+str(countLoss), "("+str((countLoss/len(lines))*100)+"% of total ATAC windows), ("+str((countLoss/(countLossWGS+countLossAtac+countLoss))*100)+"% of all losses)")
print("Common Disomy:"+str(countDisomy), "("+str((countDisomy/len(lines))*100)+"% of total ATAC windows)")
print("Unique Gains ATAC:"+str(countGainAtac), "("+str((countGainAtac/len(lines))*100)+"% of total ATAC windows), ("+str((countGainAtac/(countGainWGS+countGainAtac+countGain))*100)+"% of all gains)")
print("Unique Gains WGS:"+str(countGainWGS), "("+str((countGainWGS/len(lines))*100)+"% of total ATAC windows), ("+str((countGainWGS/(countGainWGS+countGainAtac+countGain))*100)+"% of all gains)")
print("Unique Losses ATAC:"+str(countLossAtac), "("+str((countLossAtac/len(lines))*100)+"% of total ATAC windows), ("+str((countLossAtac/(countLossWGS+countLossAtac+countLoss))*100)+"% of all losses)")
print("Unique Losses WGS:"+str(countLossWGS), "("+str((countLossWGS/len(lines))*100)+"% of total ATAC windows), ("+str((countLossWGS/(countLossWGS+countLossAtac+countLoss))*100)+"% of all losses)")
