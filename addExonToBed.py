import sys
fin=open(sys.argv[1],"r")
count=0
lines=fin.readlines()
for line in lines:
    chr=line.strip().split("\t")[0]
    start = line.strip().split("\t")[1]
    end = line.strip().split("\t")[2]
    count=count+1
    print(chr+"\t"+start+"\t"+end+"\t"+"exon"+str(count))