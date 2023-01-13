import sys

fin=open(sys.argv[1])
lines=fin.readlines()
count=0

for line in lines:
    count=count+1
    chr=line.strip().split("\t")[0]
    start=line.strip().split("\t")[1]
    end=line.strip().split("\t")[2]
    strand=line.strip().split("\t")[5]
    score=line.strip().split("\t")[4]
    print(chr+"\t"+start+"\t"+end+"\t"+str(count)+"\t"+score+"\t"+strand)
