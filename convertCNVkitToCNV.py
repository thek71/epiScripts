import sys

fin=open(sys.argv[1],"r") #input cnr file
line=fin.readline()
lines=fin.readlines()

for line in lines:
    chr=line.strip().split("\t")[0]
    start = line.strip().split("\t")[1]
    end = line.strip().split("\t")[2]
    logValue = line.strip().split("\t")[5]
    if float(logValue)<=float(-0.4):
        cnvstate="1"
        strand="-"
    if float(logValue)>-0.4 and  float(logValue)<=0.3 :
        cnvstate="2"
        strand="."
    if float(logValue)>0.3:
        cnvstate="3"
        strand="+"
    print (chr+"\t"+start+"\t"+end+"\t"+cnvstate+"\t0\t"+strand)
