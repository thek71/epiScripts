import sys

fin=open(sys.argv[1],'r')
line=fin.readline()
lines=fin.readlines()
for line in lines:
    #print(line)
    if line.strip()[0]=="t":
        print(line.strip())
    else:
        chrom=line.strip().split(",")[0]
        start = int(line.strip().split(",")[1])
        end = int(line.strip().split(",")[2])
        somy = float(line.strip().split(",")[-1])
        new_start=start
        while start<end:
            new_end=start+100000
            if start==0:
                print(chrom+"\t"+str(start)+"\t"+str(new_end)+"\t"+str(somy))
            else:
                print(chrom + "\t" + str(start) + "\t" + str(new_end-1) + "\t" + str(somy))
            start=start+100000

