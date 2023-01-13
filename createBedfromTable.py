import sys

fin=open(sys.argv[1],"r")
line=fin.readline()
lines=fin.readlines()
print("Chr\tStart\tEnd")
for line in lines:
    print(line.strip().split(" ")[1]+"\t"+line.strip().split(" ")[2]+"\t"+line.strip().split(" ")[3])