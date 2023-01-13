import sys
import csv
import pandas as pd

tsv_file=(sys.argv[1])
filename=tsv_file.split(".")[0]
csv_table=pd.read_table(tsv_file,sep='\t')

csv_table.to_csv((filename+'.csv'),index=False)