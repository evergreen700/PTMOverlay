#extract GCA codes from tsv file, convert to GCF codes and save as IDs.txt
import pandas as pd
import sys

outfile = sys.argv[1]
infile = sys.argv[2]

data = pd.read_csv(infile, sep="\t")

def modify_string(s):
    return s.replace("A", "F", 1)

data["Assembly"] = data["Assembly"].apply(modify_string)

data["Assembly"].to_csv(outfile, index=False, header=False)