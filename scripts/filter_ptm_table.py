import pandas as pd
import sys

inTable = sys.argv[1]
outTable = sys.argv[2]

cutoff = 15
nameReplace = lambda x: "Methyl" if x.endswith("Methyl") else x

bigTable = pd.read_csv(inTable, header=[0,1], index_col=0, dtype=str)

bigTable.fillna("",inplace=True)

bigTable = bigTable.map(nameReplace)
mask = (bigTable.map(lambda x: x not in ["", "-"])).astype(int).sum()>cutoff
bigTable = bigTable[bigTable.columns[mask]]

bigTable.to_csv(outTable)
