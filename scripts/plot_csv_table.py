import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys

inCSV = sys.argv[1]
outSVG = sys.argv[2]
dot_size=300

bigTable = pd.read_csv(inCSV, header=[0,1], index_col=0, dtype=str)
bigTable.fillna("",inplace=True)
npVals = bigTable.to_numpy()
npVals = np.flip(npVals, 0)

ptms = np.unique(npVals)
fig,ax = plt.subplots(1,1, figsize=(npVals.shape*np.array([.5, 1]))[::-1])
ax.set_yticks(np.arange(0,npVals.shape[0]), labels=bigTable.index[::-1])
ax.set_xticks(np.arange(0,npVals.shape[1]), labels=bigTable.columns.droplevel(level=0))
ax.set_xmargin(1/npVals.shape[1])
ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
sec = ax.secondary_xaxis(location="top")
proteins, idxs = np.unique(bigTable.columns.droplevel(level=1),return_index=True)
sec.set_xticks(idxs, labels=[s+"\n\n" for s in proteins])
ax.set_axisbelow(True)
ax.grid(True, axis='y')
for p in ptms:
    if len(p)>1:
        y,x = np.where(npVals == p)
        ax.scatter(x,y, label=p, marker="o", s=dot_size)
y,x = np.where(npVals == "")
ax.scatter(x,y, label="No PTM", color='white', edgecolors=["black"], marker="o", s=dot_size)

fig.tight_layout()
fig.legend(loc="upper left")

fig.savefig(outSVG, format="svg")
