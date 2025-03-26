import pandas as pd
import numpy as np
import json
import sys

inSeqDir = sys.argv[1] #"intermediates/rawAlign/"
inPTMDir = sys.argv[2] #"intermediates/ptm/"
ptms = sys.argv[3].split("@@") #["Phospho","Acetyl"]
kids = sys.argv[4].split("@@") #["K00016", "K01689"]
names = sys.argv[5].split("@@") #["L-Lactase dehydrogenase","Enolayse"]
outFile = sys.argv[6] #'raw_ptms.csv'

namesDict = dict(zip(kids,names))

kidFrames = dict()
for k in kids:
    inSeqFile = inSeqDir+"/"+k+".faa"
    inSeq = dict()
    with open(inSeqFile, "r") as inFile:
        line = inFile.readline().strip()
        while line:
            if line[0] == ">":
                pid = line[1:]
                inSeq[pid] = []
            else:
                inSeq[pid].append(line)
            line = inFile.readline().strip()
    inSeq = {i:list("".join(j)) for i,j in inSeq.items()}
    conSeq = "".join(pd.DataFrame(inSeq).T.mode().iloc[0])
    ptm_mask_list = []
    for p in ptms:
        inPTMFile = inPTMDir+"/"+k+"_"+p+"_aligned.json"
        with open(inPTMFile,"r") as inFile:
            ptm_sites = json.load(inFile)

        for i in ptm_sites.keys():
            locs = ptm_sites[i]["mod_site"]
            sites = np.zeros(len(conSeq))
            sites[locs] = 1
            ptm_sites[i] = sites

        ptm_sites = pd.DataFrame(ptm_sites).T
        ptm_sites["org"] = ptm_sites.index.map(lambda x: x.split(", ", maxsplit=1)[1])
        ptm_mask = ptm_sites.groupby("org").sum() > 0
        ptm_mask_list.append(ptm_mask)
    mods = pd.DataFrame(np.full(ptm_mask_list[0].shape, ""), index=ptm_mask_list[0].index)
    for i in range(len(ptms)-1,-1,-1):
        mods.mask(ptm_mask_list[i],ptms[i], inplace=True)
    mods.rename(columns=lambda x: conSeq[x]+str(x+1), inplace=True)
    kidFrames[namesDict[k]] = mods
bigTable = pd.concat(kidFrames, axis=1, names=["Protein", "Site"])
bigTable.to_csv(outFile)
