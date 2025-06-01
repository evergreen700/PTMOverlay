import json
import sys

inPTMFile = sys.argv[1]
inSeqFile = sys.argv[2]
outPTMFile = sys.argv[3]

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

with open(inPTMFile, "r") as inFile:
    inPTM = json.load(inFile)

outPTM = dict()
for item,entries in inPTM.items():
    starts = {int(i):j for i,j in entries["read_start"].items()}
    ends = {int(i):j for i,j in entries["read_end"].items()}
    ptms = entries["mod_sites"]
    seq = "".join(inSeq[item])
    seq_progress=-1
    ptm_iterators = {k:0 for k in ptms.keys()} 
    lo_starts = dict()
    lo_ends = dict()
    lo_ptms = {k:[] for k in ptms.keys()}
    for i in ptms.keys():
        ptms[i].sort()
    for i in range(len(seq)):
        if seq[i] != "-":
            seq_progress+=1
            for p in ptms.keys():
                if ptm_iterators[p] < len(ptms[p]) and seq_progress == ptms[p][ptm_iterators[p]]:
                    lo_ptms[p].append(i)
                    ptm_iterators[p]+=1
            if seq_progress in starts:
                lo_starts[i] = starts[seq_progress]
            if seq_progress in ends:
                lo_ends[i] = ends[seq_progress]
            if len(lo_ends) == len(ends):
                break
    outPTM[item] = {"read_start":lo_starts,"read_end":lo_ends,"mod_sites":lo_ptms}

with open(outPTMFile, "w") as outFile:
    json.dump(outPTM, outFile, indent=4)
