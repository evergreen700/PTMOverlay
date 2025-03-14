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
    ptms = entries["mod_site"]
    seq = "".join(inSeq[item])
    seq_progress=-1
    ptm_iterator = 0
    lo_starts = dict()
    lo_ends = dict()
    lo_ptms = []
    ptms.sort()
    for i in range(len(seq)):
        if seq[i] != "-":
            seq_progress+=1
            if ptm_iterator < len(ptms) and seq_progress == ptms[ptm_iterator]:
                lo_ptms.append(i)
                ptm_iterator+=1
            if seq_progress in starts:
                lo_starts[i] = starts[seq_progress]
            if seq_progress in ends:
                lo_ends[i] = ends[seq_progress]
            if len(lo_ends) == len(ends):
                break
    outPTM[item] = {"read_start":lo_starts,"read_end":lo_ends,"mod_site":lo_ptms}

with open(outPTMFile, "w") as outFile:
    json.dump(outPTM, outFile, indent=4)
