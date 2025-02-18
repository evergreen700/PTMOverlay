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
for item,ptms in inPTM.items():
    seq = "".join(inSeq[item])
    seq_progress=-1
    ptm_iterator = 0
    lo_ptms = []
    ptms.sort()
    for i in range(len(seq)):
        if seq[i] != "-":
            seq_progress+=1
        if seq_progress == ptms[ptm_iterator]:
            lo_ptms.append(i)
            ptm_iterator+=1
            if ptm_iterator == len(ptms):
                break
    outPTM[item] = lo_ptms

with open(outPTMFile, "w") as outFile:
    json.dump(outPTM, outFile, indent=4)
