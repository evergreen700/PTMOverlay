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
print(inSeq)

with open(inPTMFile, "r") as inFile:
    inPTM = json.load(inFile)

print(inPTM)
