import os
import sys
from collections import deque
import json

inFile = sys.argv[1]
inPTM = sys.argv[2]
outFile = sys.argv[3]
sequences = []
with open(inFile,"r") as readIn:
    seq = deque([readIn.readline()[1:]])
    line = readIn.readline()
    while line:
        if line[0] == ">":
            sequences.append(seq)
            seq = deque([line[1:]])
        else:
            seq.append(line)
        line=readIn.readline()
sequences.append(seq)

with open(inPTM,"r") as readIn:
    ptms = json.load(readIn)

with open(outFile, "w") as writer:
    #print header
    keyOrder = []
    for i in range(len(sequences)):
        seq = sequences[i].popleft()
        writer.write(f"{i:>2} "+seq)
        keyOrder.append(seq.strip())
    print(keyOrder)
    writer.write("\n\n")

    while len(sequences[0]) != 0:
        iters = 0
        for i in range(len(sequences)):
            writer.write(f"{i:>2} ")
            seq = sequences[i].popleft().strip()
            idx = iters
            for a in seq:
                if a != "-" and idx in ptms[keyOrder[i]]:
                    writer.write(a.lower())
                else:
                    writer.write(a)
                idx+=1
            writer.write('\n')
        iters=idx
        if len(sequences[0]) != 0:
            writer.write("\n\n")
