import os
import sys
from collections import deque

inFile = sys.argv[1]
outFile = sys.argv[2]
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

with open(outFile, "w") as writer:
    while len(sequences[0]) != 0:
        for i in range(len(sequences)):
            writer.write(f"{i:>2} "+sequences[i].popleft())
        if len(sequences[0]) != 0:
            writer.write("\n\n")
