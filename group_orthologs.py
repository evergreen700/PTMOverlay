import sys
import glob
import os

PROTEOME_PATH = sys.argv[1]
PRE_ALIGN_PATH = sys.argv[2]
covered_kos = set()
proteomes = glob.glob(os.path.join(PROTEOME_PATH,"*.faa"))
if not os.path.exists(PRE_ALIGN_PATH): 
    os.makedirs(PRE_ALIGN_PATH) 

for p in proteomes:
    ka = p[:-4]+".kegg.txt"
    assembly = os.path.basename(p)[:-4]
    orthologs = dict()
    with open(ka,"r") as inFile:
        for l in inFile:
            pair = l.split()
            if len(pair) == 2:
                orthologs[pair[0]]=pair[1]
    with open(p,"r") as inFile:
        line = inFile.readline()
        pid = line.split(maxsplit=1)[0][1:]
        while pid not in orthologs:
            line = inFile.readline()
            while line and line[0] != ">":
                line = inFile.readline()
            pid = line.split(maxsplit=1)[0][1:]
        ko = orthologs[pid]
        if ko in covered_kos:
            outFile = open(os.path.join(PRE_ALIGN_PATH,ko+".faa"),"a")
        else:
            covered_kos.add(ko)
            outFile = open(os.path.join(PRE_ALIGN_PATH,ko+".faa"),"w")
        line = ">"+pid+", "+assembly
        outFile.write(line)
        line = inFile.readline()
        while line:
            if line[0] == ">":
                outFile.close()
                pid = line.split(maxsplit=1)[0][1:]
                if pid not in orthologs:
                    line = inFile.readline()
                    while line and line[0] != ">":
                        line = inFile.readline()
                    continue
                ko = orthologs[pid]
                if ko in covered_kos:
                    outFile = open(os.path.join(PRE_ALIGN_PATH,ko+".faa"),"a")
                else:
                    covered_kos.add(ko)
                    outFile = open(os.path.join(PRE_ALIGN_PATH,ko+".faa"),"w")
                line = ">"+pid+", "+assembly+"\n"
            outFile.write(line)
            line = inFile.readline()
        outFile.close()
