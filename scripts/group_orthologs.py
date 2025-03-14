import sys
import glob
import os
import csv

PROTEOME_PATH = sys.argv[1]
ORTHOLOG = sys.argv[2]
SPECIES_INFO = sys.argv[3]
PRE_ALIGN_PATH = sys.argv[4]
covered_kos = set()
proteomes = glob.glob(os.path.join(PROTEOME_PATH,"*.faa"))

os.makedirs(PRE_ALIGN_PATH, exist_ok=True)

GCA_to_species = {}
with open(SPECIES_INFO, mode="r", newline="") as file:
    reader = csv.reader(file, delimiter="\t")
    headers = next(reader)
    for row in reader:
        key = row[3]
        value = (row[2])
        GCA_to_species[key] = value

for p in proteomes:
    ka = p[:-4]+".kegg.txt"
    assembly = os.path.basename(p)[:-4]
    orthologs = dict()
    with open(ka,"r") as inFile:
        for l in inFile:
            pair = l.split()
            if pair[-1] == ORTHOLOG:
                orthologs[pair[0]]=pair[1]
    with open(p,"r") as inFile:
        line = inFile.readline()
        pid = line.split(maxsplit=1)[0][1:]
        while pid not in orthologs:
            line = inFile.readline()
            while line and line[0] != ">":
                line = inFile.readline()
            if not line:
                break
            pid = line.split(maxsplit=1)[0][1:]
        if not line:
            continue
        ko = orthologs[pid]
        if ko in covered_kos:
            outFile = open(os.path.join(PRE_ALIGN_PATH,ko+".faa"),"a")
        else:
            covered_kos.add(ko)
            outFile = open(os.path.join(PRE_ALIGN_PATH,ko+".faa"),"w")
         
        line = ">"+pid+", "+assembly+", "+GCA_to_species[assembly]+"\n"

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
                line = ">"+pid+", "+assembly+", "+GCA_to_species[assembly]+"\n"
            outFile.write(line)
            line = inFile.readline()
        outFile.close()
