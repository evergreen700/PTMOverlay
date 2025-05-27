from pyteomics import pepxml
import platform
import glob
import os
import multiprocessing
from collections import defaultdict
from Bio import SeqIO
from pathlib import Path
import re
import json
import sys
import yaml
import csv

#argument processing
modification_dir = sys.argv[1] 
in_fasta = sys.argv[2]
ptm_info = sys.argv[3]
species = sys.argv[4]
outdir = sys.argv[5]
strain = sys.argv[6]
#selection = sys.argv[8].split("@@")

assembly, _ = os.path.splitext(os.path.basename(in_fasta))
modification_files = glob.glob(os.path.join(modification_dir, "*"+strain+"*.pepXML"), recursive=True)
os.makedirs(outdir, exist_ok=True)

with open(ptm_info, "r") as inFile:
    aa_mass = yaml.load(inFile, Loader=yaml.CLoader)

with open("config.yaml") as inFile:
    ptm_mass = yaml.load(inFile, Loader=yaml.CLoader)["ptm_types"]

shift_sizes = {j:i for i,j in ptm_mass.items()}
aa_mass = aa_mass["base"]

#load the full protein sequences into memory
def load_genebank_sequences():
    """Preloads all .faa sequences into a nested dictionary {GCA: {protein_id: sequence}} for fast lookups."""
    genebank_dict = dict()
    for record in SeqIO.parse(in_fasta, "fasta"):
        genebank_dict[record.id] = str(record.seq)
    return genebank_dict

genebank_dict = load_genebank_sequences()

def get_index(peptide, sequence, loc):
    """Finds peptide start index and adds location offset."""
    try:
        i = sequence.index(peptide)
        #filter out locs that are 0/-1
        mods = []
        for l in loc:
            pos = int(l["position"]) - 1
            shift_size = int(l['mass']) - aa_mass[sequence[int(l['position'])]]
            mod_type = shift_sizes.get(mod_type,"")
            if pos >= 0 and mod_type:
                mods.append((mod_type,pos))

        return i, i+len(peptide), mods #[i + int(m["position"]) - 1 for m in loc if int(m["position"])>0 and aa_mass[sequence[int(m['position'])]]+shift_size == int(m['mass'])]
    except ValueError:
        return None

def search_peptide(peptide, protein_id, locations):
    """Searches for a peptide within a protein sequence stored in memory."""
    sequence = genebank_dict.get(protein_id)
    if sequence:
        return get_index(peptide, sequence, locations)
    return None

def process_file(files):
    """Extracts unique proteins and site indexes from a .pepXML file."""
    ka = in_fasta.removesuffix(".faa")+".kegg.txt"
    orthologs = dict()
    ortho_org = dict()
    ## test vals
    entry_suffix = ", ".join([strain, species, assembly])
    with open(ka,"r") as inFile:
        for l in inFile:
            pair = l.split()
            if len(pair) == 2: # and pair[1] in selection:
                orthologs[pair[0]]=pair[1]
                if pair[1] not in ortho_org:
                    ortho_org[pair[1]]={pair[0]+", "+entry_suffix: {'read_start':dict(),'read_end':dict(),'mod_sites':{p:set() for p in ptm_mass.keys()}}}
                else:
                    ortho_org[pair[1]][pair[0]+", "+entry_suffix] = {'read_start':dict(),'read_end':dict(),'mod_sites':{p:set() for p in ptm_mass.keys()}}}
    try:
        for file in files:
            print(file)
            for entry in pepxml.read(file):
                search_hit = entry.get('search_hit', [{}])[0]
                protein = search_hit.get('proteins', [{}])[0].get('protein', "")
                if not protein or protein not in orthologs or protein.startswith("rev"):
                    continue
                modifications = search_hit.get('modifications', [])
                if modifications:
                    peptide = search_hit.get('peptide', "")
                    startIdx, endIdx, mod_indices = search_peptide(peptide, protein, modifications)
                    o = orthologs[protein]
                    key = protein+", "+entry_suffix
                    ortho_org[o][key]['read_start'][startIdx]=ortho_org[o][key]['read_start'].get(startIdx,0)+1
                    ortho_org[o][key]['read_end'][endIdx]=ortho_org[o][key]['read_end'].get(endIdx,0)+1
                    for p, l in mod_indices:
                        ortho_org[o][key]['mod_site'][p].add(l)
    except Exception as e:
        print(f"Error processing {file}: {e}")

    return ortho_org

if __name__ == "__main__":
    results = process_file(modification_files)
    for kid in results.keys():
        a = results[kid]
        for i in a.keys():
            a[i]["mod_site"] = {p:list(sites) for p,sites in a[i]["mod_site"].items()}
    with open(os.path.join(outdir,strain+outsuffix+".json"), "w") as file:
        json.dump(results, file, indent=4)
