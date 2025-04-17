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

modification_dir = sys.argv[1] 
modification_files = glob.glob(os.path.join(modification_dir, "**", "*.pepXML"), recursive=True)

genebank_dir = sys.argv[2]
genebank_files = glob.glob(os.path.join(genebank_dir, "**", "*.faa"), recursive=True)

ptm_info = sys.argv[3]
species_info = sys.argv[4]
outdir = sys.argv[5]
os.makedirs(outdir, exist_ok=True)

system = platform.system()

with open(ptm_info, "r") as inFile:
    aa_mass = yaml.load(inFile, Loader=yaml.CLoader)

with open("config.yaml") as inFile:
    ptm_mass = yaml.load(inFile, Loader=yaml.CLoader)["ptm_types"]

if len(sys.argv) > 6:
    outsuffix = "_"+sys.argv[6]
    shift_size=ptm_mass[sys.argv[6]]
else:
    outsuffix = ""
    shift_size = 0

aa_mass = aa_mass["base"]
ptm_mass = {j:i for i,j in ptm_mass.items()}


UMB_to_GCA = {}
# read data from tsv config into dictionary of umb: tuple(species, assembly)
with open(species_info, mode="r", newline="") as file:
    reader = csv.reader(file, delimiter="\t")
    headers = next(reader)
    for row in reader:
        key = row[1]
        value = (row[2], row[3])
        UMB_to_GCA[key] = value

#load the full protein sequences into memory
def load_genebank_sequences():
    """Preloads all .faa sequences into a nested dictionary {GCA: {protein_id: sequence}} for fast lookups."""
    genebank_dict = defaultdict(dict)
    for file in genebank_files:
        GCA_name = Path(file).name
        for record in SeqIO.parse(file, "fasta"):
            genebank_dict[GCA_name][record.id] = str(record.seq)
    return genebank_dict

genebank_dict = load_genebank_sequences()

def get_index(peptide, sequence, loc):
    """Finds peptide start index and adds location offset."""
    try:
        i = sequence.index(peptide)
        #filter out locs that are 0/-1
        sites = dict()
        for m in loc:
            if int(m["position"])>0:
                shift_size = int(m['mass']) - aa_mass[sequence[int(m['position'])]]
                if shift_size in ptm_mass:
                    sites[i + int(m["position"]) - 1] = ptm_mass[shift_size]
        return sites
    except ValueError:
        return None

def search_peptide(peptide, protein_id, UMB, locations):
    """Searches for a peptide within a protein sequence stored in memory."""
    GCA = UMB_to_GCA.get(UMB, 'NA')[1] + ".faa"
    sequence = genebank_dict.get(GCA, dict()).get(protein_id, None)
    if sequence:
        return get_index(peptide, sequence, locations)
    return None

UMB_pattern = re.compile(r"UMB\d{4}")

def process_file(file):
    """Extracts unique proteins and site indexes from a .pepXML file."""
    UMB = UMB_pattern.search(os.path.basename(file)).group()
    species, assembly = UMB_to_GCA.get(UMB, 'NA')
    ortho_org = defaultdict(dict)
#    try:
    for entry in pepxml.read(file):
        search_hit = entry.get('search_hit', [{}])[0]
        protein = search_hit.get('proteins', [{}])[0].get('protein', "")
        if not protein or  protein.startswith("rev"):
            continue
        UMB_match = UMB_pattern.search(entry.get('spectrum', ""))
        if not UMB_match:
            continue
        UMB_string = UMB_match.group()
        modifications = search_hit.get('modifications', [])
        if modifications:
            peptide = search_hit.get('peptide', "")
            mod_indices = search_peptide(peptide, protein, UMB_string, modifications)
            if mod_indices:
                if protein not in ortho_org[assembly]:
                    ortho_org[assembly][protein] = dict()
                ortho_org[assembly][protein].update(mod_indices)
#    except Exception as e:
#        print(f"Error processing {file}: {e}")

    return ortho_org

protein_ids = set()
subsequence_indexes = dict()

if __name__ == "__main__":
    context = "spawn" if system == "Windows" else "forkserver"
    ctx = multiprocessing.get_context(context)

    with ctx.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = pool.imap_unordered(process_file, modification_files)
        for index in results:
            for key, value in index.items():
                if key not in subsequence_indexes:
                    subsequence_indexes[key] = dict()
                for k, v in value.items():
                    if k not in subsequence_indexes[key]:
                        subsequence_indexes[key][k] = dict()
                    subsequence_indexes[key][k].update(v)

    genebank_dict = {k.removesuffix(".faa"):v for k,v in genebank_dict.items()}

    with open(os.path.join(outdir,"seqs.json"), "w") as file:
        json.dump(genebank_dict, file, indent=4)
    with open(os.path.join(outdir,"ptms.json"), "w") as file:
        json.dump(subsequence_indexes, file, indent=4)
    for kid, a in subsequence_indexes.items():
        with open(os.path.join(outdir,kid+".json"), "w") as file:
            json.dump(a, file, indent=4)
