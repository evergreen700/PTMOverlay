from pyteomics import pepxml
import glob
import os
import multiprocessing
from collections import defaultdict
from Bio import SeqIO
from pathlib import Path
import re
import json
import sys

modification_dir = sys.argv[1] 
modification_files = glob.glob(os.path.join(modification_dir, "**", "*.pepXML"), recursive=True)

genebank_dir = sys.argv[2]
genebank_files = glob.glob(os.path.join(genebank_dir, "**", "*.faa"), recursive=True)

outdir = sys.argv[3]
os.makedirs(outdir, exist_ok=True)

if len(sys.argv) > 4:
    outsuffix = "_"+sys.argv[4]
else:
    outsuffix = ""

UMB_to_GCA = {
    'UMB0490': 'GCA_002847685.2', 'UMB0064': 'GCA_002847765.1', 'UMB1298': 'GCA_002847825.1', 
    'UMB0085': 'GCA_002861815.1', 'UMB0386': 'GCA_002861965.1', 'UMB0388': 'GCA_002871555.1',
    'UMB0411': 'GCA_002871635.1', 'UMB0089': 'GCA_002871815.1', 'UMB0049': 'GCA_002884585.1',
    'UMB0663': 'GCA_002884815.1', 'UMB0088': 'GCA_002884955.3', 'UMB0683': 'GCA_002940945.1',
    'UMB3669': 'GCA_003286645.3', 'UMB0574': 'GCA_003286795.2', 'UMB0607': 'GCA_007785995.1',
    'UMB0637': 'GCA_008726885.1', 'UMB0572': 'GCA_008727025.1', 'UMB0248': 'GCA_008727075.1',
    'UMB0639': 'GCA_008728115.1', 'UMB0613': 'GCA_030218505.1', 'UMB0725': 'GCA_030218525.1',
    'UMB0612': 'GCA_030218565.1', 'UMB0029': 'GCA_030218815.1', 'UMB0016': 'GCA_030218925.1',
    'UMB2445': 'GCA_030222485.2', 'UMB0026': 'GCA_030225245.1', 'UMB0025': 'GCA_030225265.1',
    'UMB0734': 'GCA_030226495.1', 'UMB0005': 'GCA_030226535.1', 'UMB0036': 'GCA_030230945.1',
    'UMB0434': 'GCA_019890915.1'
}

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
        return sequence.index(peptide) + loc - 1
    except ValueError:
        return None

def search_peptide(peptide, protein_id, UMB, location):
    """Searches for a peptide within a protein sequence stored in memory."""
    GCA = UMB_to_GCA.get(UMB, 'NA') + ".faa"
    sequence = genebank_dict.get(GCA, {}).get(protein_id)
    if sequence:
        return get_index(peptide, sequence, location)
    return None

UMB_pattern = re.compile(r"UMB\d{4}")

def process_file(file):
    """Extracts unique proteins and site indexes from a .pepXML file."""
    UMB = UMB_pattern.search(os.path.basename(file)).group()
    assembly = UMB_to_GCA.get(UMB, 'NA')
    ka = os.path.join(genebank_dir, assembly+".kegg.txt")
    orthologs = dict()
    with open(ka,"r") as inFile:
        for l in inFile:
            pair = l.split()
            if len(pair) == 2:
                orthologs[pair[0]]=pair[1]

    ortho_org = dict()

    try:
        for entry in pepxml.read(file):
            search_hit = entry.get('search_hit', [{}])[0]
            protein = search_hit.get('proteins', [{}])[0].get('protein', "")
            if not protein or protein not in orthologs or protein.startswith("rev"):
                continue
            UMB_match = UMB_pattern.search(entry.get('spectrum', ""))
            if not UMB_match:
                continue
            UMB_string = UMB_match.group()
            modifications = search_hit.get('modifications', [])
            if modifications:
                mod_index = int(modifications[-1]['position'])
                peptide = search_hit.get('peptide', "")
                mod_index = search_peptide(peptide, protein, UMB_string, mod_index)
                if mod_index is not None:
                    o = orthologs[protein]
                    if o not in ortho_org:
                        ortho_org[o] = defaultdict(set)
                    ortho_org[o][protein+", "+assembly].add(mod_index)
    except Exception as e:
        print(f"Error processing {file}: {e}")

    return ortho_org

protein_ids = set()
subsequence_indexes = dict()

if __name__ == "__main__":
    ctx = multiprocessing.get_context("spawn")

    with ctx.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = pool.imap_unordered(process_file, modification_files)
        for index in results:
            for key, value in index.items():
                if key not in subsequence_indexes:
                    subsequence_indexes[key] = defaultdict(set)
                for k, v in value.items():
                    subsequence_indexes[key][k].update(v)
    for kid, a in subsequence_indexes.items():
        indexes = {i:list(j) for i,j in a.items()}
        with open(os.path.join(outdir,kid+outsuffix+".json"), "w") as file:
            json.dump(indexes, file, indent=4)
