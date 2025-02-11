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

ptm_file = sys.argv[3]

UMB_to_GCA = {
    'UMB0490': 'GCA_002847685.2', 'UMB0064': 'GCA_002847765.1', 'UMB1298': 'GCA_002847825.1', 
    'UMB0085': 'GCA_002861815.1', 'UMB0386': 'GCA_002861965.1', 'UMB0388': 'GCA_002871555.1',
    'UMB0411': 'GCA_002871635.1', 'UMB0089': 'GCA_002871815.1', 'UMB0049': 'GCA_002884585.1',
    'UMB0663': 'GCA_002884815.1', 'UMB0088': 'GCA_002884955.2', 'UMB0683': 'GCA_002940945.1',
    'UMB3669': 'GCA_003286645.3', 'UMB0574': 'GCA_003286795.1', 'UMB0607': 'GCA_007785995.1',
    'UMB0637': 'GCA_008726885.1', 'UMB0572': 'GCA_008727025.1', 'UMB0248': 'GCA_008727075.1',
    'UMB0639': 'GCA_008728115.1', 'UMB0613': 'GCA_030218505.1', 'UMB0725': 'GCA_030218525.1',
    'UMB0612': 'GCA_030218565.1', 'UMB0029': 'GCA_030218815.1', 'UMB0016': 'GCA_030218925.1',
    'UMB2445': 'GCA_030222485.2', 'UMB0026': 'GCA_030225245.1', 'UMB0025': 'GCA_030225265.1',
    'UMB0734': 'GCA_030226495.1', 'UMB0005': 'GCA_030226535.1', 'UMB0036': 'GCA_030230945.1',
    'UMB0434': 'NA'
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
    proteins = set()
    mod_subsequence_indexes = defaultdict(set)

    try:
        for entry in pepxml.read(file):
            search_hit = entry.get('search_hit', [{}])[0]
            protein = search_hit.get('proteins', [{}])[0].get('protein', "")
            if not protein or protein.startswith("rev"):
                continue
            proteins.add(protein)
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
                    mod_subsequence_indexes[protein].add(mod_index)
    except Exception as e:
        print(f"Error processing {file}: {e}")

    return proteins, mod_subsequence_indexes

protein_ids = set()
subsequence_indexes = defaultdict(set)

if __name__ == "__main__":
    ctx = multiprocessing.get_context("forkserver")
    with ctx.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = pool.imap_unordered(process_file, modification_files)
        for result, index in results:
            protein_ids.update(result)
            for key, value in index.items():
                subsequence_indexes[key].update(value)
    subsequence_indexes = {i:list(j) for i,j in subsequence_indexes.items()}
#    with open("protein_ids.txt", "w") as file:
#        for item in protein_ids:d
#            file.write(f"{item}\nd")
    with open(ptm_file, "w") as file:
        json.dump(subsequence_indexes, file, indent=4)
