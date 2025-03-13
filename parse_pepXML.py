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
import yaml

modification_dir = sys.argv[1] 
modification_files = glob.glob(os.path.join(modification_dir, "**", "*.pepXML"), recursive=True)

genebank_dir = sys.argv[2]
genebank_files = glob.glob(os.path.join(genebank_dir, "**", "*.faa"), recursive=True)

ptm_info = sys.argv[3]

outdir = sys.argv[4]
os.makedirs(outdir, exist_ok=True)

with open(ptm_info, "r") as inFile:
    aa_mass = yaml.load(inFile, Loader=yaml.CLoader)

with open("config.yaml") as inFile:
    ptm_mass = yaml.load(inFile, Loader=yaml.CLoader)["ptm_types"]

if len(sys.argv) > 5:
    outsuffix = "_"+sys.argv[5]
    shift_size=ptm_mass[sys.argv[5]]
else:
    outsuffix = ""
    shift_size = 0

aa_mass = aa_mass["base"]


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
        i = sequence.index(peptide)
        #filter out locs that are 0/-1
        return i, i+len(peptide), [i + int(m["position"]) - 1 for m in loc if int(m["position"])>0 and aa_mass[sequence[int(m['position'])]]+shift_size == int(m['mass'])]
    except ValueError:
        return None

def search_peptide(peptide, protein_id, UMB, locations):
    """Searches for a peptide within a protein sequence stored in memory."""
    GCA = UMB_to_GCA.get(UMB, 'NA') + ".faa"
    sequence = genebank_dict.get(GCA, {}).get(protein_id)
    if sequence:
        return get_index(peptide, sequence, locations)
    return None

UMB_pattern = re.compile(r"UMB\d{4}")

def process_file(file):
    """Extracts unique proteins and site indexes from a .pepXML file."""
    UMB = UMB_pattern.search(os.path.basename(file)).group()
    assembly = UMB_to_GCA.get(UMB, 'NA')
    ka = os.path.join(genebank_dir, assembly+".kegg.txt")
    orthologs = dict()
    ortho_org = dict()
    with open(ka,"r") as inFile:
        for l in inFile:
            pair = l.split()
            if len(pair) == 2:
                orthologs[pair[0]]=pair[1]
                if pair[1] not in ortho_org:
                    ortho_org[pair[1]]={pair[0]+", "+assembly: {'read_start':dict(),'read_end':dict(),'mod_site':set()}}
                else:
                    ortho_org[pair[1]][pair[0]+", "+assembly] = {'read_start':dict(),'read_end':dict(),'mod_site':set()}
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
                peptide = search_hit.get('peptide', "")
                startIdx, endIdx, mod_indices = search_peptide(peptide, protein, UMB_string, modifications)
                o = orthologs[protein]
                ortho_org[o][protein+", "+assembly]['read_start'][startIdx]=ortho_org[o][protein+", "+assembly]['read_start'].get(startIdx,0)+1
                ortho_org[o][protein+", "+assembly]['read_end'][endIdx]=ortho_org[o][protein+", "+assembly]['read_end'].get(endIdx,0)+1
                ortho_org[o][protein+", "+assembly]['mod_site'].update(mod_indices)
    except Exception as e:
        print(f"Error processing {file}: {e}")

    return ortho_org

protein_ids = set()
subsequence_indexes = dict()

if __name__ == "__main__":
    ctx = multiprocessing.get_context("forkserver")

    with ctx.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = pool.imap_unordered(process_file, modification_files)
        for index in results:
            for key, value in index.items():
                if key not in subsequence_indexes:
                    subsequence_indexes[key] = dict()
                for k, v in value.items():
                    if k not in subsequence_indexes[key]:
                        subsequence_indexes[key][k] = {'read_start':dict(),'read_end':dict(),'mod_site':set()}

                    for ki,vi in v['read_start'].items():
                        if ki not in subsequence_indexes[key][k]['read_start']:
                            subsequence_indexes[key][k]['read_start'][ki] = 0
                        subsequence_indexes[key][k]['read_start'][ki]+=vi

                    for ki,vi in v['read_end'].items():
                        if ki not in subsequence_indexes[key][k]['read_end']:
                            subsequence_indexes[key][k]['read_end'][ki] = 0
                        subsequence_indexes[key][k]['read_end'][ki]+=vi

                    subsequence_indexes[key][k]['mod_site'].update(v['mod_site'])
    for kid, a in subsequence_indexes.items():
        for i in a.keys():
            a[i]["mod_site"] = list(a[i]["mod_site"])
        #indexes = {i:list(j) for i,j in a.items()}
        with open(os.path.join(outdir,kid+outsuffix+".json"), "w") as file:
            json.dump(a, file, indent=4)
