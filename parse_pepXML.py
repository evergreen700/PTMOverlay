import glob
import os
from collections import defaultdict
from pathlib import Path
import re
import sys
import csv

modification_dir = sys.argv[1]
modification_files = glob.glob(os.path.join(modification_dir, "**", "*.pepXML"), recursive=True)

genebank_dir = sys.argv[2]
genebank_files = glob.glob(os.path.join(genebank_dir, "**", "*.faa"), recursive=True)

outdir = sys.argv[3]
os.makedirs(outdir, exist_ok=True)

if len(sys.argv) > 4:
    outsuffix = "_" + sys.argv[4]
else:
    outsuffix = ""

UMB_to_GCA = {}
# read data from tsv config into dictionary of umb: tuple(species, assembly)
with open("index_umb_taxa_gca.tsv", mode="r", newline="") as file:
    reader = csv.reader(file, delimiter="\t")
    headers = next(reader)
    for row in reader:
        key = row[1]
        value = (row[2], row[3])
        UMB_to_GCA[key] = value


# load the full protein sequences into memory
def load_genebank_sequences():
    """Preloads all .faa sequences into a nested dictionary {GCA: {protein_id: sequence}} for fast lookups."""
    from Bio import SeqIO

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
<<<<<<< Updated upstream
        return sequence.index(peptide) + loc - 1
    except ValueError:
        return None

def search_peptide(peptide, protein_id, UMB, location):
=======
        i = sequence.index(peptide)
        # filter out locs that are 0/-1
        return i, i + len(peptide), [i + int(m["position"]) - 1 for m in loc if int(m["position"]) > 0]
    except ValueError:
        return None


def search_peptide(peptide, protein_id, UMB, locations):
>>>>>>> Stashed changes
    """Searches for a peptide within a protein sequence stored in memory."""
    GCA = UMB_to_GCA.get(UMB, 'NA')[1] + ".faa"
    sequence = genebank_dict.get(GCA, {}).get(protein_id)
    if sequence:
        return get_index(peptide, sequence, location)
    return None


UMB_pattern = re.compile(r"UMB\d{4}")


def process_file(file):
    """Extracts unique proteins and site indexes from a .pepXML file."""
    from pyteomics import pepxml

    UMB = UMB_pattern.search(os.path.basename(file)).group()
    species, assembly = UMB_to_GCA.get(UMB, 'NA')
    ka = os.path.join(genebank_dir, assembly + ".kegg.txt")
    orthologs = dict()
    with open(ka, "r") as inFile:
        for l in inFile:
            pair = l.split()
            if len(pair) == 2:
                orthologs[pair[0]] = pair[1]

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
<<<<<<< Updated upstream
                        ortho_org[o] = defaultdict(set)
                    ortho_org[o][protein+", "+assembly].add(mod_index)
=======
                        ortho_org[o] = dict()
                    if protein + ", " + assembly not in ortho_org[o]:
                        ortho_org[o][protein + ", " + assembly + ", " + species] = {'read_start': dict(),
                                                                                    'read_end': dict(),
                                                                                    'mod_site': set()}
                    ortho_org[o][protein + ", " + assembly + ", " + species]['read_start'][startIdx] = \
                    ortho_org[o][protein + ", " + assembly + ", " + species]['read_start'].get(startIdx, 0) + 1
                    ortho_org[o][protein + ", " + assembly + ", " + species]['read_end'][endIdx] = \
                    ortho_org[o][protein + ", " + assembly + ", " + species]['read_end'].get(endIdx, 0) + 1
                    ortho_org[o][protein + ", " + assembly + ", " + species]['mod_site'].update(mod_indices)
>>>>>>> Stashed changes
    except Exception as e:
        print(f"Error processing {file}: {e}")

    return ortho_org


protein_ids = set()
subsequence_indexes = dict()

if __name__ == "__main__":
<<<<<<< Updated upstream
    import multiprocessing
    import json

=======
>>>>>>> Stashed changes
    ctx = multiprocessing.get_context("spawn")

    with ctx.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = pool.imap_unordered(process_file, modification_files)
        for index in results:
            for key, value in index.items():
                if key not in subsequence_indexes:
                    subsequence_indexes[key] = defaultdict(set)
                for k, v in value.items():
<<<<<<< Updated upstream
                    subsequence_indexes[key][k].update(v)
    for kid, a in subsequence_indexes.items():
        indexes = {i:list(j) for i,j in a.items()}
        with open(os.path.join(outdir,kid+outsuffix+".json"), "w") as file:
            json.dump(indexes, file, indent=4)
=======
                    if k not in subsequence_indexes[key]:
                        subsequence_indexes[key][k] = {'read_start': dict(), 'read_end': dict(), 'mod_site': set()}

                    for ki, vi in v['read_start'].items():
                        if ki not in subsequence_indexes[key][k]['read_start']:
                            subsequence_indexes[key][k]['read_start'][ki] = 0
                        subsequence_indexes[key][k]['read_start'][ki] += vi

                    for ki, vi in v['read_end'].items():
                        if ki not in subsequence_indexes[key][k]['read_end']:
                            subsequence_indexes[key][k]['read_end'][ki] = 0
                        subsequence_indexes[key][k]['read_end'][ki] += vi

                    subsequence_indexes[key][k]['mod_site'].update(v['mod_site'])
    for kid, a in subsequence_indexes.items():
        for i in a.keys():
            a[i]["mod_site"] = list(a[i]["mod_site"])
        # indexes = {i:list(j) for i,j in a.items()}
        with open(os.path.join(outdir, kid + outsuffix + ".json"), "w") as file:
            json.dump(a, file, indent=4)
>>>>>>> Stashed changes
