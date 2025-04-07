from Bio import Phylo
import sys
import matplotlib.pyplot as plt
import re

infile = sys.argv[1]
outfile = sys.argv[2]

os.makedirs(os.path.dirname(outfile), exist_ok=True)

species_dict = {}

def extract_species_names_from_tre(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            matches = re.findall(r'([A-Za-z]+(?:\s[A-Za-z]+)?)', line)
            for match in matches:
                name_parts = match.split()
                if len(name_parts) > 1:
                    second_word = name_parts[1]
                    species_dict[second_word] = match
                else:
                    species_dict[name_parts[0]] = match

extract_species_names_from_tre(infile)

tree = Phylo.read(infile, "newick")

def format_label(clade):
    if clade.is_terminal():
        if clade.name:
            full_species_name = species_dict.get(clade.name, clade.name)
            full_species_name = full_species_name.replace(" ", "~")
            return r"$\mathit{" + full_species_name + "}$"
    return None

fig, ax = plt.subplots(figsize=(12, 10))
Phylo.draw(tree, do_show=False, axes=ax, label_func=format_label)

plt.tight_layout()
plt.savefig(outfile, format="pdf", dpi=300, bbox_inches="tight")
plt.close()
