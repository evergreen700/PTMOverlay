from Bio import AlignIO
import os
import numpy as np
import matplotlib.pyplot as plt
from Bio import Phylo
import matplotlib.patches as mpatches
import sys
import json
import csv

infile = sys.argv[1] #html file
fasta_infile = sys.argv[2] #fasta with sequences to graph
nh_infile = sys.argv[3] #tree
json_infile = sys.argv[4] #species data and modifications dictionaries
names_file = sys.argv[5] #tsv file
outfile = sys.argv[6]

alignment = AlignIO.read(fasta_infile, "fasta")

filename = os.path.basename(os.path.dirname(infile))
split_title = filename.split("_")
ptm_types = split_title

with open(json_infile, "r") as json_file:
    data = json.load(json_file)
species_data = data["species_data"]
modifications = data["modifications"]
index_to_GCA = data["GCA_codes"]

GCA_to_name = {}

with open(names_file, "r", encoding="utf-8") as file:
    reader = csv.DictReader(file, delimiter="\t")
    for row in reader:
        GCA_to_name[row["Assembly"]] = row["Taxa"]

# ---- Load and Plot Phylogenetic Tree ---- #
tree = Phylo.read(nh_infile, "newick")
fig, axes = plt.subplots(1, 2, figsize=(15, 8), gridspec_kw={"width_ratios": [1, 2]})
for ax in axes:
    for spine in ax.spines.values():
        spine.set_visible(False)
for clade in tree.get_nonterminals():
    clade.name = None
axes[0].set_xticks([])
axes[0].set_yticks([])
Phylo.draw(tree, do_show=False, axes=axes[0])

# ---- Order sequences to match phylogenetic tree ---- #
leaf_order = [leaf.name for leaf in tree.get_terminals()][::-1] 

# Reorder sequences and modifications based on tree order
ordered_alignment = sorted(alignment, key=lambda rec: leaf_order.index(rec.id))
ordered_modifications = {species: modifications[species] for species in leaf_order}

# ---- Multiple Sequence Alignment Visualization ---- #
colors = [
    "red", "blue", "green", "magenta", "orange", "yellow", "purple", "pink", "lime", "teal",
    "cyan", "gold", "brown", "gray", "navy", "olive", "maroon", "turquoise", "indigo", "violet",
    "chartreuse", "salmon", "coral", "skyblue", "lavender", "peru", "orchid", "plum", "seagreen",
    "sienna", "steelblue", "tomato", "wheat", "khaki", "lightgreen", "darkblue", "darkred", "darkcyan",
    "darkmagenta", "darkorchid", "firebrick", "forestgreen", "dodgerblue", "goldenrod"
]
mod_colors = {label: colors[i % len(colors)] for i, label in enumerate(ptm_types)}
default_color = "black"

# Convert ordered alignment to a sequence matrix
seq_matrix = np.array([list(rec.seq) for rec in ordered_alignment])

# Insert blank lines at the top and bottom
empty_seq = np.array([" "] * len(alignment[0].seq))
seq_matrix = np.vstack([empty_seq, seq_matrix, empty_seq])
leaf_order = [" "] + leaf_order + [" "]

for i, (species, seq) in enumerate(zip(leaf_order, seq_matrix)):
    if species not in ordered_modifications:
        continue
    for j, letter in enumerate(seq):
        if letter == ".":
            letter = " … "
        mod_type = ordered_modifications[species].get(str(j), None)
        if mod_type:
            color = mod_colors[mod_type]
            axes[1].add_patch(plt.Rectangle((j - 0.4, i - 0.3), 0.8, 0.8, color=color, alpha=1.0))
            text_color = "white" if color in mod_colors.values() else "black"
        else:
            text_color = "black"

        axes[1].text(j, i, letter, ha="center", va="center", color=text_color, fontsize=8, fontweight="bold")

# Label axes to match tree order
axes[1].set_yticks(range(len(leaf_order)))
leaf_order = [
    index_to_GCA.get(species, species) if species.strip().isdigit() else species
    for species in leaf_order
]
leaf_order = [GCA_to_name.get(name, name) for name in leaf_order]
axes[1].set_yticklabels(leaf_order, fontsize=10, fontstyle='italic')
axes[1].tick_params(axis="x", which="both", bottom=False, top=False)
axes[1].tick_params(axis="y", which="both", left=False, right=False)

axes[1].set_xlim(-1, max(map(len, seq_matrix)))

legend_patches = [mpatches.Patch(color=color, label=mod) for mod, color in mod_colors.items()]
axes[1].legend(handles=legend_patches, title="Modifications", loc="upper left", bbox_to_anchor=(1.05, 1), fontsize=10)

plt.tight_layout()
plt.savefig(outfile, bbox_inches="tight")
