from bs4 import BeautifulSoup
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import matplotlib.pyplot as plt
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.patches as mpatches

GCA_to_name = {"GCA_002847685.2":"Dermabacter hominis", 'GCA_002847765.1':'Alloscardovia omnicolens', 
               'GCA_002847825.1':'Kytococcus schroeteri', 'GCA_002861815.1':'Lactobacillus crispatus', 
               'GCA_002861965.1':'Gardnerella vaginalis', 'GCA_002871555.1':'Limosilactobacillus portuensis', 
               'GCA_002871635.1':'Gardnerella leopoldii', 'GCA_002871815.1':'Bifidobacterium breve', 
               'GCA_002884585.1':'Streptococcus agalactiae', 'GCA_002884815.1':'Gardnerella piotii', 
               'GCA_002884955.2':'Aerococcus loyolae', 'GCA_002940945.1':'Limosilactobacillus pontis', 
               'GCA_003286645.3':'Aerococcus tenax', 'GCA_003286795.1':'Aerococcus urinae', 
               'GCA_007785995.1':'Lactobacillus gasseri', 'GCA_008726885.1':'Aerococcus mictus', 
               'GCA_008727025.1':'Lactobacillus jensenii', 'GCA_008727075.1':'Streptococcus urinae', 
               'GCA_008728115.1':'Lactobacillus mulieris', 'GCA_030218505.1':'Lactobacillus gasseri', 
               'GCA_030218525.1':'Lactobacillus paragasseri', 'GCA_030218565.1':'Streptococcus oralis', 
               'GCA_030218815.1':'Streptococcus mitis', 'GCA_030218925.1':'Bifidobacterium dentium', 
               'GCA_030222485.2':'Globicatella sanguinis', 'GCA_030225245.1':'Staphylococcus warneri', 
               'GCA_030225265.1':'Staphylococcus epidermidis', 'GCA_030226495.1':'Lactobacillus iners', 
               'GCA_030226535.1':'Lactobacillus delbrueckii', 'GCA_030230945.1':'Enterococcus faecalis', 
               'GCA_019890915.1':'Micrococcus luteus', 'GCA_003286795.2':'Aerococcus urinae', 
               'GCA_002884955.3':'Aerococcus loyolae'
               }

index_to_GCA = {1: "GCA_002861815.1", 2: "GCA_030226535.1", 3: "GCA_030226495.1", 4: "GCA_002861815.1",
    5: "GCA_007785995.1", 6: "GCA_030218505.1", 7: "GCA_030218525.1", 8: "GCA_002871635.1", 9: "GCA_002861965.1",
    10: "GCA_002884815.1", 11: "GCA_002847765.1", 12: "GCA_002871815.1", 13: "GCA_030218925.1", 14: "GCA_019890915.1",
    15: "GCA_002847685.2", 16: "GCA_002847825.1", 17: "GCA_030225265.1", 18: "GCA_030225245.1", 19: "GCA_030230945.1",
    20: "GCA_002884585.1", 21: "GCA_030218565.1", 22: "GCA_030218815.1", 23: "GCA_008727075.1", 24: "GCA_002871555.1",
    25: "GCA_002940945.1", 26: "GCA_003286795.2", 27: "GCA_002884955.3", 28: "GCA_003286645.3", 29: "GCA_008726885.1",
    30: "GCA_030222485.2", 31: "GCA_008728115.1", 32: "GCA_008727025.1", 33: "GCA_030218505.1", 34: "GCA_030218525.1",
    35: "GCA_007785995.1", 36: "GCA_030226495.1"
}

#html_file = file path to html file with Multiple Sequence Alignment already completed
with open(html_file, "r", encoding="utf-8") as f:
    soup = BeautifulSoup(f, "html.parser")

# ---- Extract species names and sequences ---- #
species_data = {}
modifications = {}

table_rows = soup.find_all("table")[1].find_all("tr")[1:]

for row in table_rows:
    cols = row.find_all("td")
    if len(cols) < 2:
        continue
    species_id = cols[0].text.strip()
    sequence_tag = cols[1].find("pre")
    if not sequence_tag:
        continue  
    sequence = sequence_tag.get_text(strip=True)
    mod_positions = {}
    index = 0
    for span in cols[1].find_all("span"):
        residue = span.text.strip()
        if residue:
            for mod in ["AcetylCarbamyl", "Phospho", "MonoMethyl", "DiMethyl", "TriMethyl"]:
                if mod in span.get("class", []):
                    mod_positions[index] = mod
            index += 1 
    species_data[species_id] = sequence
    modifications[species_id] = mod_positions

# ---- Convert to Biopython format ---- #
alignment_records = [
    SeqRecord(Seq(seq.replace("…", ".")), id=name, description="")
    for name, seq in species_data.items()
]

fasta_filename = "extracted_sequences.fasta"
with open(fasta_filename, "w") as fasta_file:
    SeqIO.write(alignment_records, fasta_file, "fasta")

# ---- Read Multiple Sequence Alignment ---- #
alignment = AlignIO.read(fasta_filename, "fasta")

# ---- Compute and Plot Phylogenetic Tree ---- #
calculator = DistanceCalculator("identity")
dm = calculator.get_distance(alignment)
constructor = DistanceTreeConstructor()
tree = constructor.upgma(dm)
Phylo.write(tree, "generated_tree.nh", "newick")
tree = Phylo.read("generated_tree.nh", "newick")

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
ordered_alignment = sorted(alignment, key=lambda rec: leaf_order.index(rec.id))
ordered_modifications = {species: modifications[species] for species in leaf_order}

# ---- Multiple Sequence Alignment Visualization ---- #
mod_colors = {
    "AcetylCarbamyl": "purple",
    "Phospho": "red",
    "MonoMethyl": "blue",
    "DiMethyl": "green",
    "TriMethyl": "orange"
}
default_color = "black"
seq_matrix = np.array([list(rec.seq) for rec in ordered_alignment])
empty_seq = np.array([" "] * len(alignment[0].seq))
seq_matrix = np.vstack([empty_seq, seq_matrix, empty_seq])
leaf_order = [" "] + leaf_order + [" "]

for i, (species, seq) in enumerate(zip(leaf_order, seq_matrix)):
    if species not in ordered_modifications:
        continue
    for j, letter in enumerate(seq):
        if letter == ".":
            letter = " … "
        mod_type = ordered_modifications[species].get(j, None)
        if mod_type:
            color = mod_colors[mod_type]
            axes[1].add_patch(plt.Rectangle((j - 0.4, i - 0.3), 0.8, 0.8, color=color, alpha=1.0))
            text_color = "white" if color in mod_colors.values() else "black"
        else:
            text_color = "black"

        axes[1].text(j, i, letter, ha="center", va="center", color=text_color, fontsize=8, fontweight="bold")

axes[1].set_yticks(range(len(leaf_order)))
leaf_order = [
    index_to_GCA.get(int(species), species) if species.strip().isdigit() else species
    for species in leaf_order
]
leaf_order = [GCA_to_name.get(name, name) for name in leaf_order]
axes[1].set_yticklabels(leaf_order, fontsize=10, fontstyle='italic')
axes[1].tick_params(axis="x", which="both", bottom=False, top=False)
axes[1].tick_params(axis="y", which="both", left=False, right=False)
axes[1].set_xlim(-1, max(map(len, seq_matrix)))

custom_labels = {
    "AcetylCarbamyl": "Acetylation",
    "Phospho": "Phosphorylation",
    "MonoMethyl": "Monomethylation",
    "DiMethyl": "Dimethylation",
    "TriMethyl": "Trimethylation"
}

legend_patches = [mpatches.Patch(color=color, label=custom_labels[mod]) for mod, color in mod_colors.items()]
axes[1].legend(handles=legend_patches, title="Modifications", loc="upper left", bbox_to_anchor=(1.05, 1), fontsize=10)

plt.tight_layout()
# plt.savefig("alignment_plot.pdf", bbox_inches="tight")
plt.show()
