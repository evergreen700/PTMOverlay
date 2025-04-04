from Bio import Phylo
import sys
import matplotlib.pyplot as plt

infile = sys.argv[1]
outfile = sys.argv[2]

tree = Phylo.read(infile, "newick")

fig, ax = plt.subplots(figsize=(12, 10))
Phylo.draw(tree, do_show=False, axes=ax)

plt.tight_layout()
plt.savefig(outfile, format="pdf", dpi=300, bbox_inches="tight")
plt.close()