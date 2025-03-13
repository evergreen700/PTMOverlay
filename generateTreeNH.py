import sys
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

fasta_infile = sys.argv[1]
outfile = sys.argv[2]

alignment = AlignIO.read(fasta_infile, "fasta")
calculator = DistanceCalculator("identity")
dm = calculator.get_distance(alignment)
constructor = DistanceTreeConstructor()
tree = constructor.upgma(dm)
Phylo.write(tree, outfile, "newick")
