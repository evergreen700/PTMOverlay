import sys
import os
import shutil
import pandas as pd

#This script tidys up outputs used in our figures and makes them look a bit nicer

inPDF = sys.argv[1]
inCSV = sys.argv[2]

outPDF = sys.argv[3]
outCSV = sys.argv[4]

os.makedirs('figures',exist_ok=True)

#Fig2: just copy over the file into figures directory
shutil.copy(inPDF, outPDF)

#Fig4-csv: reorder the rows according to the tree and the columns according to the pathway order
colOrder = ["GAPDH","PGK","ENO1_2_3","pfkA","adhE","GPI","PK","FBA","TPI","glk"]
rowOrder = """GCA_002847685.2, Dermabacter hominis
GCA_019890915.1, Micrococcus luteus
GCA_002847825.1, Kytococcus schroeteri
GCA_002847765.1, Alloscardovia omnicolens
GCA_030218925.1, Bifidobacterium dentium
GCA_002871815.1, Bifidobacterium breve
GCA_002871635.1, Gardnerella leopoldii
GCA_002861965.1, Gardnerella vaginalis
GCA_002884815.1, Gardnerella piotii
GCA_030225265.1, Staphylococcus epidermidis
GCA_030225245.1, Staphylococcus warneri
GCA_030230945.1, Enterococcus faecalis
GCA_030218815.1, Streptococcus mitis
GCA_030218565.1, Streptococcus oralis
GCA_008727075.1, Streptococcus anginosus
GCA_002884585.1, Streptococcus agalactiae
GCA_030222485.2, Globicatella sanguinis
GCA_008726885.1, Aerococcus mictus
GCA_003286795.2, Aerococcus urinae
GCA_003286645.3, Aerococcus tenax
GCA_002884955.3, Aerococcus loyolae
GCA_002940945.1, Limosilactobacillus pontis
GCA_030226535.1, Lactobacillus delbrueckii
GCA_030226495.1, Lactobacillus iners
GCA_030218525.1, Lactobacillus paragasseri
GCA_008728115.1, Lactobacillus mulieris
GCA_008727025.1, Lactobacillus jensenii
GCA_007785995.1, Lactobacillus gasseri
GCA_002861815.1, Lactobacillus crispatus""".split("\n")

ptmTable = pd.read_csv(inCSV, header=[0,1], index_col=0, dtype=str)
ptmTable = ptmTable.reindex(rowOrder)
ptmTable = ptmTable[colOrder]
ptmTable.to_csv(outCSV)

#Fig4-tree/pdf: copy over
