from pathlib import Path
import os
import json
import glob
import sys

#---------------CONFIG--------------------
configfile: "config.yaml"

ORTHOLOGS=config["orthologs_to_align"]
PROTEOMES=config["proteome_dir"]
PRE_ALIGN_FASTAS=config["pre_align_fasta_dir"]
RAW_ALIGNMENTS=config["raw_alignment_dir"]
MUSCLE='runMUSCLE.py'

#----------HANDLING FOR "ALL" KO'S----------
if type(ORTHOLOGS)==str and OHRTHOLOGS[:4].upper() == "ALL>":
  kofiles = glob.glob(os.path.join(PROTEOMES,"*.kegg.txt"))
  threshold = int(ORTHOLOGS[4:])
  ORTHOLOGS = dict()
  for k in kofiles:
    with open(k,"r") as inFile:
      for l in inFile:
        pair = l.split()
        if len(pair)==2:
          ORTHOLOGS[pair[1]] = ORTHOLOGS.get(pair[1],0)+1
  ORTHOLOGS = [i for i,j in ORTHOLOGS.items() if j > threshold] 
if type(ORTHOLOGS)!=list:
  print("ERROR: bad ortholog selection")
  print("Orthologs should be a list of kegg orthology id's (K#####)")
  print("or a minimum number of proteomes that the ortholog is found in (all>15)")
  sys.exit()
#------------------RULES--------------------

rule preAlignBenchmark:
  input:
    fastas=expand(RAW_ALIGNMENTS+'/{ko}.faa', ko=ORTHOLOGS)

rule muscle:
  input:
    fasta=PRE_ALIGN_FASTAS+'/{ko}.faa'
  output:
    alignment=RAW_ALIGNMENTS+'/{ko}.faa'
  shell:
    '''
    mkdir -p {RAW_ALIGNMENTS}
    python3 {MUSCLE} {input.fasta} {output.alignment}
    '''

rule group_orthologs:
  output:
    fastas=expand(PRE_ALIGN_FASTAS+'/{ko}.faa', ko=ORTHOLOGS)
  shell:
    '''
    python3 group_orthologs.py {PROTEOMES} {PRE_ALIGN_FASTAS}
    '''    
