from pathlib import Path
import os
import json
import glob
import sys

#---------------CONFIG--------------------
configfile: "config.yaml"

ORTHOLOGS=config["orthologs_to_align"]
PTM_TYPES=config["ptm_types"]
PROTEOMES=config["proteome_dir"]
PEPXML_DIR=config["pepXML_dir"]
PTM_DIR="ptm"
PRE_ALIGN_FASTAS=config["pre_align_fasta_dir"]
RAW_ALIGNMENTS=config["raw_alignment_dir"]
MUSCLE='muscle-linux-x86.v5.3' #TODO: determine muscle based on OS or move to configfile

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
    clws=expand(RAW_ALIGNMENTS+'/{ko}.faa', ko=ORTHOLOGS),
    jsons=expand(PTM_DIR+'/{ptm}.json', ptm=PTM_TYPES)

rule muscle:
  input:
    fasta=PRE_ALIGN_FASTAS+'/{ko}.faa'
  output:
    alignment=RAW_ALIGNMENTS+'/{ko}.faa'
  shell:
    '''
    mkdir -p {RAW_ALIGNMENTS}
    ./{MUSCLE} -align {input.fasta} -output {output.alignment}
    '''

rule extract_ptms:
  input:
    pepXML_dir=PEPXML_DIR+'/{PTM_TYPE}'
  output:
    ptms=PTM_DIR+'/{PTM_TYPE}.json'
  shell:
    '''
    mkdir -p {PTM_DIR}
    python3 parse_pepXML.py {input.pepXML_dir} {PROTEOMES} {output.ptms}
    '''

rule group_orthologs:
  output:
    fastas=expand(PRE_ALIGN_FASTAS+'/{ko}.faa', ko=ORTHOLOGS)
  shell:
    '''
    python3 group_orthologs.py {PROTEOMES} {PRE_ALIGN_FASTAS}
    '''    
