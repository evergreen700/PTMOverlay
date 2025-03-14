from pathlib import Path
import os
import json
import glob
import sys
import requests

#---------------CONFIG--------------------
configfile: "config.yaml"

#-----------ORTHOLOGS/PATHWAYS------------
if "orthologs_to_align" in config:
	RORTHOLOGS = config["orthologs_to_align"]
else:
	RORTHOLOGS = []
if "pathways" in config:
	for p in config["pathways"]:
		url="https://rest.kegg.jp/link/ko/"+p
		kos=requests.get(url).text
		RORTHOLOGS+=[i.split("\t")[-1][3:] for i in kos.split("\n")[:-1]]
#----------------PTM TYPES---------------
PTM_TYPES=config["ptm_types"].keys()

#-----------DIR SETUP------------------
PROTEOMES=config["proteome_dir"]
PEPXML_DIR=config["pepXML_dir"]
PTM_DIR=config["ptm_dir"]
PRE_ALIGN_FASTAS=config["pre_align_fasta_dir"]
RAW_ALIGNMENTS=config["raw_alignment_dir"]
MUSCLE='runMUSCLE.py'
TREE_DIR="tree"
TREE_ALIGN_DIR="phyloAlign"
SPECIES_INFO=config["species_info"]
FINAL_ALIGNMENTS=config["final_alignment_dir"]

#----------HANDLING FOR "ALL" KO'S----------
kofiles = glob.glob(os.path.join(PROTEOMES,"*.kegg.txt"))
ORTHOLOGS = dict()
for k in kofiles:
  with open(k,"r") as inFile:
    for l in inFile:
      pair = l.split()
      if len(pair)==2:
        ORTHOLOGS[pair[1]] = ORTHOLOGS.get(pair[1],0)+1
if type(RORTHOLOGS)==str and ROHRTHOLOGS[:4].upper() == "ALL>":
  threshold = int(ORTHOLOGS[4:])
  ORTHOLOGS = [i for i,j in ORTHOLOGS.items() if j > threshold] 
else:
  ORTHOLOGS = {i for i,j in ORTHOLOGS.items() if j > 2}
  for i in set(RORTHOLOGS) - ORTHOLOGS:
    print(i,"will not be aligned due to an insufficent number of sequences (n < 2)")
  ORTHOLOGS = ORTHOLOGS.intersection(set(RORTHOLOGS))

if type(ORTHOLOGS)!=list and type(ORTHOLOGS)!=set:
  print("ERROR: bad ortholog selection")
  print("Orthologs should be a list of kegg orthology id's (K#####)")
  print("or a minimum number of proteomes that the ortholog is found in (all>15)")
  sys.exit()
#------------------RULES--------------------

wildcard_constraints:
  ko="K[0-9]+",
  ptm_type="[A-Za-z]+",
  ptm_types="[A-Za-z_]+"

rule preAlignBenchmark:
  input:
    html=expand(FINAL_ALIGNMENTS+'/{ko}__{ptm_type}.html', ko=ORTHOLOGS, ptm_type="_".join(PTM_TYPES)),
    jsons=expand(PTM_DIR+'/{ko}_{ptm_type}_aligned.json', ko=ORTHOLOGS, ptm_type=PTM_TYPES),
    pdfs=expand(TREE_ALIGN_DIR+'/{ko}__{ptm_types}.pdf', ko=ORTHOLOGS, ptm_types="_".join(PTM_TYPES))

rule alignPTMs:
  input:
    fasta=RAW_ALIGNMENTS+'/{ko}.faa',
    ptms=PTM_DIR+'/{ko}_{ptm_type}.json'
  output:
    ptms=PTM_DIR+'/{ko}_{ptm_type}_aligned.json'
  shell:
    '''
    python3 ptm_liftover.py {input.ptms} {input.fasta} {output.ptms}
    '''

rule fastaAnnotate:
  input:
    alignment=RAW_ALIGNMENTS+'/{ko}.faa',
    ptms=expand(PTM_DIR+'/{{ko}}_{pt}_aligned.json', pt=lambda w: w.ptm_types.split("_"))
  output:
    html=FINAL_ALIGNMENTS+'/{ko}__{ptm_types}.html'
  shell:
    '''
    python3 reFormatFasta.py {input.alignment} {output.html} {input.ptms}
    '''

rule muscle:
  input:
    fasta=PRE_ALIGN_FASTAS+'/{ko}.faa'
  output:
    alignment=RAW_ALIGNMENTS+'/{ko}.faa'
  shell:
    '''
    python3 {MUSCLE} {input.fasta} {output.alignment}
    '''

rule extract_ptms:
  input:
    pepXML_dir=PEPXML_DIR+'/{ptm_type}',
    mass='ptm_mass.yaml'
  output:
    ptms=expand(PTM_DIR+'/{ko}_{{ptm_type}}.json', ko=ORTHOLOGS)
  shell:
    '''
    python3 parse_pepXML.py {input.pepXML_dir} {PROTEOMES} {input.mass} {PTM_DIR} {wildcards.ptm_type}
    '''

rule group_orthologs:
  output:
    fastas=PRE_ALIGN_FASTAS+'/{ko}.faa'
  shell:
    '''
    python3 group_orthologs.py {PROTEOMES} {wildcards.ko} {PRE_ALIGN_FASTAS}
    '''    

rule generate_tree_fasta:
  input:
    html = FINAL_ALIGNMENTS+'/{ko}__{ptm_types}.html'
  output:
    fasta=TREE_DIR+'/{ko}__{ptm_types}.faa',
    json=TREE_DIR+'/{ko}__{ptm_types}.json'
  shell:
    '''
    python3 generateTreeFasta.py {input.html} {output.fasta} {output.json}
    '''

rule generate_tree_file:
  input:
    fasta=TREE_DIR+'/{ko}__{ptm_types}.faa'
  output:
    nh=TREE_DIR+'/{ko}__{ptm_types}.nh'
  shell:
    '''
    python3 generateTreeNH.py {input.fasta} {output.nh}
    '''

rule generate_tree:
  input:
    html = FINAL_ALIGNMENTS+'/{ko}__{ptm_types}.html',
    fasta=TREE_DIR+'/{ko}__{ptm_types}.faa',
    nh=TREE_DIR+'/{ko}__{ptm_types}.nh',
    json=TREE_DIR+'/{ko}__{ptm_types}.json',
    tsv= SPECIES_INFO
  output:
    pdf=TREE_ALIGN_DIR+'/{ko}__{ptm_types}.pdf'
  shell:
    '''
    python3 generateTree.py {input.html} {input.fasta} {input.nh} {input.json} {input.tsv} {output.pdf}
    '''
