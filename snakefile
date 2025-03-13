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
MUSCLE='runMUSCLE.py'
TREE_DIR="tree"
TREE_ALIGN_DIR="phyloAlign"
SPECIES_INFO=config["species_info"]

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

wildcard_constraints:
  ko="K[0-9]+",
  ptm_type="[A-Za-z]+",
  ptm_types="[A-Za-z_]+"

rule preAlignBenchmark:
  input:
    html=expand(RAW_ALIGNMENTS+'/{ko}__{ptm_type}.html', ko=ORTHOLOGS, ptm_type="_".join(PTM_TYPES)),
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
    python ptm_liftover.py {input.ptms} {input.fasta} {output.ptms}
    '''

rule fastaAnnotate:
  input:
    alignment=RAW_ALIGNMENTS+'/{ko}.faa',
    ptms=expand(PTM_DIR+'/{{ko}}_{pt}_aligned.json', pt=lambda w: w.ptm_types.split("_"))
  output:
    html=RAW_ALIGNMENTS+'/{ko}__{ptm_types}.html'
  shell:
    '''
    python reFormatFasta.py {input.alignment} {output.html} {input.ptms}
    '''

rule muscle:
  input:
    fasta=PRE_ALIGN_FASTAS+'/{ko}.faa'
  output:
    alignment=RAW_ALIGNMENTS+'/{ko}.faa'
  shell:
    '''
    python {MUSCLE} {input.fasta} {output.alignment}
    '''

rule extract_ptms:
  input:
    pepXML_dir=PEPXML_DIR+'/{ptm_type}',
    mass='ptm_mass.yaml'
  output:
    ptms=expand(PTM_DIR+'/{ko}_{{ptm_type}}.json', ko=ORTHOLOGS)
  shell:
    '''
    python parse_pepXML.py {input.pepXML_dir} {PROTEOMES} {input.mass} {PTM_DIR} {wildcards.ptm_type}
    '''

rule group_orthologs:
  output:
    fastas=expand(PRE_ALIGN_FASTAS+'/{ko}.faa', ko=ORTHOLOGS)
  shell:
    '''
    python group_orthologs.py {PROTEOMES} {PRE_ALIGN_FASTAS}
    '''    
rule generate_tree_fasta:
  input:
    html = RAW_ALIGNMENTS+'/{ko}__{ptm_types}.html'
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
    html = RAW_ALIGNMENTS+'/{ko}__{ptm_types}.html',
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
