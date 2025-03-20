from pathlib import Path
import os
import json
import glob
import sys
import requests

#---------------CONFIG--------------------
configfile: "config.yaml"
PYTHON = config["python"]

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
TREE_DIR=config["tree_intermediates"]
TREE_ALIGN_DIR=config["tree_final_alignments"]
SPECIES_INFO=config["species_info"]
FINAL_ALIGNMENTS=config["final_alignment_dir"]

#----------HANDLING FOR "ALL" KO'S----------
kofiles = glob.glob(os.path.join(PROTEOMES,"*.kegg.txt"))
ORTHOLOGS = dict()
CUTOFF = config["cutoff"]
for k in kofiles:
  with open(k,"r") as inFile:
    for l in inFile:
      pair = l.split()
      if len(pair)==2:
        ORTHOLOGS[pair[1]] = ORTHOLOGS.get(pair[1],0)+1
if RORTHOLOGS=="ALL":
  ORTHOLOGS = [i for i,j in ORTHOLOGS.items() if j > CUTOFF] 
else:
  ORTHOLOGS = {i for i,j in ORTHOLOGS.items() if j > CUTOFF}
  for i in set(RORTHOLOGS) - ORTHOLOGS:
    print(i,"will not be aligned due to an insufficent number of sequences (n < "+str(CUTOFF)+")")
  ORTHOLOGS = ORTHOLOGS.intersection(set(RORTHOLOGS))

if type(ORTHOLOGS)!=list and type(ORTHOLOGS)!=set:
  print("ERROR: bad ortholog selection")
  print("Orthologs should be a list of kegg orthology id's (K#####)")
  print("or a minimum number of proteomes that the ortholog is found in (all>15)")
  sys.exit()

#-------LOOKUP KEGG SYMBOLS AND NAMES-------
ORTHOLOGS = list(ORTHOLOGS)
SYMBOLS = []
NAMES = []
FULLNAMES = dict()
SN_MATCHUPS = dict()
for i in ORTHOLOGS:
  url="https://rest.kegg.jp/get/"+i
  koInfo=requests.get(url).text
  koInfo = koInfo.split("\n")
  symbol = koInfo[1][12:]
  name = koInfo[2][12:]
  FULLNAMES[i] = name
  symbol = symbol.split(", ")[0].strip()
  symbol = re.sub("\.","_",symbol)
  name = re.sub("\[.*]","",name)
  name = re.sub("\(.*\) ", "", name)
  name = re.sub("[\(\),]","", name)
  name = re.sub(" /.*$","", name).strip()
  name = re.sub("[/ :]","_", name)
  SYMBOLS.append(symbol)
  NAMES.append(name)
  SN_MATCHUPS[i] = (symbol, name)

#------------------RULES--------------------

wildcard_constraints:
  ko="K[0-9]+",
  ptm_type="[A-Za-z]+",
  ptm_types="[A-Za-z_]+"

rule preAlignBenchmark:
  input:
    html=expand(FINAL_ALIGNMENTS+'/'+"_".join(PTM_TYPES)+'/{ko}__{symbol}__{name}.html',zip, ko=ORTHOLOGS, symbol=SYMBOLS, name=NAMES),
    pdfs=expand(TREE_ALIGN_DIR+'/'+"_".join(PTM_TYPES)+'/{ko}__{symbol}__{name}.tree.pdf',zip, ko=ORTHOLOGS, symbol=SYMBOLS, name=NAMES),
    csv=PTM_DIR+'/raw_ptms.csv'

rule makeBigCSV:
  input:
    ptm=expand(PTM_DIR+'/{k}_{p}_aligned.json',k=ORTHOLOGS,p=PTM_TYPES),
    faa=expand(RAW_ALIGNMENTS+'/{k}.faa',k=ORTHOLOGS)
  output:
    csv=PTM_DIR+'/raw_ptms.csv'
  params:
    ptm=" ".join(PTM_TYPES),
    ko=" ".join(ORTHOLOGS)
  shell:
    '''
    {PYTHON} scripts/make_ptm_csv.py {PTM_DIR} {RAW_ALIGNMENTS} '{params.ptm}' '{params.ko}' {output.csv}
    '''

rule alignPTMs:
  input:
    fasta=RAW_ALIGNMENTS+'/{ko}.faa',
    ptms=PTM_DIR+'/{ko}_{ptm_type}.json'
  output:
    ptms=PTM_DIR+'/{ko}_{ptm_type}_aligned.json'
  shell:
    '''
    {PYTHON} scripts/ptm_liftover.py {input.ptms} {input.fasta} {output.ptms}
    '''

rule fastaAnnotate:
  input:
    alignment=RAW_ALIGNMENTS+'/{ko}.faa',
    ptms=expand(PTM_DIR+'/{{ko}}_{pt}_aligned.json', pt=lambda w: w.ptm_types.split("_"))
  output:
    html=FINAL_ALIGNMENTS+'/{ptm_types}/{ko}__{symbol}__{name}.html'
  params:
    title=lambda w: FULLNAMES[w.ko]
  shell:
    '''
    {PYTHON} scripts/reFormatFasta.py {input.alignment} {output.html} '{params.title}' {input.ptms}
    '''

rule muscle:
  input:
    fasta=PRE_ALIGN_FASTAS+'/{ko}.faa'
  output:
    alignment=RAW_ALIGNMENTS+'/{ko}.faa'
  shell:
    '''
    {PYTHON} scripts/runMUSCLE.py {input.fasta} {output.alignment}
    '''

rule extract_ptms:
  input:
    pepXML_dir=PEPXML_DIR+'/{ptm_type}',
    mass='scripts/ptm_mass.yaml'
  output:
    ptms=expand(PTM_DIR+'/{ko}_{{ptm_type}}.json', ko=ORTHOLOGS)
  shell:
    '''
    {PYTHON} scripts/parse_pepXML.py {input.pepXML_dir} {PROTEOMES} {input.mass} {SPECIES_INFO} {PTM_DIR} {wildcards.ptm_type}
    '''

rule group_orthologs:
  output:
    fastas=PRE_ALIGN_FASTAS+'/{ko}.faa'
  shell:
    '''
    {PYTHON} scripts/group_orthologs.py {PROTEOMES} {wildcards.ko} {SPECIES_INFO} {PRE_ALIGN_FASTAS}
    '''    

rule generate_tree_fasta:
  input:
    html = expand(FINAL_ALIGNMENTS+'/{{ptm_types}}/{{ko}}__{symbol}__{name}.html', symbol=lambda w: SN_MATCHUPS[w.ko][0], name=lambda w: SN_MATCHUPS[w.ko][1])
  output:
    fasta=TREE_DIR+'/{ko}__{ptm_types}.faa',
    json=TREE_DIR+'/{ko}__{ptm_types}.json'
  params:
    name=lambda w: SN_MATCHUPS[w.ko][1],
    symbol=lambda w: SN_MATCHUPS[w.ko][0]
  shell:
    '''
    {PYTHON} scripts/generateTreeFasta.py {input.html} {output.fasta} {output.json}
    '''

rule generate_tree_file:
  input:
    fasta=TREE_DIR+'/{ko}__{ptm_types}.faa'
  output:
    nh=TREE_DIR+'/{ko}__{ptm_types}.nh'
  shell:
    '''
    {PYTHON} scripts/generateTreeNH.py {input.fasta} {output.nh}
    '''

rule generate_tree:
  input:
    html = FINAL_ALIGNMENTS+'/{ptm_types}/{ko}__{symbol}__{name}.html',
    fasta=TREE_DIR+'/{ko}__{ptm_types}.faa',
    nh=TREE_DIR+'/{ko}__{ptm_types}.nh',
    json=TREE_DIR+'/{ko}__{ptm_types}.json',
    tsv= SPECIES_INFO
  output:
    pdf=TREE_ALIGN_DIR+'/{ptm_types}/{ko}__{symbol}__{name}.tree.pdf'
  shell:
    '''
    {PYTHON} scripts/generateTree.py {input.html} {input.fasta} {input.nh} {input.json} {input.tsv} {output.pdf}
    '''

rule download_example_data:
  output:
    zip=ancient("mass_spec.zip")
  shell:
    '''
    {PYTHON} scripts/download_example_data.py
    '''

rule extract_example_data:
  output:
    Phospho=directory("mass_spec/Phospho"),
    Acetyl=directory("mass_spec/Acetyl"),
    Carbamyl=directory("mass_spec/Carbamyl"),
    MonoMethyl=directory("mass_spec/MonoMethyl"),
    DiMethyl=directory("mass_spec/DiMethyl"),
    TriMethyl=directory("mass_spec/TriMethyl")
  input:
    zip=ancient("mass_spec.zip")
  shell:
    '''
    {PYTHON} scripts/download_example_data.py
    '''

  
