from pathlib import Path
import os
import json
import glob
import sys
import requests
import pandas as pd


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
SPECIES_INFO=config["species_info"]

FINAL_ALIGNMENTS=config["final_alignment_dir"]

PTM_DIR="intermediates/ptm"
PRE_ALIGN_FASTAS="intermediates/preAlign"
RAW_ALIGNMENTS="intermediates/rawAlign"
TREE_DIR="intermediates/tree"
TAXON_TREE_DIR="intermediates/taxon_tree"
if "batch_name" in config:
  BATCH_PREFIX = config["batch_name"]+"_"
else:
  BATCH_PREFIX = ""

#---------STRAIN/ASSEMBLY MATCHUP---------------
org_table = pd.read_csv("index_umb_taxa_gca.tsv", sep="\t", index_col="Inde\
x")
strain_to_assembly = dict(zip(org_table["UMB"],org_table["Assembly"]))
strain_to_species = dict(zip(org_table["UMB"],org_table["Taxa"]))

#----------HANDLING FOR "ALL" KO'S----------
kofiles = glob.glob(os.path.join(PROTEOMES,"*.kegg.txt"))
ORTHOLOGS = dict()
CUTOFF = 2 if "ortholog_sequence_cutoff" not in config else max(config["ortholog_sequence_cutoff"],2)
MIN_PTMS_CUTOFF = .2 if "min_ptms_freq" not in config else config["min_ptms_freq"]
for k in kofiles:
  with open(k,"r") as inFile:
    for l in inFile:
      pair = l.split()
      if len(pair)==2:
        ORTHOLOGS[pair[1]] = ORTHOLOGS.get(pair[1],0)+1
if RORTHOLOGS=="ALL":
  ORTHOLOGS = [i for i,j in ORTHOLOGS.items() if j >= CUTOFF] 
else:
  ORTHOLOGS = {i for i,j in ORTHOLOGS.items() if j >= CUTOFF}
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
    pdfs=expand(FINAL_ALIGNMENTS+'/'+"_".join(PTM_TYPES)+'/{ko}__{symbol}__{name}.tree.pdf',zip, ko=ORTHOLOGS, symbol=SYMBOLS, name=NAMES),
    csv=FINAL_ALIGNMENTS+'/'+BATCH_PREFIX+'filtered_ptms.csv',
#    svg=FINAL_ALIGNMENTS+'/'+BATCH_PREFIX+'filtered_ptms.svg',
    taxonomic_tree=FINAL_ALIGNMENTS+'/Taxonomy_Tree.pdf'

rule plotCSV:
  input:
    csv=FINAL_ALIGNMENTS+'/'+BATCH_PREFIX+'filtered_ptms.csv'
  output:
    svg=FINAL_ALIGNMENTS+'/'+BATCH_PREFIX+'filtered_ptms.svg'
  shell:
    '''
    {PYTHON} scripts/plot_csv_table.py {input.csv} {output.svg}
    '''

rule filterCSV:
  input:
    csv=FINAL_ALIGNMENTS+'/'+BATCH_PREFIX+'unfiltered_ptms.csv'
  output:
    csv=FINAL_ALIGNMENTS+'/'+BATCH_PREFIX+'filtered_ptms.csv'
  shell:
    '''
    {PYTHON} scripts/filter_ptm_table.py {input.csv} {output.csv} {MIN_PTMS_CUTOFF}
    '''

rule makeBigCSV:
  input:
    ptm=expand(PTM_DIR+'/{k}__aligned.json',k=ORTHOLOGS),
    faa=expand(RAW_ALIGNMENTS+'/{k}.faa',k=ORTHOLOGS)
  output:
    csv=FINAL_ALIGNMENTS+'/'+BATCH_PREFIX+'unfiltered_ptms.csv'
  params:
    ptm="@@".join(PTM_TYPES),
    ko="@@".join(ORTHOLOGS),
    kn="@@".join(SYMBOLS)
  shell:
    '''
    {PYTHON} scripts/ptm_table.py {RAW_ALIGNMENTS} {PTM_DIR} '{params.ptm}' '{params.ko}' '{params.kn}' {output.csv}
    '''

rule alignPTMs:
  input:
    fasta=RAW_ALIGNMENTS+'/{ko}.faa',
    ptms=PTM_DIR+'/{ko}__ortholog.json'
  output:
    ptms=PTM_DIR+'/{ko}__aligned.json'
  shell:
    '''
    {PYTHON} scripts/ptm_liftover.py {input.ptms} {input.fasta} {output.ptms}
    '''

rule fastaAnnotate:
  input:
    alignment=RAW_ALIGNMENTS+'/{ko}.faa',
    ptm=PTM_DIR+'/{ko}__aligned.json'
  output:
    html=FINAL_ALIGNMENTS+'/{ptm_types}/{ko}__{symbol}__{name}.html'
  params:
    title=lambda w: FULLNAMES[w.ko]
  shell:
    '''
    {PYTHON} scripts/reFormatFasta.py {input.alignment} {output.html} '{params.title}' {input.ptm}
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

rule gather_ptms:
  input:
    ptm_jsons=expand(PTM_DIR+'/{strain_name}__strain.json', strain_name=strain_to_assembly.keys())
  output:
    ptm_json=PTM_DIR+'/{ko}__ortholog.json'
  shell:
    '''
    {PYTHON} scripts/gather_ptm.py {input.ptm_jsons} {wildcards.ko} {output.ptm_json}
    '''

rule split_ptms:
  input:
    pepXML_dir=ancient(PEPXML_DIR+"/{strain_name}"),
    fasta=ancient(expand(PROTEOMES+'/{proteome_name}', proteome_name=lambda w: strain_to_assembly[w.strain_name])),
    mass=ancient('scripts/ptm_mass.yaml')
  output:
    ptms=PTM_DIR+'/{strain_name}__strain.json'
  params:
    species=lambda w: strain_to_species[w.strain_name],
  resources:
    mem_gb=3
  shell:
    '''
    {PYTHON} scripts/split_pepXML.py {input.pepXML_dir} {input.fasta} {input.mass} '{params.species}' {PTM_DIR} {wildcards.strain_name}
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
    pdf=FINAL_ALIGNMENTS+'/{ptm_types}/{ko}__{symbol}__{name}.tree.pdf'
  shell:
    '''
    {PYTHON} scripts/generateTree.py {input.html} {input.fasta} {input.nh} {input.json} {input.tsv} {output.pdf}
    '''

rule download_proteome_ftp:
  input:
    cred="ftp_credentials.yaml"
  output:
    proteome=PROTEOMES+"/{proteome}"
  params:
    path="sequence/FASTA_Files/{proteome}"
  shell:
    '''
    {PYTHON} scripts/download_ftp.py {input.cred} {params.path} {output.proteome}
    '''

rule download_mass_spec_ftp:
  input:
    cred="ftp_credentials.yaml"
  output:
    proteome=directory(PEPXML_DIR+"/{strain}")
  params:
    path="other/pepXML_Files/*{strain}*"
  shell:
    '''
    {PYTHON} scripts/download_ftp.py {input.cred} {params.path} {output.proteome}
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
    {PYTHON} scripts/extract_example_data.py
    '''

######## Code from Matthew Cloward for generating proteome phylogenies ###############


rule get_species_ids:
  input:
    info = SPECIES_INFO
  output:
    ids = TAXON_TREE_DIR+'/IDs.txt'
  shell:
    '''
    {PYTHON} scripts/extract_species_ID.py {output.ids} {input.info}
    '''

rule get_taxonomy_objects_from_ids:
  input:
    ids = TAXON_TREE_DIR+'/IDs.txt'
  output:
    taxons = TAXON_TREE_DIR+'/TaxonomyObjs.txt'
  shell:
    '''
    bash scripts/getTaxonomyObjsFromAccessionIDs.sh {input.ids} {output.taxons}
    '''

rule get_taxons_from_objects:
  input:
    ids = TAXON_TREE_DIR+'/IDs.txt',
    taxons = TAXON_TREE_DIR+'/TaxonomyObjs.txt'
  output:
    tax = TAXON_TREE_DIR+'/Taxonomy.tsv'
  shell:
    '''
    {PYTHON} scripts/getTaxonomiesFromTaxonomyObjs.py {input.ids} {input.taxons} {output.tax}
    '''

rule convert_taxonomy_to_tree:
  input:
    tax = TAXON_TREE_DIR+'/Taxonomy.tsv'
  output:
    tree_json = TAXON_TREE_DIR+'/Tree.json'
  shell:
    '''
    {PYTHON} scripts/convertTaxononmyFileToTree.py {input.tax} {output.tree_json}
    '''

rule tree_to_nh:
  input:
    tree_json = TAXON_TREE_DIR+'/Tree.json'
  output:
    tree_nh = TAXON_TREE_DIR+'/Tree.tre'
  shell:
    '''
    {PYTHON} scripts/jsonToNewick.py {input.tree_json} {output.tree_nh}
    '''

rule generate_taxon_tree:
  input:
    tree_nh = TAXON_TREE_DIR+'/Tree.tre'
  output:
    final_tree = FINAL_ALIGNMENTS+'/Taxonomy_Tree.pdf'
  shell:
    '''
    {PYTHON} scripts/generateTaxonTree.py {input.tree_nh} {output.final_tree}
    '''

  
########## Hardcoded Rules for reproducibility #################

rule figures:
  input:
    pdf=FINAL_ALIGNMENTS+'/Phospho_Acetyl_MonoMethyl_DiMethyl_TriMethyl/K01689__ENO1_2_3__enolase_1_2_3.tree.pdf',
    csv=FINAL_ALIGNMENTS+'/'+BATCH_PREFIX+'filtered_ptms.csv',
    tree=FINAL_ALIGNMENTS+'/Taxonomy_Tree.pdf'
  output:
    pdf='figures/figure2_Enolase.pdf',
    svg='figures/figure4_Glycolysis_ptms.svg',
    tree='figures/figure4_Glycolysis_tree.pdf',
    combined='figures/figure4_combined.svg'
  params:
    csv=FINAL_ALIGNMENTS+'/'+BATCH_PREFIX+'filtered_organized_ptms.csv',
    tree_svg=TAXON_TREE_DIR+'/taxon_tree.svg',
    scaled_tree=TAXON_TREE_DIR+'/scaled_taxon_tree.svg'
  shell:
    '''
    {PYTHON} scripts/gatherFigures.py {input.pdf} {input.csv} {input.tree} {output.pdf} {params.csv} {output.tree}
    {PYTHON} scripts/plot_csv_table.py {params.csv} {output.svg}
    pdf2svg {input.tree} {params.tree_svg}
    {PYTHON} scripts/scale_svg.py {params.tree_svg} {params.scaled_tree}
    {PYTHON} scripts/combine_taxon_pathway.py {params.scaled_tree} {output.svg} {output.combined}
    '''

