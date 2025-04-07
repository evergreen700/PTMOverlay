# PTMOverlay
Capstone PTM alignment

## Installation
### Docker (Recommended)
This workflow can be run from docker. To set it up, first clone the git repository:
```
git clone https://github.com/evergreen700/PTMOverlay
```
Then, build the docker image:
```
docker build -t PTMOverlayImg .
```

### Native

Install dependencies using pip:
```
pip install matplotlib pandas pyteomics biopython lxml snakemake ncbi-taxonomist
```
Clone the git repository:
```
git clone https://github.com/evergreen700/PTMOverlay
```

We have not extensively tested all the packages this tool requires. These are the versions we have used:
- pyteomics: 4.7.5
- biopython: 1.85
- lxml: 5.3.1
- snakemake: 9.1.1
- ncbi-taxonomist: 1.2.1

These versions are not the only ones that will work, but are the ones that were used by the authors.

## Execution
To run the workflow, place proteomes and kegg annotation files in the folder designated as `proteome_dir` in the config file. Place .pepXML files in folders sorted by ptm type within the folder designated as `pepXML_dir` in the config file. Below is an example of the file structure:
```
PTMOverlay
+ proteome
| + GCA_002847685.2.faa
| + GCA_002847685.2.kegg.txt
| + ...
+ mass_spec
| + Phospho
| | + BioD_urine_UMB0005_01_12Apr24_Arwen_WBEH-23-02-03.pepXML
| | + ... (other phospho pepXMLs)
| + ... (other PTM types)
+ index_umb_taxa_gca.tsv
+ README.md
+ snakefile
+ ...
```

`index_umb_taxa_gca.tsv` is a tab-separated file that isused to match between mass spec strain IDs (UMB####), species name, and proteome assembly (GCA). If you are using your own mass spec and proteome files, make sure that the names are in this tsv.

### Docker Execution
To run the workflow from the docker image, run
```
cd PTMOverlay
docker run -v ./:/PTMOverlay PTMOverlayImg /bin/bash -c "cd /PTMOverlay && snakemake"
```

### Native Execution
To run the workflow on your operating system:
```
cd PTMOverlay
snakemake
```

## Example run
The proteome and kegg annotation files are included as an example. Until we graduate and the files are no longer hosted on BYU Box, the proteome files can be downloaded and installed from Box automatically when running snakemake.
We recommend keeping the config file the way it is for the first run. If there are other orthologs or pathways you want to look at on the 31 species, rerunning with modified parameters will run faster if intermediates are already generated.
