# PTM
Capstone PTM alignment

## Dependencies
This workflow uses the python packages pyteomics and snakemake
```
pip install pyteomics biopython lxml snakemake
```

## Execution
To run the workflow, place proteomes and kegg annotation files in the folder designated as `proteome_dir` in the config file. Place .pepXML files in folders sorted by ptm type within the folder designated as `pepXML_dir` in the config file. Below is an example of the file structure:
```
PTM
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
Then run by running snakemake from the root directory:
```
snakemake
```
`index_umb_taxa_gca.tsv` is a tab-separated file that isused to match between mass spec strain IDs (UMB####), species name, and proteome assembly (GCA). If you are using your own mass spec and proteome files, make sure that the names are in this tsv.

## Example run
The proteome and kegg annotation files are included as an example. To install the corresponding mass spec data, run:
```
python3 scripts/download_example_data.py
```
Once the mass spec data is installed, you should have a `mass_spec` directory with 5 subdirectories. Try the tool out by running snakemake:
```
snakemake
```
