# PTM
Capstone PTM alignment

## Dependencies
This workflow uses the python packages pyteomics and snakemake
```
pip install pyteomics biopython snakemake
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
+ README.md
+ snakefile
+ ...
```
Then run by running snakemake from the root directory:
```
snakemake
```
