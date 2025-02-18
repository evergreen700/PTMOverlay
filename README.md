# PTM
Capstone PTM alignment

## Dependencies:
This workflow uses the python packages pyteomics and snakemake
```
pip install pyteomics snakemake
```

## Execution:
To run the workflow, place proteomes and kegg annotation files in the folder designated as `PROTEOMES` in the snakefile. Then run by running snakemake from the root directory:
```
snakemake
```
