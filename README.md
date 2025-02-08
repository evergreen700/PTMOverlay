# PTM
Capstone PTM alignment

## Dependencies:
This workflow uses MUSCLE sequence alignment. Install muscle ([https://drive5.com/muscle/](https://drive5.com/muscle/)) and make sure the path in the snakefile points to the muscle executable.\
This workflow also uses the python packages pyteomics and snakemake
```
pip install pyteomics snakemake
```

## Execution:
To run the workflow, place proteomes and kegg annotation files in the folder designated as `PROTEOMES` in the snakefile. Then run by running snakemake from the root directory:
```
snakemake
```
