# PTMOverlay
Capstone PTM alignment

## Dependencies
This workflow uses the python packages pyteomics and snakemake
```
pip install pyteomics biopython lxml snakemake ncbi-taxonomist
```
#### Tested Package versions:
- pyteomics: 4.7.5
- biopython: 1.85
- lxml: 5.3.1
- snakemake: 9.1.1
- ncbi-taxonomist: 1.2.1

These versions are not the only ones that will work, but are the ones that were used by the authors.

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
The proteome and kegg annotation files are included as an example. Until we graduate and the files are no longer hosted on BYU Box, the proteome files can be downloaded and installed from Box automatically when running snakemake:
```
snakemake
```
We recommend keeping the config file the way it is for the first run. If there are other orthologs or pathways you want to look at on the 31 species, rerunning with modified parameters will run faster if intermediates are already generated.

## Figure regeneration (Capstone Reproducibility)
To generate the figures used in BIO 465: Bioinformatics Capstone, run
```
snakemake figures
```
The figures on our poster were as follows:
1. PTMOverlay flowchart (made by hand)
2. Enolase PTMs conservation on select intervals (generated thru code and edited in Adobe Illustrator)
3. Conserved PTM sites on enolase structure (made by Youngki You at PNNL)
4. Significant PTM sites in glycolysis pathway (generated thru code and edited in Adobe Illustrator)

The pre-Adobe Illustrator components of figures 2 and 4 will be in the `figures` folder after running.
