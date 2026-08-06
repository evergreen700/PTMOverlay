# PTMOverlay
Capstone PTM alignment
This pipeline has been successfully tested on Linux, Windows and some Mac systems. It is known to have issues on ARM-based M4 systems.

## UV Environment (Recommended)
To set up the workflow, first clone the git repository:
```
git clone https://github.com/evergreen700/PTMOverlay
cd PTMOverlay
```
If UV is not installed:
```
pip install uv
```
Create the UV environment:
```
uv sync
```
Run the workflow:
```
uv run snakemake --cores all
```
The final output files will be in the folder `align_reports`

## Native

There are extensive packages to install, so a conda environment is recommended. Instructions for creating a new conda environment can be found here: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

Install dependencies using pip:
```
pip install matplotlib pandas pyteomics biopython lxml snakemake ncbi-taxonomist svgutils six
brew install pdf2svg (Mac)
sudo apt-get install pdf2svg (Linux)
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

### MUSCLE Executable
In the executables directory of the repository, check to see if a MUSCLE executable that works with your operating system is present. If not, go to this website and download the correct executable and place it in in executables directory.

https://drive5.com/muscle/downloads_v3.htm

Then edit the runMUSCLE.py file accordingly.

If on Windows:
- Edit line 13 to point to the correct executable

If on Linux:
- Edit lines 16, 17 or 19, 20 depending on your machine
  
If on Mac: 
- Edit line 23 or 25 depending on your machine

### Execution
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

`index_umb_taxa_gca.tsv` is a tab-separated file that is used to match between mass spec strain IDs (UMB####), species name, and proteome assembly (GCA). If you are using your own mass spec and proteome files, make sure that the names are in this tsv.

To run the workflow on your operating system:
```
cd PTMOverlay
snakemake --cores all
```

## Example run
The kegg annotation files are included as an example. The proteomes and the mass spec search results
are downloaded automatically from the MassIVE ftp server (see `ftp_credentials.yaml`) when running
snakemake: `download_sequence_ftp` fetches one `.fasta` per assembly into `proteome/`, and
`download_mass_spec_ftp` fetches one directory per search into `mass_spec/`, unpacking the `.pepXML`
files out of the downloaded zip archives. Until the data is public, edit the `ftp_credentials.yaml` file with your username and password to access the data before running the workflow. Budget ~10 GB of free disk for the unpacked mass spec data.
We recommend keeping the config file the way it is for the first run. If there are other orthologs or pathways you want to look at on the 31 species, rerunning with modified parameters will run faster if intermediates are already generated.
