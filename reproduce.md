## Figure regeneration (Capstone Reproducibility)
This can be reproduced on Mac and Linux machines.\
To generate the figures used in BIO 465: Bioinformatics Capstone, follow these steps: 
1. Install docker from https://docs.docker.com
2. Clone this repository by running this in your terminal:
   ```
   git clone https://github.com/evergreen700/PTMOverlay
   ```
3. Navigate inside of this repository and build the docker image:
   ```
   cd PTMOverlay
   docker build -t ptm-overlay .
   ```
4. Run snakemake from inside the docker image:  
  ```
  docker run -v ./:/PTMOverlay ptm-overlay /bin/bash -c "cd /PTMOverlay && snakemake figures"
  ```
If the workflow fails because of memory overload during the extract_ptms rule, specify the size of your docker container in the snakemake command:
  ```
  docker run -v ./:/PTMOverlay ptm-overlay /bin/bash -c "cd /PTMOverlay && snakemake figures --resources mem_gb=<container-size>"
  ```
Where `<container-size>` is the size of the docker container in GB. This can be as low as 7 but that will significantly slow things down.\
If you do not want to use docker, use the following commands:
   ```
   cd PTMOverlay
   pip install matplotlib pandas pyteomics biopython lxml snakemake ncbi-taxonomist svgutils six
   brew install pdf2svg
   snakemake figures
   ```

The figures on our poster were as follows:
1. PTMOverlay flowchart (made by hand)
2. Enolase PTMs conservation on select intervals (generated thru code and edited in Adobe Illustrator)
3. Conserved PTM sites on enolase structure (made by Youngki You at PNNL)
4. Significant PTM sites in glycolysis pathway (generated thru code and edited in Adobe Illustrator)

The pre-Adobe Illustrator components of figures 2 and 4 will be in the `figures` folder after running.
