## Figure regeneration (Capstone Reproducibility)
This can be reproduced on Mac and Linux machines.\
To generate the figures used in BIO 465: Bioinformatics Capstone, follow these steps: 
1. Install docker from https://docs.docker.com
2. Clone this repository by running this in your terminal: \
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

The figures on our poster were as follows:
1. PTMOverlay flowchart (made by hand)
2. Enolase PTMs conservation on select intervals (generated thru code and edited in Adobe Illustrator)
3. Conserved PTM sites on enolase structure (made by Youngki You at PNNL)
4. Significant PTM sites in glycolysis pathway (generated thru code and edited in Adobe Illustrator)

The pre-Adobe Illustrator components of figures 2 and 4 will be in the `figures` folder after running.
