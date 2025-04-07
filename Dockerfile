FROM python:3.12-bookworm

RUN apt-get update -y
RUN apt-get install git -y
RUN pip install matplotlib pandas pyteomics biopython lxml snakemake ncbi-taxonomist
#RUN git clone https://github.com/evergreen700/PTMOverlay && cd PTMOverlay && \ 
#git reset --hard 4fc6458
#ENV HOME=/PTMOverlay
