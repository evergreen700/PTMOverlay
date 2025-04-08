FROM python:3.12-bookworm

RUN apt-get update -y
RUN apt-get install -y git pdf2svg
RUN pip install matplotlib pandas pyteomics biopython lxml bs4 snakemake ncbi-taxonomist svgutils six
