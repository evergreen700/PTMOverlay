from pathlib import Path
import os
import json

ORTHOLOGS=["K07560","K07566"]
PROTEOMES="proteome"
PRE_ALIGN_FASTAS="preAlign"
RAW_ALIGNMENTS="rawAlign"
MUSCLE='muscle-linux-x86.v5.3'

rule preAlignBenchmark:
  input:
    fastas=expand(RAW_ALIGNMENTS+'/{ko}.faa', ko=ORTHOLOGS)

rule muscle:
  input:
    fasta=PRE_ALIGN_FASTAS+'/{ko}.faa'
  output:
    alignment=RAW_ALIGNMENTS+'/{ko}.faa'
  shell:
    '''
    mkdir -p {RAW_ALIGNMENTS}
    ./{MUSCLE} -align {input.fasta} -output {output.alignment}
    '''

rule group_orthologs:
  output:
    fastas=expand(PRE_ALIGN_FASTAS+'/{ko}.faa', ko=ORTHOLOGS)
  shell:
    '''
    python3 group_orthologs.py {PROTEOMES} {PRE_ALIGN_FASTAS}
    '''    
