from bs4 import BeautifulSoup
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
import json
import re

infile = sys.argv[1] #html
outfile = sys.argv[2] #fasta
json_outfile = sys.argv[3] #json

# --- Load HTML --- #
with open(infile, "r", encoding="utf-8") as f:
    soup = BeautifulSoup(f, "html.parser")

table_rows = soup.find_all("table")[0].find_all("tr")

species_dict = {}

for row in table_rows:
    cols = row.find_all("td")
    if len(cols) < 2:
        continue
    species_index = cols[0].text.strip().replace(":", "")
    species_name = cols[1].text.strip()
    match = re.search(r"GCA_\d+\.\d+", species_name)
    if match:
        species_dict[int(species_index)] = match.group(0)

# ---- Extract species names and sequences ---- #
species_data = {}
modifications = {}

table_rows = soup.find_all("table")[1].find_all("tr")[1:]

species_set = set()
for row in table_rows:
    cols = row.find_all("td")
    if len(cols) < 2:
        continue
    species_id = cols[0].text.strip()
    species_set.add(species_id)
species_count = len(species_set)

filename = os.path.basename(os.path.dirname(infile))
split_title = filename.split("_")
ptm_types = split_title

for row in table_rows[:species_count]:
    cols = row.find_all("td")
    if len(cols) < 2:
        continue
    species_id = cols[0].text.strip()
    sequence_tag = cols[1].find("pre")
    if not sequence_tag:
        continue
    sequence = sequence_tag.get_text(strip=True)
    mod_positions = {}
    index = 0
    for span in cols[1].find_all("span"):
        residue = span.text.strip()
        if residue:
            for mod in ptm_types:
                if mod in span.get("class", []):
                    mod_positions[index] = mod
            index += 1
    species_data[species_id] = sequence
    modifications[species_id] = mod_positions

# ---- Convert to Biopython format ---- #
alignment_records = [SeqRecord(Seq(seq), id=name, description="") for name, seq in species_data.items()]

with open(outfile, "w") as fasta_file:
    SeqIO.write(alignment_records, fasta_file, "fasta")

# Save species_data and modifications to a JSON file
species_info = {
    "species_data": species_data,
    "modifications": modifications,
    "GCA_codes": species_dict
}

with open(json_outfile, "w") as json_file:
    json.dump(species_info, json_file, indent=4)

