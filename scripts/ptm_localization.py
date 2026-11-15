
import pandas as pd
import json
import re
import glob
import os
import sys

from Bio import SeqIO
from collections import defaultdict

# ============================================================
# 1. INPUT ARGUMENTS
# ============================================================

localization_file = sys.argv[1] #xlsx file with PTM localization information
infile = sys.argv[2] #fasta files with full-length protein sequences
threshold = float(sys.argv[3]) #localization threshold
ptm_locals = sys.argv[4] #JSON file with PTM information
output_file = sys.argv[5] #output file to save filtered JSON file
output_dir = sys.argv[6] #directory to save filtering reports

# The fourth argument is the path to one JSON file
json_files = [ptm_locals]
if not os.path.isfile(ptm_locals):
    raise FileNotFoundError(
        f"Input JSON file not found: {ptm_locals}"
    )

# ============================================================
# 2. READ LOCALIZATION SPREADSHEET
# ============================================================

# Read the Excel workbook
all_sheets = pd.read_excel(localization_file, sheet_name=None)

# Specify the three tabs you want to include
sheets_to_include = [ "Acetyl_Carbamyl","Methylation","Phosphorylation"]

# Combine only those three tabs
localization_df = pd.concat([all_sheets[sheet] for sheet in sheets_to_include],ignore_index=True)

required_columns = ["Strain","PTM","Site_Probabilities", "Protein"]
missing_columns = [col for col in required_columns if col not in localization_df.columns]
if missing_columns:
    raise ValueError(
        f"Missing spreadsheet columns: {missing_columns}"
    )

localization_df["Protein"] = (localization_df["Protein"].astype(str).str.strip())
localization_df["Strain"] = (localization_df["Strain"].astype(str).str.strip())

# ============================================================
# 3. LOAD ORIGINAL PROTEIN SEQUENCES
# ============================================================

fasta_files = glob.glob(os.path.join(infile, "**", "*.fasta"),recursive=True)

proteins = {}
for fasta_file in fasta_files:
    for record in SeqIO.parse(fasta_file, "fasta"):
        protein_id = record.id
        sequence = str(record.seq)
        if protein_id in proteins:
            if proteins[protein_id] != sequence:
                raise ValueError(
                    f"Conflicting sequences for {protein_id}"
                )
        proteins[protein_id] = sequence

# ============================================================
# 4. PARSE SITE PROBABILITIES
# ============================================================

def parse_site_probabilities(site_string):
    """
    Extract peptide sequence and individual PTM
    localization probabilities.
    Example:
    T(0.0214)GPS(0.0226)T(0.9561)KADAAKEEVK
    Returns:
        peptide:
            TGPSTKADAAKEEVK

        sites:
            [
                {
                    position: 0,
                    residue: T,
                    probability: 0.0214
                },
                ...
            ]
    Positions are zero-based within the peptide.
    """
    pattern = r"([STYKR])\((0(?:\.\d+)?|1(?:\.0+)?)\)"
    sites = []
    for match in re.finditer(pattern, site_string):
        position = sum(char.isalpha() for char in site_string[:match.start()])
        sites.append({
            "position": position,
            "residue": match.group(1),
            "probability": float(match.group(2))
        })
    peptide = re.sub(
        r"\((?:0(?:\.\d+)?|1(?:\.0+)?)\)",
        "",
        site_string
    )
    return peptide, sites

# ============================================================
# 5. FIND PEPTIDE IN FULL-LENGTH PROTEIN
# ============================================================

def map_peptide(protein_id, peptide):
    protein_sequence = proteins.get(protein_id)
    if protein_sequence is None:
        return []
    positions = []
    start = 0
    while True:
        index = protein_sequence.find(peptide, start)
        if index == -1:
            break
        positions.append(index)
        # Allow overlapping matches
        start = index + 1
    return positions

# ============================================================
# 6. BUILD LOCALIZATION LOOKUP
# ============================================================

localization_lookup = defaultdict(list)
unmapped_results = []
ambiguous_results = []
for _, row in localization_df.iterrows():
    protein_id = row["Protein"]
    strain = row["Strain"]
    ptm_type = row["PTM"]
    site_string = row["Site_Probabilities"]
    if pd.isna(site_string):
        continue
    peptide, sites = parse_site_probabilities(
        str(site_string)
    )
    if not sites:
        continue
    peptide_starts = map_peptide(
        protein_id,
        peptide
    )
    if len(peptide_starts) == 0:
        unmapped_results.append({
            "strain": strain,
            "protein": protein_id,
            "peptide": peptide,
            "reason": "Peptide not found in protein"
        })
        continue
    if len(peptide_starts) > 1:
        ambiguous_results.append({
            "strain": strain,
            "protein": protein_id,
            "peptide": peptide,
            "reason": "Multiple peptide matches"
        })
        continue
    peptide_start = peptide_starts[0]
    for site in sites:
        protein_position = ( peptide_start + site["position"])
        key = (strain,protein_id,ptm_type,protein_position)
        localization_lookup[key].append(site["probability"])

# ============================================================
# STEP 7: FILTER EXISTING JSON FILES USING LOCALIZATION SCORES
# ============================================================

def get_ptm_type_from_filename(json_file):
    filename = os.path.basename(json_file).lower()
    ptm_types = {
        "acetyl": "Acetylation (K)",
        "carbamyl": "Carbamylation (K)",
        "dimethyl": "Dimethylation (K/R)",
        "monomethyl": "Monomethylation (K/R)",
        "phospho": "Phosphorylation (STY)",
        "trimethyl": "Trimethylation (K/R)"
    }
    for name, ptm_type in ptm_types.items():
        if name in filename:
            return ptm_type
    raise ValueError(
        f"Could not identify PTM type from filename: {filename}"
    )

retained_report = []
removed_report = []
unmatched_report = []

for json_file in json_files:
    ptm_type = get_ptm_type_from_filename(json_file)
    with open(json_file, "r") as f:
        data = json.load(f)
    filtered_data = {}
    for protein_key, entries in data.items():
        parts = protein_key.split(", ")
        if len(parts) < 4:
            unmatched_report.append({
                "json_file": json_file,
                "protein_key": protein_key,
                "reason": "Could not parse JSON protein key"
            })
            # Preserve the entire protein entry
            filtered_data[protein_key] = entries.copy()
            continue
        protein_id = parts[0]
        strain = parts[1]
        existing_ptms = entries["mod_site"]
        retained_ptms = []
        for position in existing_ptms:
            key = (
                strain,
                protein_id,
                ptm_type,
                int(position)
            )
            matching_scores = localization_lookup.get(key, [])
            if not matching_scores:
                retained_ptms.append(position)
                unmatched_report.append({
                    "json_file": json_file,
                    "protein": protein_id,
                    "strain": strain,
                    "ptm_type": ptm_type,
                    "position": position,
                    "reason": "No matching localization result"
                })
                continue
            best_probability = max(matching_scores)
            if best_probability >= threshold:
                retained_ptms.append(position)
                retained_report.append({
                    "json_file": json_file,
                    "protein": protein_id,
                    "strain": strain,
                    "ptm_type": ptm_type,
                    "position": position,
                    "probability": best_probability
                })
            else:
                # Only sites in this category are removed
                removed_report.append({
                    "json_file": json_file,
                    "protein": protein_id,
                    "strain": strain,
                    "ptm_type": ptm_type,
                    "position": position,
                    "probability": best_probability
                })
        filtered_data[protein_key] = entries.copy()
        # Keep retained AND unmatched sites
        filtered_data[protein_key]["mod_site"] = retained_ptms
    with open(output_file, "w") as f:
        json.dump(filtered_data, f, indent=4)

# ============================================================
# STEP 8: SAVE FILTERING REPORTS
# ============================================================

os.makedirs(output_dir, exist_ok=True)

pd.DataFrame(retained_report).to_csv(
    os.path.join(output_dir, "retained_ptms.csv"),
    mode="a",
    header=not os.path.exists(os.path.join(output_dir, "retained_ptms.csv")),
    index=False
)

pd.DataFrame(removed_report).to_csv(
    os.path.join(output_dir, "removed_ptms.csv"),
    mode="a",
    header=not os.path.exists(os.path.join(output_dir, "removed_ptms.csv")),
    index=False
)

pd.DataFrame(unmatched_report).to_csv(
    os.path.join(output_dir, "unmatched_ptms.csv"),
    mode="a",
    header=not os.path.exists(os.path.join(output_dir, "unmatched_ptms.csv")),
    index=False
)