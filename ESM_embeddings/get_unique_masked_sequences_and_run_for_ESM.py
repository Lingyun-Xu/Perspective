import re
import requests
import pandas as pd
import ast
import numpy as np
import math
import os
import subprocess

def is_dna_rna(sequence):
    return all(char in 'ACGT' for char in sequence) or all(char in 'ACGU' for char in sequence) or all(char in '(DA)(DC)(DT)(DG)(DU)(UNK)(PED)(C49)(5CM)(GTP)' for char in sequence)

df_only_AA = df[~df['sequence'].apply(is_dna_rna)]
df_only_AA_noncanical = df_only_AA[df_only_AA['sequence'].str.contains(r'\(')]
pattern = r"\(([^\)]*)\)"
# Extract all strings inside parentheses from the 'sequence' column
matches = df_only_AA_noncanical["sequence"].str.extractall(pattern)
# 'matches' is a DataFrame with match groups. The actual text is in column 0.
unique_noncanonical = matches[0].unique()
has_parent = {}
no_parent_found = []
request_errors = []
for code in unique_noncanonical:
    # Skip empty or weird strings
    code_stripped = code.strip()
    if not code_stripped:
        no_parent_found.append(code)
        continue
    
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{code_stripped}"
    try:
        resp = requests.get(url, timeout=5)
        if resp.status_code == 200:
            data = resp.json()
            # Fetch the field that typically indicates the parent residue
            parent_id = data.get("chem_comp", {}).get("mon_nstd_parent_comp_id", None)
            if parent_id:
                has_parent[code_stripped] = parent_id
            else:
                no_parent_found.append(code_stripped)
        else:
            # If the code doesn't exist in the RCSB dictionary or is not found
            no_parent_found.append(code_stripped)
    except requests.exceptions.RequestException as e:
        print(f"Error fetching {code_stripped} from RCSB: {e}")
        request_errors.append((code_stripped, str(e)))

# Standard 3-letter -> single-letter dictionary
standard_3to1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    "SEC": "U",  # selenocysteine
    "PYL": "O",  # pyrrolysine
}

no_parent_set = set(no_parent_found)
request_error_set = set(x[0] for x in request_errors)  # extract just the code

def map_code_to_single_letter(code):
    # 1) If this code is recognized in has_parent
    if code in has_parent:
        parent_list = has_parent[code]
        # If there's exactly one parent, try to map it
        if len(parent_list) == 1:
            parent_3letter = parent_list[0]
            # Does it map to a standard single letter?
            if parent_3letter in standard_3to1:
                return standard_3to1[parent_3letter]
            else:
                # parent_3letter not recognized in standard_3to1
                return "X"
        else:
            # multiple possible parents → ambiguous, return 'X'
            return "X"
    # 2) If code is explicitly in no_parent or request_errors, default to 'X'
    elif code in no_parent_set or code in request_error_set:
        return "X"
    else:
        # 3) Anything else not covered → 'X'
        return "X"

# Regex to find parenthesized text
pattern = re.compile(r"\((.*?)\)")

def replace_noncanonical(sequence):
    """
    For each '(CODE)' in the sequence, replace with the
    mapped single-letter or 'X', as determined by map_code_to_single_letter().
    """
    def _replacer(match):
        code = match.group(1).strip()  # e.g. "MSE", "BEB", etc.
        return map_code_to_single_letter(code)

    # Perform substitution across the entire string
    return pattern.sub(_replacer, sequence)

df_only_AA["sequence"] = df_only_AA["sequence"].apply(replace_noncanonical)

#check type of "unobserved_residue_xyz"
row= df_only_AA.iloc[3] # any position
content=row['unobserved_residue_xyz']
inside_content=content[0]
print(content)
print(type(content))
print(inside_content)
print(type(inside_content))

#if content is array
def flatten_unobserved_ndarray(val):
    if isinstance(val, np.ndarray) and val.shape == (1,):
        # Now val[0] is another array, e.g. shape (N,) of dict
        inner = val[0]
        if isinstance(inner, np.ndarray) and inner.dtype == object:
            # Convert to a real Python list of dicts
            return inner.tolist()
    # Fallback for empty or unexpected structures
    return []
def extract_seq_ids_final(row):
    flattened = flatten_unobserved_ndarray(row['unobserved_residue_xyz'])
    if flattened:
        return flattened, len(flattened)
    else:
        return np.nan, np.nan



# if content is str
def extract_seq_ids_final(row):
    raw_val = row["unobserved_residue_xyz"]
    # If it's not even a string, skip
    if not isinstance(raw_val, str) or not raw_val.strip():
        return np.nan, np.nan
    # Safely parse the Python-like string
    try:
        parsed = ast.literal_eval(raw_val)  # e.g. [[{'beg_seq_id':..., 'end_seq_id':...}]]
    except (SyntaxError, ValueError):
        # If parsing fails, treat as no data
        return np.nan, np.nan
    # Now 'parsed' should be a list (outer) that might contain one or more sub-lists
    if isinstance(parsed, list) and len(parsed) > 0:
        # For example, parsed[0] might be [ {'beg_seq_id':..., 'end_seq_id':...}, ...]
        first_elem = parsed[0]
        if isinstance(first_elem, list):
            # That inner list presumably contains the dicts
            return first_elem, len(first_elem)
    # Fallback if structure is missing or not as expected
    return np.nan, np.nan




# if content is a list
def extract_seq_ids_final(row):
    """
    Given a row whose 'unobserved_residue_xyz' column might look like:
      [[{'beg_seq_id': 1, 'end_seq_id': 17}, {'beg_seq_id': 117, 'end_seq_id': 122}]]
    This function returns (list_of_dicts, length_of_list).
    If it's empty or not a list, returns (NaN, NaN).
    """
    val = row['unobserved_residue_xyz']  # e.g. [[{...}, {...}], ...] or []
    
    # Check if val is a non-empty list
    if isinstance(val, list) and len(val) > 0:
        # Typically, val[0] is the actual list of dicts
        first_element = val[0]
        if isinstance(first_element, list) and len(first_element) > 0:
            return first_element, len(first_element)
    
    # Fallback if structure doesn't match or is empty
    return np.nan, np.nan




df_only_AA02 = df_only_AA.copy()

df_only_AA02['positions_to_be_masked'], \
df_only_AA02['number_of_positions_to_be_masked'] = zip(
    *df_only_AA02.apply(extract_seq_ids_final, axis=1)
)


def mask_amino_acid(sequence, positions_dict, offset=0):
    # Convert positions from 1-based to 0-based index
    beg_seq_id = positions_dict['beg_seq_id'] - 1 + offset
    end_seq_id = positions_dict['end_seq_id'] - 1 + offset
    # Ensure the indices are within the bounds of the sequence
    if beg_seq_id < 0 or end_seq_id >= len(sequence) or beg_seq_id > end_seq_id:
        raise ValueError("Invalid indices")
    # Replace the entire range with <mask>
    mask_length = end_seq_id - beg_seq_id + 1
    masked_sequence = sequence[:beg_seq_id] + "<mask>" * mask_length + sequence[end_seq_id + 1:]
    # Calculate the new offset after inserting <mask> tags
    new_offset = offset + (mask_length * 6) - mask_length
    return masked_sequence, new_offset
# Apply the masking function to each row
def apply_masking(row):
    sequence = row['sequence']
    positions_to_be_masked = row['positions_to_be_masked']
    if not isinstance(positions_to_be_masked, list):
        return sequence
    offset = 0
    for position in positions_to_be_masked:
        sequence, offset = mask_amino_acid(sequence, position, offset)
    return sequence

df_only_AA02['sequence_finished_mask'] = df_only_AA02.apply(apply_masking, axis=1)

df_only_AA02 = df_only_AA02.drop_duplicates(subset='sequence_finished_mask')
df_only_AA02.to_csv('xxxxxx.csv',index=False) #fill the blank xxxxxx to get "2021-2024_published_four_plant_species_in_PDB_012425_mask_unique.csv" and "sequence_finished_mask_filtered_PDB_012525_mask_unique.csv" 

"""
let_me_know_what_to_exclude_seq = pd.merge(pdb_df_012525_only_AA02, plant_df_012525_only_AA02_unique_mask, on='sequence', how='inner', suffixes=('_pdb', '_plantpdb'))
sequences_to_exclude = let_me_know_what_to_exclude_seq["sequence"].unique()
pdb_df_excluded = pdb_df_012525_only_AA02[
    ~pdb_df_012525_only_AA02["sequence"].isin(sequences_to_exclude)
]
pdb_df_excluded.to_csv('sequence_finished_mask_filtered_PDB_012525_mask_unique_exclude_allseq_from2021-20250117_plants.csv',index=False)
"""




output_dir="/home/lx5/LYX/redo_PDB/ESM_embedding"
def write_sequences_to_fasta(sequences, sequence_ids, output_file):
    """
    Save sequences to a FASTA file, excluding any sequence containing "(" or ")".
    """
    with open(output_file, 'w') as f:
        for sequence_id, sequence in zip(sequence_ids, sequences):
            if "(" not in sequence and ")" not in sequence:  # Check for parentheses
                f.write(f'>{sequence_id}\n{sequence}\n')
            else:
                print(f'Skipped sequence with ID {sequence_id} due to presence of parentheses.')
csv_files = [
"2021-2024_published_four_plant_species_in_PDB_012425_mask_unique.csv","sequence_finished_mask_filtered_PDB_012525_mask_unique.csv"
,"sequence_finished_mask_filtered_PDB_012525_mask_unique_exclude_allseq_from2021-20250117_plants.csv"]

for csv_file in csv_files:
    df = pd.read_csv(csv_file)

    # Extract sequences and ids
    sequences = df['sequence_finished_mask'].tolist()
    sequence_ids = df['instance_id'].tolist()

    # Construct new FASTA file name
    base_name = os.path.basename(csv_file)
    new_file_name = base_name.replace('.csv', '.fasta')
    new_file_path = os.path.join(output_dir, new_file_name)

    # Write sequences and ids to a new FASTA file
    write_sequences_to_fasta(sequences, sequence_ids, new_file_path)
    print(f'Processed {csv_file} and saved to {new_file_path}')


def remove_redundancy(input_file_path, output_file_path):
    unique_sequences = {}
    label_counts = {}
    current_sequence_name = ''
    current_sequence = ''

    with open(input_file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):  # It's a sequence label
                if current_sequence:
                    # If the sequence is not already added, add it to the dictionary with its label
                    if current_sequence not in unique_sequences:
                        unique_sequences[current_sequence] = current_sequence_name
                    current_sequence = ''  # Reset the current sequence for the next one
                # Handle the sequence name for uniqueness
                label_base = line.split()[0]  # Assuming the unique part of the label is before any spaces
                label_counts[label_base] = label_counts.get(label_base, 0) + 1
                # Modify label if it's a duplicate
                current_sequence_name = f"{label_base}_{label_counts[label_base]}" if label_counts[label_base] > 1 else label_base
            else:
                current_sequence += line.strip()

        # Handle the last sequence
        if current_sequence and current_sequence not in unique_sequences:
            unique_sequences[current_sequence] = current_sequence_name

    # Write the unique sequences to the output file
    with open(output_file_path, 'w') as output_file:
        for sequence, label in unique_sequences.items():
            output_file.write(f'>{label}\n')
            output_file.write(f'{sequence}\n')
directory_path = '/home/lx5/LYX/redo_PDB/ESM_embedding'
fasta_files=["sequence_finished_mask_filtered_PDB_012525_mask_unique_exclude_allseq_from2021-20250117_plants.fasta",
"sequence_finished_mask_filtered_PDB_012525_mask_unique.fasta",
"2021-2024_published_four_plant_species_in_PDB_012425_mask_unique.fasta"]
for fasta_file in fasta_files:
    input_file_path = os.path.join(directory_path, fasta_file)
    output_file_path = os.path.join(directory_path, f'cleaned_{fasta_file}')
    remove_redundancy(input_file_path, output_file_path)
    print(f'Processed {fasta_file} and saved to {output_file_path}')

#extract.py and esm2_t33_650M_UR50D are downloaded from Evolutionary Scale Modeling (https://github.com/facebookresearch/esm); environment set up accordingly.
subprocess.run([
        "python3", "extract.py", "esm2_t33_650M_UR50D",
        f'cleaned_{fasta_file}', "output_pt_files",
        "--repr_layer", "33", "--include", "mean"
    ])

