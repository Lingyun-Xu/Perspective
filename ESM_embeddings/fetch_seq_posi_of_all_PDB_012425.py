#!/usr/bin/env python3
import os
import time
import re
import pandas as pd
import requests
import subprocess
import json
import numpy as np
def fetch_protein_pdb_ids(batch_size=1000):
    """
    Fetch all PDB IDs for entries where Polymer Entity Type is Protein.

    Args:
        batch_size (int): Number of rows to fetch per request.

    Returns:
        list: List of PDB IDs.
    """
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    all_pdb_ids = []
    start = 0

    while True:
        # Construct the query payload
        query = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "entity_poly.rcsb_entity_polymer_type",
                    "operator": "exact_match",
                    "value": "Protein"
                }
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {
                    "start": start,
                    "rows": batch_size
                },
                "results_content_type": ["experimental"]
            }
        }

        # Send the POST request
        response = requests.post(url, json=query)
        response.raise_for_status()
        data = response.json()

        # Extract PDB IDs from the results
        pdb_ids = [result['identifier'] for result in data.get('result_set', [])]
        all_pdb_ids.extend(pdb_ids)

        # Break if no more results
        if len(pdb_ids) < batch_size:
            break

        # Update the start for the next batch
        start += batch_size

    return all_pdb_ids


def fetch_pdb_data(query, variables=None, retries=3):
    url = "https://data.rcsb.org/graphql"
    headers = {'Content-Type': 'application/json'}
    payload = {'query': query}
    if variables:
        payload['variables'] = variables

    for attempt in range(retries):
        try:
            response = requests.post(url, json=payload, headers=headers, timeout=10)
            response.raise_for_status()
            data = response.json()
            if 'data' in data:
                return data['data']
            else:
                print(f"Unexpected response format: {data}")
                return None
        except requests.exceptions.RequestException as e:
            print(f"Request failed on attempt {attempt + 1} for query: {query} | Error: {e}")
            if attempt < retries - 1:
                time.sleep(2 ** attempt)  # Exponential backoff
            else:
                return None

def process_pdb_id(pdb_id):
    entity_query = f"""
    {{
      entry(entry_id: "{pdb_id}") {{
        polymer_entities {{
          rcsb_id
          entity_poly {{
            pdbx_seq_one_letter_code
            rcsb_entity_polymer_type
          }}
          rcsb_polymer_entity_container_identifiers {{
            asym_ids
          }}
        }}
      }}
    }}
    """
    entry_data = fetch_pdb_data(entity_query)
    if entry_data is None:
        print(f"fetch_pdb_data returned None for PDB ID: {pdb_id}")
        return []
    elif not entry_data or 'entry' not in entry_data or not isinstance(entry_data['entry'], dict):
        print(f"Invalid response structure: {entry_data}")
        return []

    elif 'polymer_entities' not in entry_data['entry']:
        print(f"Error: Missing or invalid polymer entities for PDB ID: {pdb_id}")
        print(f"Debug: Raw entry_data -> {entry_data}")
        return []

    instance_ids = []
    data_rows = []

    for entity in entry_data['entry']['polymer_entities']:
        if entity['entity_poly']['rcsb_entity_polymer_type'] != "Protein":
            continue  # Skip non-protein sequences
        for asym_id in entity['rcsb_polymer_entity_container_identifiers']['asym_ids']:
            instance_id = f"{pdb_id}.{asym_id}"
            instance_ids.append(instance_id)
            data_rows.append({
                "pdb_id": pdb_id,
                "entity_id": entity['rcsb_id'],
                "sequence": entity['entity_poly']['pdbx_seq_one_letter_code'],
                "instance_id": instance_id,
                "unobserved_residue_xyz": []
            })

    if instance_ids:
        features_query = """
        query getFeatures($instance_ids: [String!]!) {
          polymer_entity_instances(instance_ids: $instance_ids) {
            rcsb_id
            rcsb_polymer_instance_feature {
              type
              feature_positions {
                beg_seq_id
                end_seq_id
              }
            }
          }
        }
        """
        features_data = fetch_pdb_data(features_query, variables={"instance_ids": instance_ids})
        if features_data and 'polymer_entity_instances' in features_data:
            for instance in features_data['polymer_entity_instances']:
                if 'rcsb_polymer_instance_feature' in instance and instance['rcsb_polymer_instance_feature'] is not None:
                    unobserved_features = [
                        feature['feature_positions'] for feature in instance['rcsb_polymer_instance_feature']
                        if feature['type'] == 'UNOBSERVED_RESIDUE_XYZ'
                    ]
                    for row in data_rows:
                        if row['instance_id'] == instance['rcsb_id']:
                            row['unobserved_residue_xyz'].extend(unobserved_features)
        else:
            print(f"No feature data available for instances of PDB ID: {pdb_id}")

    return data_rows

# Batch processing
def batch_process_pdb_ids(pdb_ids, batch_size=1000):
    all_data = []
    total_batches = len(pdb_ids) // batch_size + (1 if len(pdb_ids) % batch_size > 0 else 0)

    for i in range(total_batches):
        batch_ids = pdb_ids[i * batch_size: (i + 1) * batch_size]
        print(f"Processing batch {i + 1} of {total_batches}")
        for pdb_id in batch_ids:
            print(pdb_id)
            data_rows=process_pdb_id(pdb_id)
            if not data_rows:
                continue
            all_data.extend(data_rows)
        time.sleep(1)  # Throttle requests to avoid hitting rate limits

    return pd.DataFrame(all_data)

def main():
    # Step 1: Fetch PDB IDs
    if not os.path.exists("protein_pdb_ids_012425.json"):
        pdb_ids = fetch_protein_pdb_ids()
        with open("protein_pdb_ids_012425.json", "w") as json_file:
            json.dump(pdb_ids, json_file, indent=4)
    else:
        with open("protein_pdb_ids_012425.json", "r") as json_file:
            pdb_ids = json.load(json_file)

    # Step 2: Process PDB IDs
    pdb_ids = [id for id in pdb_ids if id is not None]

    if not pdb_ids:
        print("No valid plant PDB IDs found in the input list.")
    else:
        pdb_df_012425 = batch_process_pdb_ids(pdb_ids)
        pdb_df_012425.to_parquet('fetch_API_for_PDB20250124_seq_and_posi.parquet')

if __name__ == "__main__":
    main()
