import requests
import json

def fetch_pdb_ids_by_organism(organism_name):
    # Define the API endpoint
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    # Maximum rows per request
    max_rows = 1000
    start = 0
    all_pdb_ids = []
    
    while True:
        # Define the query payload with pagination
        payload = {
            "query": {
                "type": "group",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "full_text",
                        "parameters": {
                            "value": organism_name
                        }
                    },
                    {
                        "type": "group",
                        "nodes": [
                            {
                                "type": "terminal",
                                "service": "text",
                                "parameters": {
                                    "attribute": "rcsb_entry_info.structure_determination_methodology",
                                    "value": "experimental",
                                    "operator": "exact_match"
                                }
                            },
                            {
                                "type": "terminal",
                                "service": "text",
                                "parameters": {
                                    "attribute": "rcsb_entity_source_organism.ncbi_scientific_name",
                                    "value": organism_name,
                                    "operator": "exact_match"
                                }
                            },
                            {
                                "type": "terminal",
                                "service": "text",
                                "parameters": {
                                    "attribute": "entity_poly.rcsb_entity_polymer_type",
                                    "value": "Protein",
                                    "operator": "exact_match"
                                }
                            },
                            {
                                "type": "terminal",
                                "service": "text",
                                "parameters": {
                                    "attribute": "rcsb_accession_info.initial_release_date",
                                    "value": {
                                        "from": "2021-01-01",
                                        "to": "2025-01-24",
                                        "include_lower": True,
                                        "include_upper": True
                                    },
                                    "operator": "range"
                                }
                            }
                        ],
                        "logical_operator": "and"
                    }
                ],
                "logical_operator": "and"
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {
                    "start": start,
                    "rows": max_rows
                },
                "results_content_type": ["experimental"],
                "sort": [
                    {
                        "sort_by": "score",
                        "direction": "desc"
                    }
                ],
                "scoring_strategy": "combined"
            }
        }
        
        # Make the POST request
        try:
            response = requests.post(url, json=payload, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            # Extract PDB IDs from the response
            pdb_ids = [result["identifier"] for result in data.get("result_set", [])]
            all_pdb_ids.extend(pdb_ids)
            
            print(f"Retrieved {len(pdb_ids)} PDB IDs in this batch. Total so far: {len(all_pdb_ids)}")
            
            # Break the loop if fewer results are returned than max_rows
            if len(pdb_ids) < max_rows:
                break
            
            # Update the start parameter for the next batch
            start += max_rows
        
        except requests.exceptions.RequestException as e:
            print(f"Failed to fetch PDB IDs for {organism_name}: {e}")
            break

    print(f"Retrieved a total of {len(all_pdb_ids)} PDB IDs for {organism_name}.")
    return all_pdb_ids

# Run the function for multiple organisms
if __name__ == "__main__":
    organisms = ["Arabidopsis thaliana", "Zea mays", "Glycine max", "Oryza sativa"]
    all_pdb_data = {}

    for organism in organisms:
        pdb_ids = fetch_pdb_ids_by_organism(organism)
        all_pdb_data[organism] = pdb_ids

        # Save each organism's results to a separate file
        if pdb_ids:
            filename = f"{organism.replace(' ', '_').replace('.', '')}_pdb_ids_2021_2024_012425.txt"
            with open(filename, "w") as f:
                f.write("\n".join(pdb_ids))
            print(f"PDB IDs for {organism} saved to {filename}.")
        else:
            print(f"No PDB IDs found for {organism}.")

    # Optionally, save all results to a single JSON file
    with open("all_plants_pdb_ids_2021_2024_012425.json", "w") as json_file:
        json.dump(all_pdb_data, json_file, indent=4)
    print("All results saved to all_organisms_pdb_ids_2021_2024.json.")


df = pd.DataFrame([(organism, pdb_id) for organism, pdb_ids in all_pdb_data.items() for pdb_id in pdb_ids],
                  columns=["Organism", "PDB_ID"])

def fetch_pdb_data(query, variables=None):
    url = "https://data.rcsb.org/graphql"
    headers = {'Content-Type': 'application/json'}
    payload = {'query': query}
    if variables:
        payload['variables'] = variables
    try:
        response = requests.post(url, json=payload, headers=headers, timeout=10)
        response.raise_for_status()  # This will raise an HTTPError for bad requests (4XX or 5XX)
        return response.json()['data']
    except requests.exceptions.RequestException as e:
        print(f"Request failed: {e}")
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
#     print(entry_data)
    if not entry_data or 'entry' not in entry_data or not entry_data['entry'].get('polymer_entities'):
        print(f"No polymer entities found or invalid data for PDB ID: {pdb_id}")
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
            all_data.extend(process_pdb_id(pdb_id))
        time.sleep(1)  # Throttle requests to avoid hitting rate limits

    return pd.DataFrame(all_data)


pdb_ids = df['PDB_ID'].values.tolist()
# pdb_ids=['1EG0']
# Filter out None values from the list
pdb_ids = [id for id in pdb_ids if id is not None]

if not pdb_ids:
    print("No valid plant PDB IDs found in the input list.")
else:
    plant_df = batch_process_pdb_ids(pdb_ids)
    plant_df.to_csv('fetch_API_for_2021_to_2024_newly_crystallized_plant_protein_PDB_ID_and_organism_name_seq_and_posi_012405.csv',index=False)
plant_df



