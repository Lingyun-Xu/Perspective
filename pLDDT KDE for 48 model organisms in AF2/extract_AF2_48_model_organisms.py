#/home/soumajit/sda/LYX/EXTRACT_AF2_Info/extract_AF2_model_organisms/extract_AF2_48_model_organisms.py
#564,446 entries were extracted from 48 model organisms. Same as summarized in https://academic.oup.com/nar/article/52/D1/D368/7337620

import tarfile
import gzip
import io
import pandas as pd

# Paths to local .tar files
tar_file_paths =[
'UP000000429_85962_HELPY_v4.tar',   'UP000001940_6239_CAEEL_v4.tar',      'UP000008524_185431_TRYB2_v4.tar',
'UP000000437_7955_DANRE_v4.tar',    'UP000002059_502779_PARBA_v4.tar',    'UP000008816_93061_STAA8_v4.tar',
'UP000000535_242231_NEIG1_v4.tar',  'UP000002195_44689_DICDI_v4.tar',    'UP000008827_3847_SOYBN_v4.tar',
'UP000000559_237561_CANAL_v4.tar',  'UP000002296_353153_TRYCC_v4.tar',    'UP000008854_6183_SCHMA_v4.tar',
'UP000000579_71421_HAEIN_v4.tar',   'UP000002311_559292_YEAST_v4.tar' ,   'UP000018087_1391915_SPOS1_v4.tar',
'UP000000586_171101_STRR6_v4.tar',  'UP000002438_208964_PSEAE_v4.tar',  'UP000020681_1299332_MYCUL_v4.tar',
'UP000000589_10090_MOUSE_v4.tar',   'UP000002485_284812_SCHPO_v4.tar', 'UP000024404_6282_ONCVO_v4.tar',
'UP000000625_83333_ECOLI_v4.tar',   'UP000002494_10116_RAT_v4.tar',       'UP000030665_36087_TRITR_v4.tar',
'UP000000799_192222_CAMJE_v4.tar',  'UP000002716_300267_SHIDS_v4.tar',    'UP000035681_6248_STRER_v4.tar',
'UP000000803_7227_DROME_v4.tar',    'UP000005640_9606_HUMAN_v4.tar',      'UP000053029_1442368_9EURO2_v4.tar',
'UP000000805_243232_METJA_v4.tar',  'UP000006304_1133849_9NOCA1_v4.tar',  'UP000059680_39947_ORYSJ_v4.tar',
'UP000000806_272631_MYCLE_v4.tar',  'UP000006548_3702_ARATH_v4.tar',      'UP000078237_100816_9PEZI1_v4.tar',
'UP000001014_99287_SALTY_v4.tar',   'UP000006672_6279_BRUMA_v4.tar',      'UP000094526_86049_9EURO1_v4.tar',
'UP000001450_36329_PLAF7_v4.tar',   'UP000007305_4577_MAIZE_v4.tar',      'UP000270924_6293_WUCBA_v4.tar',
'UP000001584_83332_MYCTU_v4.tar',   'UP000007841_1125630_KLEPH_v4.tar',   'UP000274756_318479_DRAME_v4.tar',
'UP000001631_447093_AJECG_v4.tar',  'UP000008153_5671_LEIIN_v4.tar',      'UP000325664_1352_ENTFC_v4.tar']

# Function to parse the CIF file and extract the desired information
def extract_cif_info(cif_content):
    extracted_info = {}
    lines = cif_content.split('\n')
    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if '_struct_ref_seq.pdbx_PDB_id_code' in line:
            extracted_info['pdbx_PDB_id_code'] = line.split()[-1]
        elif '_ma_target_ref_db_details.db_accession' in line:
            extracted_info['db_accession'] = line.split()[-1]
        elif '_ma_target_ref_db_details.db_code' in line:
            extracted_info['db_code'] = line.split()[-1]
        elif '_ma_target_ref_db_details.gene_name' in line:
            extracted_info['gene_name'] = line.split()[-1]
        elif '_ma_target_ref_db_details.ncbi_taxonomy_id' in line:
            extracted_info['ncbi_taxonomy_id'] = line.split()[-1]
        elif '_ma_target_ref_db_details.seq_db_sequence_checksum' in line:
            extracted_info['EMBI_EBI_sequence_checksum'] = line.split()[-1]
        elif '_ma_qa_metric_global.metric_value' in line:
            extracted_info['plddt'] = line.split()[-1]
        elif '_entity_poly.pdbx_seq_one_letter_code' in line and '_can' not in line:
            # Check if the sequence starts on the same line or on the next line
            parts = line.split()
            if len(parts) > 1 and parts[-1] != '_entity_poly.pdbx_seq_one_letter_code':
                # Single-line sequence
                sequence = ''.join(parts[1:])
                extracted_info['AAseq_one_letter_code'] = sequence
            elif parts[-1] == '_entity_poly.pdbx_seq_one_letter_code':
                # Multiline sequence
                sequence = ''
                i += 1
                while lines[i].strip() != ';':
                    sequence += lines[i].strip()
                    i += 1
                extracted_info['AAseq_one_letter_code'] = sequence
        elif '_ma_target_ref_db_details.organism_scientific' in line:
            organism_name_parts = line.split()[1:]
            organism_name = ' '.join(organism_name_parts).strip('"')
            extracted_info['organism_scientific'] = organism_name
        elif '_ma_template_ref_db_details.template_id' in line:
            i += 1
            template_ids = []
            while i < len(lines) and lines[i].strip() and not lines[i].startswith('#'):
                parts = lines[i].split()
                template_id = parts[0]
                template_ids.append(template_id)
                i += 1
            extracted_info['template_ids'] = template_ids

        i += 1
    return extracted_info

total_info = []

for tar_file_path in tar_file_paths:
    try:
        with tarfile.open(tar_file_path, 'r') as tar:
            for member in tar.getmembers():
                if member.name.endswith('.cif.gz'):
                    f = tar.extractfile(member)
                    if f is not None:
                        with gzip.open(io.BytesIO(f.read()), 'rt') as gz:
                            cif_content = gz.read()
                            info = extract_cif_info(cif_content)
                            total_info.append(info)
        print(f"{tar_file_path} finished processing")
    except Exception as e:
        print(f"Error processing {tar_file_path}: {e}")

processed_info = []

for element in total_info:
    processed_element = {}
    for key, value in element.items():
        if isinstance(value, list):
            processed_element[key] = ', '.join(map(str, value))
        else:
            processed_element[key] = value
    processed_info.append(processed_element)

# Create DataFrame from the processed list of dictionaries
df = pd.DataFrame(processed_info)
df.to_parquet('AF2_48_model_oganisms.parquet')
