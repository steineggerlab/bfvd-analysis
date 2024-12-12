import os
import json
import pandas as pd

def read_fasta(file):
    fasta = {}
    with open(file) as f:
        data = f.readlines()
        for line in data:
            line = line.strip()
            if line.startswith(">"):
                header = '_'.join(line.split()[0].split("_")[1:])
                fasta[header] = ""
            else:
                fasta[header] += line
        
    return fasta

def get_index(splitted_seq, origin_seq):
    start_idx = 0
    end_idx = 0
    for idx, aa in enumerate(origin_seq):
        if aa == splitted_seq[0]:
            start_idx = idx
            if origin_seq[start_idx:start_idx+len(splitted_seq)] == splitted_seq:
                end_idx = start_idx + len(splitted_seq) - 1
                break
    return start_idx, end_idx

### Default values
model_category = "AB-INITIO"
model_type = "single"
# oligomeric_state = "MONOMER"
# model_format = "PDB"
# provider = "BFVD"
confidence_type = "pLDDT"
sequence_identity = 1.0
created = "2024-11-01" 
experimental_method = "THEORETICAL MODEL"

metadata = "/home/seamustard52/bfvd-analysis/metadata/bfvd_logan-model_length_unk_nmsa_plddt_ptm_logan.tsv"
pathfile = "/home/seamustard52/bfvd-analysis/metadata/bfvd_logan-entry_filepath_rep.tsv"
fasta_source = "/home/seamustard52/bfvd-analysis/uniref30_2302_db_virusdb_rep_processed.fasta"
fasta_origin = "/home/seamustard52/bfvd-analysis/uniref30_2302_db_virusdb_rep.fasta" 
uniparc_list = "/home/seamustard52/bfvd-analysis/3d-beacons-client/uniparc.list"
save_dir = "/home/seamustard52/bfvd-analysis/3d-beacons-client/data/metadata"
# save_dir = "/home/seamustard52/bfvd-analysis/3d-beacons-client/data/metadata_uniparc"

bfvd_df = pd.read_csv(metadata, sep="\t", header = None, index_col=False,
                      names = ["model", "length", "unkcnt", "nmsa", "pLDDT", "ptm", "logan"])
bfvd_path = pd.read_csv(pathfile, sep="\t", header = None, index_col = False, 
                        usecols=[0,1], names = ["model","path"])
bfvd_path = bfvd_path.set_index('model').T.to_dict('records')[0]

uniparc = set()
with open(uniparc_list) as f:
    data = f.readlines()
    for line in data:
        uniparc.add(line.strip())

fasta_source_dict = read_fasta(fasta_source)
fasta_origin_dict = read_fasta(fasta_origin)

for index, row in bfvd_df.iterrows():
    model_identifier = row["model"]
    mapping_accession = model_identifier.split("_")[0]
    # model_url = ""
    # model_page_url = f"https://bfvd.foldseek.com/cluster/{mapping_accession}"
    confidence_avg_score = row["pLDDT"]

    splitted = True if len(model_identifier.split("_")) == 2 else False
    if not splitted:
        start = 1
        end = row["length"]
        coverage = 1.0
    else :
        start_idx, end_idx = get_index(fasta_source_dict[model_identifier], fasta_origin_dict[mapping_accession])
        start = start_idx + 1
        end = end_idx + 1
        coverage = round((end - start + 1) / len(fasta_origin_dict[mapping_accession]),2)

    if mapping_accession.startswith("UPI") or mapping_accession in uniparc:
        mapping_accession_type = "uniparc"
    else:
        mapping_accession_type = "uniprot"
    #     continue
    metadata = {
        "mappingAccession": mapping_accession,
        "mappingAccessionType": mapping_accession_type,
        "start": start,
        "end": end,
        "modelCategory": model_category,
        "modelType": model_type,
        "experimentalMethod": experimental_method,
        "confidenceType": confidence_type,
        "confidenceAvgLocalScore": confidence_avg_score,
        "createdDate": created,
        "sequenceIdentity": sequence_identity,
        "coverage": coverage,
    }
        
    pdbfile = bfvd_path[model_identifier].split("/")[-1]
    jsonfile = bfvd_path[model_identifier].split("/")[-1].replace(".pdb", ".json")
    with open(f"{save_dir}/{jsonfile}", "w") as f:
        json.dump(metadata, f, indent=4)