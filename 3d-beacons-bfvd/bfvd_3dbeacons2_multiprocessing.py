import os
import json
import requests
import xml.etree.ElementTree as ET
import pandas as pd
from multiprocessing import Pool, Manager

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

def retrieve_uniprot(acc):
    UNIPROT_XML_URL = "https://www.uniprot.org/uniprot"
    ac_id = None
    desc = None
    try:
        xml_root = ET.fromstring(requests.get(f"{UNIPROT_XML_URL}/{acc}.xml").content)
        namespace = "{http://uniprot.org/uniprot}"
        try:
            ac_id = xml_root.find(f"./{namespace}entry/{namespace}name").text
        except:
            ac_id = xml_root.find(f"./{namespace}entry/{namespace}accession").text

        try:
            desc = xml_root.find(
                f"./{namespace}entry/{namespace}protein/{namespace}recommendedName"
                f"/{namespace}fullName"
            ).text
        except:
            desc = xml_root.find(
                f"./{namespace}entry/{namespace}protein/{namespace}submittedName"
                f"/{namespace}fullName"
            ).text        
        print(f"LOG: Retrieved UniProt ID {ac_id} for {acc}")
    except:
        try:
            print(f"LOG: Trying again for {acc}")
            ac_id, desc = retrieve_uniprot(acc)
        except:
            print(f"Error in parsing UniProt XML for {acc}!")
            ac_id = ""
            desc = ""

    return ac_id, desc

def process_fasta_entry(args):
    """Process a single FASTA entry."""
    acc, seq, uniparc = args
    acc_data = {}
    length = len(seq)

    if acc in uniparc:
        acc_data["uniprot_entry"] = {"ac": acc, "sequence_length": length}
        desc = ""
    else:
        id, desc = retrieve_uniprot(acc)
        acc_data["uniprot_entry"] = {"ac": acc, "id": id, "sequence_length": length}

    acc_data["structures"] = []
    print(f"LOG: Processed {acc}")
    return acc, acc_data, desc

# Main code
metadata = "/home/seamustard52/bfvd-analysis/metadata/bfvd_logan-model_length_unk_nmsa_plddt_ptm_logan.tsv"
pathfile = "/home/seamustard52/bfvd-analysis/metadata/bfvd_logan-entry_filepath_rep.tsv"
fasta_origin = "/home/seamustard52/bfvd-analysis/uniref30_2302_db_virusdb_rep.fasta" 
uniparc_list = "/home/seamustard52/bfvd-analysis/3d-beacons-bfvd/uniparc.list"
save_dir = "/home/seamustard52/bfvd-analysis/3d-beacons-bfvd/summary"

bfvd_df = pd.read_csv(metadata, sep="\t", header=None, index_col=False,
                      names=["model", "length", "unkcnt", "nmsa", "pLDDT", "ptm", "logan"])
bfvd_path = pd.read_csv(pathfile, sep="\t", header=None, index_col=False, 
                        usecols=[0, 1], names=["model", "path"])
bfvd_path = bfvd_path.set_index('model').T.to_dict('records')[0]

uniparc = set()
with open(uniparc_list) as f:
    data = f.readlines()
    for line in data:
        uniparc.add(line.strip())

fasta_dict = read_fasta(fasta_origin)
# Sample 10 entries in the FASTA dictionary
# fasta_test = {k: fasta_dict[k] for k in list(fasta_dict.keys())[:10]}

summary_dict = {}
acc_to_desc = {}

# Use Manager to share `uniparc` across processes
with Manager() as manager:
    uniparc_shared = manager.list(uniparc)  # Share the uniparc set across processes

    # Prepare arguments for multiprocessing
    fasta_args = [(acc, seq, uniparc_shared) for acc, seq in fasta_dict.items()]
    # fasta_args = [(acc, seq, uniparc_shared) for acc, seq in fasta_test.items()]

    # Use multiprocessing to process the FASTA entries
    with Pool(processes=os.cpu_count()) as pool:
        results = pool.map(process_fasta_entry, fasta_args)

    # Collect results into the shared dictionary
    for acc, acc_data, desc in results:
        summary_dict[acc] = acc_data
        acc_to_desc[acc] = desc

# Process structures and save to JSON
for index, row in bfvd_df.iterrows():
    model_id = row["model"]
    mapping_acc = model_id.split("_")[0]

    if mapping_acc not in summary_dict:
        continue

    model_url = f"https://bfvd.steineggerlab.workers.dev/cif/{model_id}.cif"
    model_page_url = f"https://bfvd.foldseek.com/cluster/{mapping_acc}"

    if model_id == mapping_acc:
        start = 1
        end = row["length"]
        coverage = 1.0
    else:
        split = int(model_id.split("_")[1])
        full = summary_dict[mapping_acc]["uniprot_entry"]["sequence_length"]
        
        if split == 1:
            start = 1
            end = row["length"]
            coverage = round(row["length"] / full, 2)
        else:
            offset = summary_dict[mapping_acc]["structures"][-1]["summary"]["uniprot_end"]
            start = offset + 1
            end = offset + row["length"]
            coverage = round((end - start + 1) / full, 2)

    desc = acc_to_desc.get(mapping_acc, "")

    summary_dict[mapping_acc]["structures"].append(
        {"summary": {
            "model_identifier": model_id,
            "model_url": model_url,
            "model_page_url": model_page_url,
            "model_format": "MMCIF",
            "model_category": "AB-INITIO",
            "provider": "BFVD",
            "created": "2024-11-01",
            "sequence_identity": 1.0,
            "coverage": coverage,
            "uniprot_start": start,
            "uniprot_end": end,
            "confidence_type": "pLDDT",
            "confidence_avg_local_score": row["pLDDT"],
            "entities": [{"entity_type": "POLYMER",
                          "entity_poly_type": "POLYPEPTIDE(L)",
                          "identifier": mapping_acc,
                          "identifier_category": "UNIPARC" if (mapping_acc in uniparc) else "UNIPROT",
                          "description": desc,
                          "chain_ids": ["A"]
                          }]
            }
        }
    )
    print(f"LOG: Added structure {model_id} to {mapping_acc}")

for acc, data in summary_dict.items():
    p_out = os.path.join(save_dir, acc + ".json")
    with open(p_out, "wt") as f:
        s = json.dumps(data, indent=2)
        f.write(s)
    print(f"LOG: Written summary for {acc} to {p_out}")
