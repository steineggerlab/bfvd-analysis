import os
import json
import requests
import xml.etree.ElementTree as ET
import pandas as pd


def read_fasta(file):
    fasta = {}
    with open(file) as f:
        data = f.readlines()
        for line in data:
            line = line.strip()
            if line.startswith(">"):
                # Header looks like
                # >UniRef100_A0A2P1GMZ4 Replicase polyprotein 1ab n=1 Tax=Hainan hebius popei torovirus TaxID=2116385 RepID=A0A2P1GMZ4_9NIDO
                header = '_'.join(line.split()[0].split("_")[1:])
                fasta[header] = ""
            else:
                fasta[header] += line
        
    return fasta

def retrieve_uniprot(acc):
    UNIPROT_XML_URL = "https://www.uniprot.org/uniprot"
    ac_id = None
    description = None
    try:
        xml_root = ET.fromstring(requests.get(f"{UNIPROT_XML_URL}/{acc}.xml").content)
        namespace = "{http://uniprot.org/uniprot}"
        try:
            ac_id = xml_root.find(f"./{namespace}entry/{namespace}name").text
        except:
            ac_id = xml_root.find(f"./{namespace}entry/{namespace}accession").text

        try:
            description = xml_root.find(
                f"./{namespace}entry/{namespace}protein/{namespace}recommendedName"
                f"/{namespace}fullName"
            ).text
        except:
            description = xml_root.find(
                f"./{namespace}entry/{namespace}protein/{namespace}submittedName"
                f"/{namespace}fullName"
            ).text        
        print(f"LOG: Retrieved UniProt ID {ac_id} for {acc}")
    except:
        try:
            print(f"LOG: Trying again for {acc}")
            ac_id, description = retrieve_uniprot(acc)
        except:
            print(f"Error in parsing UniProt XML for {acc}!")
            ac_id = ""
            description = ""

    return ac_id, description

metadata = "/home/seamustard52/bfvd-analysis/metadata/bfvd_logan-model_length_unk_nmsa_plddt_ptm_logan.tsv"
pathfile = "/home/seamustard52/bfvd-analysis/metadata/bfvd_logan-entry_filepath_rep.tsv"
fasta_origin = "/home/seamustard52/bfvd-analysis/uniref30_2302_db_virusdb_rep.fasta" 
uniparc_list = "/home/seamustard52/bfvd-analysis/3d-beacons-client/uniparc.list"
save_dir = "/home/seamustard52/bfvd-analysis/3d-beacons-bfvd/summary"

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


fasta_dict = read_fasta(fasta_origin)
acc_to_desc = {}
summary_dict = {}

for (acc, seq) in fasta_dict.items():
    acc_data = {}
    length = len(seq)

    if acc.startswith("UPI"):
        uniparc.add(acc)

    if acc in uniparc:
        acc_data["uniparc_entry"] = {"ac": acc,
                                "sequence_length": length}
    else:

        id, desc = retrieve_uniprot(acc)
        acc_to_desc[acc] = desc
        acc_data["uniprot_entry"] = {"ac": acc,
                                "id": id,
                                "sequence_length": length}
    
    acc_data["structures"] = []
    summary_dict[acc] = acc_data
    print(f"LOG: Processed {acc}")

for index, row in bfvd_df.iterrows():
    model_id = row["model"]
    mapping_acc = model_id.split("_")[0]

    model_url = f"https://bfvd.steineggerlab.workers.dev/cif/{model_id}.cif"
    # model_page_url = f"https://bfvd.foldseek.com/cluster/{mapping_acc}"

    if model_id == mapping_acc:
        start = 1
        end = row["length"]
        coverage = 1.0
    else:
        # split
        split = int(model_id.split("_")[1])
        if mapping_acc in uniparc:
            full = summary_dict[mapping_acc]["uniparc_entry"]["sequence_length"]
        else:
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

    try:
        desc = acc_to_desc[mapping_acc]
    except:
        desc = "" # For uniparc entries

    summary_dict[mapping_acc]["structures"].append(
        {"summary": {
            "model_identifier": model_id,
            "model_url": model_url,
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
                          "identifier_category": "UNIPARC" if mapping_acc in uniparc else "UNIPROT",
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