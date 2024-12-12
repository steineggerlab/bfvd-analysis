import os
import json
import argparse
import pandas as pd

from uniprot_api import submit_id_mapping, check_id_mapping_results_ready, get_id_mapping_results_link, get_id_mapping_results_search
from retrieve_uniprot import retrieve_uniprot

def get_uniprot(acc_list, from_db = "UniProtKB_AC-ID", to_db = "UniProtKB"):
    job_id = submit_id_mapping(acc_list, from_db, to_db)
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_search(link)
    return results

metadata = "/home/seamustard52/bfvd-analysis/3d-beacons-bfvd/metadata/bfvd_logan-entry_acc_start_end_len_plddt_taxid_organism_src.tsv"
uniprot = "/home/seamustard52/bfvd-analysis/3d-beacons-bfvd/metadata/uniprot-acc_length_taxid_organism_src_id_description_gene.tsv"
save_dir = "/home/seamustard52/bfvd-analysis/3d-beacons-bfvd/summary"

bfvd_df = pd.read_csv(metadata, sep="\t", header = None, index_col=False,
                      names = ["model", "acc", "start", "end", "length", "pLDDT", "taxid", "organism", "src"],
                      dtype={"model": str, "acc": str, "start": int, "end": int, "length": int, "pLDDT": float, "taxid": str, "organism": str, "src":str}
                      )

uniprot_data = pd.read_csv(uniprot, sep="\t", header = None, index_col=False,
                            names = ["acc", "length", "taxid", "organism", "src", "id", "description", "gene"],
                            dtype={"acc": str, "length": int, "taxid": str, "organism": str, "src": str, "id": str, "description": str, "gene": str}
                            ).set_index("acc").T.to_dict()
# acc_data = bfvd_df[["acc", "length", "taxid", "organism", "src"]].drop_duplicates().set_index("acc").T.to_dict()

# ### Retrieve data (id, description, gene name) from UniProt
# acc_list = [acc for acc, data in acc_data.items() if data["src"] == "UNIPROT"]
# results = get_uniprot(acc_list)

summary_data = {}
# TODO Remove: Sample 10 entries
# acc_data = {k: acc_data[k] for k in list(acc_data.keys())[:10]}

for (acc, data) in uniprot_data.items():
    local_data = {}
    local_data["uniprot_entry"] = {}

    if data["src"] == "UNIPARC":
        local_data["uniprot_entry"] = {"ac": acc, "sequence_length": data["length"]}
    else:
        local_data["uniprot_entry"] = {"ac": acc, "id": id, "sequence_length": data["length"]}
    
    local_data["structures"] = []
    summary_data[acc] = local_data
    print(f"LOG: Read Uniprot {acc}")

for index, row in bfvd_df.iterrows():
    model_id = row["model"]
    mapping_acc = row["acc"]

    if mapping_acc not in summary_data:
        continue

    model_url = f"https://bfvd.steineggerlab.workers.dev/cif/{model_id}.cif"
    model_page_url = f"https://bfvd.foldseek.com/cluster/{mapping_acc}"

    desc = uniprot_data[mapping_acc]["description"]

    summary_data[mapping_acc]["structures"].append(
        {"summary": {
            "model_identifier": model_id,
            "model_url": model_url,
            "model_page_url": model_page_url,
            "model_format": "MMCIF",
            "model_category": "AB-INITIO",
            "provider": "BFVD",
            "created": "2024-11-01",
            "sequence_identity": 1.0,
            "coverage": round((row["end"] - row["start"] + 1) / row["length"], 2),
            "uniprot_start": row["start"],
            "uniprot_end": row["end"],
            "confidence_type": "pLDDT",
            "confidence_avg_local_score": row["pLDDT"],
            "entities": [{"entity_type": "POLYMER",
                          "entity_poly_type": "POLYPEPTIDE(L)",
                          "identifier": mapping_acc,
                          "identifier_category": row["src"],
                          "description": desc,
                          "chain_ids": ["A"]
                          }]
            }
        }
    )
    print(f"LOG: Added structure {model_id} to {mapping_acc}")

for acc, data in summary_data.items():
    p_out = os.path.join(save_dir, acc + ".json")
    with open(p_out, "wt") as f:
        s = json.dumps(data, indent=2)
        f.write(s)
    print(f"LOG: Written summary for {acc} to {p_out}")
    break # TODO: remove