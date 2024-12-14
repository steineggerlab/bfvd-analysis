import os
import json
import argparse
import pandas as pd
from tqdm.contrib.concurrent import process_map

from uniprot_api import submit_id_mapping, check_id_mapping_results_ready, get_id_mapping_results_link, get_id_mapping_results_search

def get_uniprot(acc_list, from_db = "UniProtKB_AC-ID", to_db = "UniProtKB"):
    job_id = submit_id_mapping(acc_list, from_db, to_db)
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_search(link)
    return results

def write_json(pools):
    acc, data, dir = pools
    path = os.path.join(dir, acc + ".json")
    with open(path, "wt") as f:
        s = json.dumps(data, indent=2)
        f.write(s)

argparser = argparse.ArgumentParser(description="Generate summary files for 3D Beacons")
argparser.add_argument("--metadata", "-m", help="Path to the metadata file", 
                       default="/home/seamustard52/bfvd-analysis/3d-beacons-bfvd/metadata/bfvd_logan-entry_acc_start_end_len_plddt_taxid_organism_src.tsv")
argparser.add_argument("--uniprot", "-u", help="Path to the UniProt metadata file", 
                       default="/home/seamustard52/bfvd-analysis/3d-beacons-bfvd/metadata/uniprot-acc_length_taxid_organism_src_id_description_gene.tsv")
argparser.add_argument("--save_dir", "-o", help="Path to the directory to save the summary files",
                          default="/home/seamustard52/bfvd-analysis/3d-beacons-bfvd/summary")
argparser.add_argument("--cpus", "-c", help="Number of CPUs to use", type=int, default=os.cpu_count())
argparser.add_argument("--test", "-t", help="Test mode", action="store_true")

args = argparser.parse_args()

bfvd_df = pd.read_csv(args.metadata, sep="\t", header = None, index_col=False,
                      names = ["model", "acc", "start", "end", "length", "pLDDT", "taxid", "organism", "src"],
                      dtype={"model": str, "acc": str, "start": int, "end": int, "length": int, "pLDDT": float, "taxid": str, "organism": str, "src":str}
                      )

uniprot_data = pd.read_csv(args.uniprot, sep="\t", header = None, index_col=False,
                            names = ["acc", "length", "taxid", "organism", "src", "id", "description", "gene"],
                            dtype={"acc": str, "length": int, "taxid": str, "organism": str, "src": str, "id": str, "description": str, "gene": str}
                            ).set_index("acc").T.to_dict()

summary_data = {}

for (acc, data) in uniprot_data.items():
    local_data = {}
    local_data["uniprot_entry"] = {}

    if data["src"] == "UNIPARC":
        local_data["uniprot_entry"] = {"ac": acc, "sequence_length": data["length"]}
    else:
        local_data["uniprot_entry"] = {"ac": acc, "id": data["id"], "sequence_length": data["length"]}
    
    local_data["structures"] = []
    summary_data[acc] = local_data

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

if not args.test:
    pools= [(acc, data, args.save_dir) for acc, data in summary_data.items()]
    process_map(write_json, pools, max_workers=args.cpus)
else:
    for acc, data in summary_data.items():
        write_json((acc, data, args.save_dir))
        break