import os
import pandas as pd
import xml.etree.ElementTree as ET
import requests
from multiprocessing import Pool, Manager
from tqdm import tqdm

def retrieve_uniprot(acc):
    UNIPROT_XML_URL = "https://www.uniprot.org/uniprot"
    ac_id = None
    description = None
    gene = None
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

        try:
            gene = xml_root.find(
                f"./{namespace}entry/{namespace}gene/{namespace}name"
            ).text
        except:
            print(f"LOG: Failed to retrieve gene name for {acc}")
            gene = ""

        print(f"LOG: Retrieved UniProt ID {ac_id} for {acc}")
    except:
        print(f"ERROR: Failed to retrieve UniProt ID for {acc}")

    return ac_id, description, gene

def process_entry(args):
    acc, data = args
    if data["src"] == "UNIPARC":
        id = ""
        desc = ""
        gene = ""
    else:
        id, desc, gene = retrieve_uniprot(acc)
    return acc, id, desc, gene

input = "/home/seamustard52/bfvd-analysis/3d-beacons-bfvd/metadata/bfvd_logan-entry_acc_start_end_len_plddt_taxid_organism_src.tsv"
output = "/home/seamustard52/bfvd-analysis/3d-beacons-bfvd/metadata/uniprot-acc_length_taxid_organism_src_id_description_gene3.tsv"

bfvd_df = pd.read_csv(input, sep="\t", header = None, index_col=False,
                      names = ["model", "acc", "start", "end", "length", "pLDDT", "taxid", "organism", "src"],
                      dtype={"model": str, "acc": str, "start": int, "end": int, "length": int, "pLDDT": float, "taxid": str, "organism": str, "src":str}
                      )
acc_data = bfvd_df[["acc", "length", "taxid", "organism", "src"]].drop_duplicates().set_index("acc").T.to_dict()

# TEST: Sample 10 entries
# acc_data = {k: acc_data[k] for k in list(acc_data.keys())[:10]}

# uniprot_data = {}
args = [(acc, data) for acc, data in acc_data.items()]

with Pool(processes= os.cpu_count()) as pool:
    # results = pool.map(process_entry, args)
    for acc, id, desc, gene in tqdm(pool.imap_unordered(process_entry, args), total=len(args)):
        # uniprot_data[acc] = {"id": id, "description": desc, "gene": gene}
        acc_data[acc]["id"] = id
        acc_data[acc]["description"] = desc
        acc_data[acc]["gene"] = gene
        print(f"LOG: Processed {acc}")

# for acc, data in uniprot_data.items():
#     acc_data[acc]["id"] = data["id"]
#     acc_data[acc]["description"] = data["description"]
#     acc_data[acc]["gene"] = data["gene"]

# for acc, id, desc, gene in results:
#     acc_data[acc]["id"] = id
#     acc_data[acc]["description"] = desc
#     acc_data[acc]["gene"] = gene
#     print(f"LOG: Saved {acc}")

acc_data_df = pd.DataFrame.from_dict(acc_data, orient='index')
acc_data_df.to_csv(output, sep="\t", header=False)
