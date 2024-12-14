import os
import pandas as pd
import xml.etree.ElementTree as ET
import requests
from multiprocessing import Pool
from tqdm import tqdm
import argparse

TIMEOUT = 60
def get_fields(xml_root):
    namespace = "{http://uniprot.org/uniprot}"
    ac_id = None
    description = None
    gene = None
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
        gene = ""

    return ac_id, description, gene

def retrieve_uniprot(acc):
    UNIPROT_XML_URL = "https://www.uniprot.org/uniprot"
    try:
        response = requests.get(f"{UNIPROT_XML_URL}/{acc}.xml", timeout=TIMEOUT)
        xml_root = ET.fromstring(response.content)
        ac_id, description, gene = get_fields(xml_root)
        print(f"LOG: Retrieved UniProt Data for {acc}")

    except:
        try: # One more try
            response = requests.get(f"{UNIPROT_XML_URL}/{acc}.xml", timeout=TIMEOUT)
            xml_root = ET.fromstring(response.content)
            ac_id, description, gene = get_fields(xml_root)
            print(f"LOG: Retrieved UniProt Data for {acc}")
        except:
            ac_id, description, gene = None, None, None
            print(f"ERROR: Failed to retrieve data for {acc}")

    return ac_id, description, gene

def process_entry(pools):
    acc, data = pools
    if data["src"] == "UNIPARC":
        id = ""
        desc = ""
        gene = ""
    else:
        id, desc, gene = retrieve_uniprot(acc)
    return acc, id, desc, gene

if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--input", "-i", type=str, default="/home/seamustard52/bfvd-analysis/3d-beacons-bfvd/metadata/bfvd_logan-entry_acc_start_end_len_plddt_taxid_organism_src.tsv")
    argparser.add_argument("--output", "-o", type=str, required=True)
    argparser.add_argument("--sample", "-s", type=int, default=None)
    args = argparser.parse_args()
    # output = "/home/seamustard52/bfvd-analysis/3d-beacons-bfvd/metadata/uniprot_sample-acc_length_taxid_organism_src_id_description_gene.tsv"

    bfvd_df = pd.read_csv(args.input, sep="\t", header = None, index_col=False,
                        names = ["model", "acc", "start", "end", "length", "pLDDT", "taxid", "organism", "src"],
                        dtype={"model": str, "acc": str, "start": int, "end": int, "length": int, "pLDDT": float, "taxid": str, "organism": str, "src":str},
                        ).fillna("")
    acc_data = bfvd_df[["acc", "length", "taxid", "organism", "src"]].drop_duplicates().set_index("acc").T.to_dict()

    if args.sample:
        acc_data = {k: acc_data[k] for k in list(acc_data.keys())[:args.sample]}

    pools = [(acc, data) for acc, data in acc_data.items()]

    with Pool(processes= os.cpu_count()) as pool:
        for (acc, id, desc, gene) in tqdm(pool.imap_unordered(process_entry, pools), total=len(pools)):
            acc_data[acc]["id"] = id
            acc_data[acc]["description"] = desc
            acc_data[acc]["gene"] = gene

    acc_data_df = pd.DataFrame.from_dict(acc_data, orient='index')
    acc_data_df.to_csv(args.output, sep="\t", header=False)
