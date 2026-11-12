#!/usr/bin/env python3
"""
=== INPUTS ===
data/VMR_MSL41.v1.20260729.xlsx                            - ICTV VMR (rank columns D..R)
data/bfvd2_taxid_rank_scientificname_lineage.tsv           - BFVD v2: id, taxid, rank, scientificname, lineage
data/bfvd3_uniprot_2025_03_1st-acc_taxid_taxname_host.tsv  - BFVD v3: id, taxid, taxname, host
=============

ICTV taxa covered per rank by BFVD v2 and v3.

Outputs (analysis/results/coverage/):
  ictv_rank_totals.tsv           rank, n_ictv_taxa            (ICTV's own distinct taxa per rank)
  lineage_by_taxid_ictv_only.tsv dataset, taxid, rank, rank_taxid, rank_name
                                 (dataset v2 / v3; rank_taxid = rank_name = the ICTV taxon name)
"""
import csv
import os

import openpyxl

from ictv_mapping import IctvMapper, VMR_XLSX

REPO_ROOT = os.environ.get("BFVD_ROOT", os.getcwd())
V2_TSV = os.path.join(REPO_ROOT, "data/bfvd2_taxid_rank_scientificname_lineage.tsv")
V3_TSV = os.path.join(REPO_ROOT, "data/bfvd3_uniprot_2025_03_1st-acc_taxid_taxname_host.tsv")
RESULTS = os.path.join(REPO_ROOT, "analysis/results/coverage")
TOTALS_PATH = os.path.join(RESULTS, "ictv_rank_totals.tsv")
LINEAGE_PATH = os.path.join(RESULTS, "lineage_by_taxid_ictv_only.tsv")

RANKS = ["realm", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
COL_IDX = {"realm": 3, "kingdom": 5, "phylum": 7, "class": 9, "order": 11, "family": 13, "genus": 15, "species": 17}


def load_taxa(path, taxid_col, name_col):
    """{taxid: organism name}, first occurrence of each taxid."""
    taxa = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= name_col or not parts[taxid_col]:
                continue
            taxa.setdefault(int(parts[taxid_col]), parts[name_col])
    return taxa


def load_vmr_ranks():
    """-> ({rank: set of ICTV taxon names}, {ictv_id: {rank: name}})"""
    wb = openpyxl.load_workbook(VMR_XLSX, read_only=True, data_only=True)
    rows = wb["VMR MSL41"].iter_rows(values_only=True)
    header = next(rows)
    for rank, idx in COL_IDX.items():
        assert header[idx] == rank.capitalize(), (rank, header[idx])
    assert header[18] == "ICTV_ID"
    totals = {rank: set() for rank in RANKS}
    ictv_ranks = {}
    for row in rows:
        for rank, idx in COL_IDX.items():
            if row[idx]:
                totals[rank].add(row[idx])
        if row[18] not in ictv_ranks:  # isolates of one species repeat the same ranks
            ictv_ranks[row[18]] = {rank: row[idx] for rank, idx in COL_IDX.items() if row[idx]}
    return totals, ictv_ranks


def main():
    os.makedirs(RESULTS, exist_ok=True)
    totals, ictv_ranks = load_vmr_ranks()
    with open(TOTALS_PATH, "w", newline="") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow(["rank", "n_ictv_taxa"])
        for rank in RANKS:
            writer.writerow([rank, len(totals[rank])])

    mapper = IctvMapper()
    mapped = {"v2": mapper.map_taxa(load_taxa(V2_TSV, 1, 3), "BFVD v2"),
              "v3": mapper.map_taxa(load_taxa(V3_TSV, 1, 2), "BFVD v3")}

    with open(LINEAGE_PATH, "w", newline="") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow(["dataset", "taxid", "rank", "rank_taxid", "rank_name"])
        for dataset, taxid_map in mapped.items():
            for taxid, ictv_id in taxid_map.items():
                ranks = ictv_ranks.get(ictv_id, {})
                for rank in RANKS:
                    name = ranks.get(rank)
                    if name:
                        writer.writerow([dataset, taxid, rank, name, name])

    print(f"{'rank':10s} {'ICTV':>6s} {'v2':>6s} {'v3':>6s}")
    for rank in RANKS:
        n = {d: len({r for t, i in m.items() for r in [ictv_ranks.get(i, {}).get(rank)] if r})
             for d, m in mapped.items()}
        print(f"{rank:10s} {len(totals[rank]):6d} {n['v2']:6d} {n['v3']:6d}")
    print(f"-> {TOTALS_PATH}\n-> {LINEAGE_PATH}")


if __name__ == "__main__":
    main()
