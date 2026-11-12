#!/usr/bin/env python3
"""Per-proteome sequence coverage by BFVD v2 / v3 of the UniProt viral reference proteomes

=== INPUTS ===
data/proteome_search/search-canonical-bfvd{2,3}.m8  (analysis/searchProteome.sh)
uniprot/viral_proteome_sequences.tsv                (protein_id, proteome_id, source)
=============

Input format: standard 12-col mmseqs/BLAST tabular output (query, target,
pident, alnlen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bits).
The search only reports full-length, exact, ungapped matches: coverage per
protein is binary (present at 100%, or absent), not a partial-alignment fraction.

Per-proteome coverage, in [0, 1]: covered_proteins / n_proteins

Output: analysis/results/proteome/proteome_coverage_canonical.tsv
"""

import os
import csv
from collections import defaultdict

REPO_ROOT = os.environ.get("BFVD_ROOT", os.getcwd())
SEARCH_DIR = f"{REPO_ROOT}/data/proteome_search"
PROTEINS_TSV = f"{REPO_ROOT}/uniprot/viral_proteome_sequences.tsv"
RESULTS = f"{REPO_ROOT}/analysis/results/proteome"

# fasta_set -> (source filter in viral_proteome_sequences.tsv, m8 filename prefix)
FASTA_SETS = {
    "canonical": "main",
}
VERSIONS = ["bfvd2", "bfvd3"]


def load_protein_table():
    """Return {combined_id: (proteome_id, source)}."""
    table = {}
    with open(PROTEINS_TSV, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            combined_id = f"{row['protein_id']}_{row['proteome_id']}"
            table[combined_id] = (row["proteome_id"], row["source"])
    return table


def m8_query_ids(path):
    """Return the set of distinct query IDs with >=1 hit (full-length exact match)."""
    ids = set()
    with open(path) as fh:
        for line in fh:
            ids.add(line.split("\t", 1)[0])
    return ids


def main():
    os.makedirs(RESULTS, exist_ok=True)
    protein_table = load_protein_table()

    for fasta_set, source_filter in FASTA_SETS.items():
        proteins = [(cid, pid) for cid, (pid, source) in protein_table.items()
                    if source_filter is None or source == source_filter]
        print(f"{fasta_set}: {len(proteins)} proteins")

        covered = {}
        for ver in VERSIONS:
            m8_path = f"{SEARCH_DIR}/search-{fasta_set}-{ver}.m8"
            covered[ver] = m8_query_ids(m8_path)
            print(f"  {ver}: {len(covered[ver])} proteins with a full-length exact hit")

        per_proteome = defaultdict(lambda: {
            "n_proteins": 0, "bfvd2_covered_proteins": 0, "bfvd3_covered_proteins": 0,
        })
        for combined_id, proteome_id in proteins:
            row = per_proteome[proteome_id]
            row["n_proteins"] += 1
            for ver in VERSIONS:
                if combined_id in covered[ver]:
                    row[f"{ver}_covered_proteins"] += 1

        out_path = f"{RESULTS}/proteome_coverage_{fasta_set}.tsv"
        cols = ["proteome_id", "n_proteins",
                "bfvd2_covered_proteins", "bfvd2_coverage_frac_by_count",
                "bfvd3_covered_proteins", "bfvd3_coverage_frac_by_count"]
        with open(out_path, "w") as out:
            out.write("\t".join(cols) + "\n")
            for proteome_id in sorted(per_proteome):
                row = per_proteome[proteome_id]
                bfvd2_frac_count = row["bfvd2_covered_proteins"] / row["n_proteins"]
                bfvd3_frac_count = row["bfvd3_covered_proteins"] / row["n_proteins"]
                out.write("\t".join(str(x) for x in [
                    proteome_id, row["n_proteins"],
                    row["bfvd2_covered_proteins"], round(bfvd2_frac_count, 6),
                    row["bfvd3_covered_proteins"], round(bfvd3_frac_count, 6),
                ]) + "\n")
        print(f"  wrote {out_path} ({len(per_proteome)} proteomes)")


if __name__ == "__main__":
    main()
