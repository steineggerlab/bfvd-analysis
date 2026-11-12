#!/usr/bin/env python3
"""
=== INPUTS ===
data/VMR_MSL41.v1.20260729.xlsx - ICTV Virus Metadata Resource (VMR), sheet "VMR MSL41"
  Columns used: Species (R), ICTV_ID (S), Virus name(s) (U), Virus name abbreviation(s) (V),
  Virus GENBANK accession (X)
data/vmr_accession_taxid.tsv    - accession, taxid: NCBI taxids of the VMR's GenBank accessions.
data/taxdump/                   - NCBI taxdump (via taxdump_utils.py)
=============

Shared ICTV mapping used by taxonomic_coverage.py. Taxon - ICTV species mapping strategy:
  1. accession_species_taxid: its NCBI species-rank taxid equals the species-rank taxid of
     the GenBank accessions the VMR cites for exactly one ICTV species;
  2. accession_species_taxid+name_match: the species-rank taxid is shared by several ICTV
     species (ICTV split one NCBI species), and the taxon's normalized name matches the
     species / virus name / abbreviation of exactly one of them;
  3. name_match: otherwise, the normalized name matches exactly one ICTV species.
Taxa matching none are "N/A".
"""
import csv
import os
import re
import sys
import time
from collections import Counter, defaultdict

import openpyxl

REPO_ROOT = os.environ.get("BFVD_ROOT", os.getcwd())
sys.path.insert(0, os.path.join(REPO_ROOT, "analysis"))
from taxdump_utils import TaxonomyResolver

VMR_XLSX = os.path.join(REPO_ROOT, "data/VMR_MSL41.v1.20260729.xlsx")
ACC_TAXID_CACHE = os.path.join(REPO_ROOT, "data/vmr_accession_taxid.tsv")

STEPS = ("accession_species_taxid", "accession_species_taxid+name_match", "name_match", "N/A")

# VMR accession cells: "(start..end)" position suffixes and "<segment>: " labels are stripped
POS_SUFFIX_RE = re.compile(r"\s*\([^)]*\)\s*$")
SEGMENT_LABEL_RE = re.compile(r"^[^:]{1,20}:\s*")
NORM_RE = re.compile(r"[^a-z0-9]+")

ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
BATCH_SIZE = 200
REQ_INTERVAL_SEC = 0.4
MAX_RETRIES = 5


def normalize(s):
    if not s:
        return ""
    return NORM_RE.sub(" ", s.lower()).strip()


def _text(v):
    return "" if v is None else str(v)


def clean_accession(token):
    token = POS_SUFFIX_RE.sub("", token.strip()).strip()
    token = SEGMENT_LABEL_RE.sub("", token).strip()
    return POS_SUFFIX_RE.sub("", token).strip()


def load_vmr_accessions():
    """One dict per (ICTV species isolate, accession): ictv_id, species, virus_name,
    virus_abbrev, accession."""
    wb = openpyxl.load_workbook(VMR_XLSX, read_only=True, data_only=True)
    rows = wb["VMR MSL41"].iter_rows(values_only=True)
    header = next(rows)
    assert header[17] == "Species"
    assert header[18] == "ICTV_ID"
    assert header[20] == "Virus name(s)"
    assert header[21] == "Virus name abbreviation(s)"
    assert header[23] == "Virus GENBANK accession"
    out = []
    for row in rows:
        if not row[23]:
            continue
        for token in str(row[23]).split(";"):
            acc = clean_accession(token)
            if acc:
                out.append({"ictv_id": _text(row[18]), "species": _text(row[17]),
                            "virus_name": _text(row[20]), "virus_abbrev": _text(row[21]),
                            "accession": acc})
    return out


def _fetch_batch(accs, session):
    import requests
    params = {"db": "nuccore", "id": ",".join(accs), "retmode": "json"}
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = session.post(ESUMMARY_URL, data=params, timeout=60)
            resp.raise_for_status()
            result = resp.json().get("result", {})
            # docsums are keyed by numeric UID; map back via caption / accessionversion
            out = {}
            for uid in result.get("uids", []):
                doc = result.get(uid, {})
                taxid = doc.get("taxid")
                caption = doc.get("caption") or doc.get("accessionversion", "")
                if not taxid:
                    continue
                out[caption] = str(taxid)
                if "." in caption:
                    out.setdefault(caption.split(".")[0], str(taxid))
            return out
        except (requests.RequestException, ValueError) as e:
            wait = 2 ** attempt
            print(f"  batch fetch failed (attempt {attempt}/{MAX_RETRIES}): {e}; retrying in {wait}s", file=sys.stderr)
            time.sleep(wait)
    print(f"  batch permanently failed for {len(accs)} accessions", file=sys.stderr)
    return {}


def accession_taxids(acc_rows):
    """accession -> NCBI taxid (int) for every VMR accession NCBI resolves."""
    accessions = list(dict.fromkeys(r["accession"] for r in acc_rows))
    cache = {}
    if os.path.exists(ACC_TAXID_CACHE):
        with open(ACC_TAXID_CACHE) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                cache[row["accession"]] = row["taxid"]
    todo = [a for a in accessions if a not in cache]
    if todo:
        import requests
        print(f"fetching taxids for {len(todo)} accessions from NCBI Entrez")
        session = requests.Session()
        write_header = not os.path.exists(ACC_TAXID_CACHE)
        with open(ACC_TAXID_CACHE, "a", newline="") as out_f:
            writer = csv.writer(out_f, delimiter="\t")
            if write_header:
                writer.writerow(["accession", "taxid"])
            for i in range(0, len(todo), BATCH_SIZE):
                batch = todo[i:i + BATCH_SIZE]
                acc_to_taxid = _fetch_batch(batch, session)
                for acc in batch:
                    taxid = acc_to_taxid.get(acc, "") or acc_to_taxid.get(acc.split(".")[0], "")
                    writer.writerow([acc, taxid])
                    cache[acc] = taxid
                out_f.flush()
                time.sleep(REQ_INTERVAL_SEC)
    return {acc: int(t) for acc, t in cache.items() if t}


class IctvMapper:
    """Maps NCBI taxids (plus organism names) to ICTV species."""

    def __init__(self, resolver=None):
        self.resolver = resolver or TaxonomyResolver()
        acc_rows = load_vmr_accessions()
        acc_taxid = accession_taxids(acc_rows)

        # ICTV species -> majority species-rank taxid of its cited accessions
        per_species = defaultdict(lambda: {"species": None, "species_taxids": []})
        for row in acc_rows:
            rec = per_species[row["ictv_id"]]
            rec["species"] = row["species"]
            taxid = acc_taxid.get(row["accession"])
            if taxid is None:
                continue
            sp_taxid = self.resolver.species_rank_taxid(taxid)
            if sp_taxid is not None:
                rec["species_taxids"].append(sp_taxid)
        self.sp_taxid_to_ictv = defaultdict(list)  # >1 ictv_id: ICTV split this NCBI species
        self.ictv_to_species = {}
        for ictv_id, rec in per_species.items():
            self.ictv_to_species[ictv_id] = rec["species"]
            if rec["species_taxids"]:
                top_taxid = Counter(rec["species_taxids"]).most_common(1)[0][0]
                self.sp_taxid_to_ictv[top_taxid].append(ictv_id)

        # normalized species / virus name / abbreviation -> ICTV species
        self.name_to_ictv = defaultdict(list)
        self.ictv_to_names = defaultdict(set)
        for row in acc_rows:
            ictv_id = row["ictv_id"]
            for field in ("species", "virus_name", "virus_abbrev"):
                norm = normalize(row[field])
                if not norm:
                    continue
                self.ictv_to_names[ictv_id].add(norm)
                if ictv_id not in self.name_to_ictv[norm]:
                    self.name_to_ictv[norm].append(ictv_id)

    def map(self, taxid, norm_name):
        """-> (ictv_id or None, mapping_step)"""
        sp_taxid = self.resolver.species_rank_taxid(taxid)
        candidates = self.sp_taxid_to_ictv.get(sp_taxid, []) if sp_taxid is not None else []
        if len(candidates) == 1:
            return candidates[0], "accession_species_taxid"
        if len(candidates) > 1:
            name_matches = [c for c in candidates if norm_name in self.ictv_to_names.get(c, ())]
            if len(name_matches) == 1:
                return name_matches[0], "accession_species_taxid+name_match"
        name_hits = self.name_to_ictv.get(norm_name, [])
        if len(name_hits) == 1:
            return name_hits[0], "name_match"
        return None, "N/A"

    def map_taxa(self, taxa, label):
        """taxa: {taxid: organism name} -> {taxid: ictv_id} for mapped taxa (input order kept)."""
        out, steps = {}, Counter()
        for taxid, taxname in taxa.items():
            ictv_id, step = self.map(taxid, normalize(taxname))
            steps[step] += 1
            if ictv_id is not None:
                out[taxid] = ictv_id
        print(f"{label}: {len(taxa)} taxids; " +
              ", ".join(f"{s} {steps[s]} ({100 * steps[s] / len(taxa):.2f}%)" for s in STEPS))
        return out
