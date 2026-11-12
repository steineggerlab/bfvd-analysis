#!/usr/bin/env python3
"""
Shared helpers for walking the NCBI taxonomy dump. Used to map taxIDs (ICTV-side VMR
accession -> species taxid, and BFVD-side organism taxid -> species taxid)
and calculate coverage (full multi-rank lineage per taxid).
"""
import os

REPO_ROOT = os.environ.get("BFVD_ROOT", os.getcwd())
TAXDUMP_DIR = os.path.join(REPO_ROOT, "data/taxdump")

NODES_DMP = os.path.join(TAXDUMP_DIR, "nodes.dmp")
MERGED_DMP = os.path.join(TAXDUMP_DIR, "merged.dmp")
DELNODES_DMP = os.path.join(TAXDUMP_DIR, "delnodes.dmp")
NAMES_DMP = os.path.join(TAXDUMP_DIR, "names.dmp")


def _split_dmp_line(line):
    # taxdump .dmp fields are separated by "\t|\t", lines end in "\t|"
    return [f.strip() for f in line.rstrip("\n").rstrip("|").split("|")]


def load_nodes():
    """taxid (int) -> (parent_taxid (int), rank (str))"""
    nodes = {}
    with open(NODES_DMP) as f:
        for line in f:
            fields = _split_dmp_line(line)
            taxid = int(fields[0])
            parent = int(fields[1])
            rank = fields[2]
            nodes[taxid] = (parent, rank)
    return nodes


def load_merged():
    """old_taxid (int) -> new_taxid (int)"""
    merged = {}
    with open(MERGED_DMP) as f:
        for line in f:
            fields = _split_dmp_line(line)
            merged[int(fields[0])] = int(fields[1])
    return merged


def load_delnodes():
    """set of deleted taxids (int)"""
    dels = set()
    with open(DELNODES_DMP) as f:
        for line in f:
            fields = _split_dmp_line(line)
            dels.add(int(fields[0]))
    return dels


def load_scientific_names():
    """taxid (int) -> scientific name (str)"""
    names = {}
    with open(NAMES_DMP) as f:
        for line in f:
            fields = _split_dmp_line(line)
            if fields[3] != "scientific name":
                continue
            names[int(fields[0])] = fields[1]
    return names


class TaxonomyResolver:
    """Resolves a raw NCBI taxid to its species-rank ancestor taxid.

    Handles merged taxids (redirect to the surviving id) and deleted taxids
    (unresolvable -> None). Walks the nodes.dmp parent chain upward until a
    node with rank == "species" is found; taxids at or above genus rank (no
    species-rank ancestor exists on the path to root) resolve to None.
    """

    def __init__(self, nodes=None, merged=None, delnodes=None):
        self.nodes = nodes if nodes is not None else load_nodes()
        self.merged = merged if merged is not None else load_merged()
        self.delnodes = delnodes if delnodes is not None else load_delnodes()

    def canonicalize(self, taxid):
        """Resolve a possibly-merged/deleted taxid to a current taxid, or None."""
        seen = set()
        while taxid in self.merged and taxid not in seen:
            seen.add(taxid)
            taxid = self.merged[taxid]
        if taxid in self.delnodes or taxid not in self.nodes:
            return None
        return taxid

    def species_rank_taxid(self, taxid, max_steps=64):
        """Return the species-rank ancestor taxid of `taxid`, or None if
        `taxid` can't be resolved or has no species-rank ancestor."""
        taxid = self.canonicalize(taxid)
        if taxid is None:
            return None
        current = taxid
        for _ in range(max_steps):
            parent, rank = self.nodes.get(current, (None, None))
            if rank == "species":
                return current
            if parent is None or parent == current:
                return None
            current = parent
        return None

    def rank_lineage(self, taxid, target_ranks, max_steps=64):
        """Walk `taxid` up to the root, returning {rank: ancestor_taxid} for
        every rank in `target_ranks` found along the way (ranks not present
        on this lineage are absent from the result).
        Returns {} if `taxid` can't be resolved."""
        taxid = self.canonicalize(taxid)
        if taxid is None:
            return {}
        target_ranks = set(target_ranks)
        found = {}
        current = taxid
        for _ in range(max_steps):
            parent, rank = self.nodes.get(current, (None, None))
            if rank in target_ranks and rank not in found:
                found[rank] = current
            if len(found) == len(target_ranks):
                break
            if parent is None or parent == current:
                break
            current = parent
        return found
