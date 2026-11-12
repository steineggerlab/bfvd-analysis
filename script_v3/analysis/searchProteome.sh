#!/bin/bash
FASTADIR=${WKDIR}/proteome_coverage

declare -A REF_PROTEOME=(
    [canonical]=viral_proteome_canonical.fasta
)

declare -A BFVD=(
    [bfvd2]=uniref30_2302_db_virusdb_rep.fasta
    [bfvd3]=uniprot_2025_03_virus.fasta
)

mkdir -p "$FASTADIR/result"

for ref_tag in "${!REF_PROTEOME[@]}"; do
    for db_tag in "${!BFVD[@]}"; do
        ref=${REF_PROTEOME[$ref_tag]}
        db=${BFVD[$db_tag]}
        out="$FASTADIR/result/search-${ref_tag}-${db_tag}.m8"

        echo "[$(date +%F_%T)] searching $ref ($ref_tag) vs $db ($db_tag) -> $out"
        mmseqs easy-search \
            "$FASTADIR/$ref" \
            "$FASTADIR/$db" \
            "$out" \
            "$SCRATCH/${ref_tag}-${db_tag}" \
            --min-seq-id 1.0 -c 1.0
    done
done
