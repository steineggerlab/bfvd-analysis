#!/bin/bash
set -euo pipefail
if [ -z ${SLURM_SUBMIT_DIR+x} ]; then
    cd $(dirname $0)
    source ../common/env.sh
else 
    cd $SLURM_SUBMIT_DIR
    source ./common/env.sh
fi

PDB_DIR=$BENCHMARK_DIR/pdbs
UNIPROT_PDB_RAW=$BENCHMARK_DIR/uniprot_2025_03_virus_with_pdb_mapping.tsv 
BFVD_PREFIX=$BENCHMARK_DIR/uniprot_2025_03_virus
uniprot_pdb_cols=entry_review_length_pdb
PDB_STATUS=$BENCHMARK_DIR/pdb_date_status.tsv

AF_CUTOFF=$(date -d "2018-04-30" +%s)

getPDBDate() {
    pdb_id=$1
    # try to get pdb first and then cif
    date_str=$(curl -s -f "https://files.rcsb.org/download/${pdb_id}.pdb" | head -1 | grep -oE '[0-9]{2}-[A-Z]{3}-[0-9]{2}')
    if [[ -z "$date_str" ]]; then
        date_str=$(curl -s -f "https://files.rcsb.org/download/${pdb_id}.cif" | grep "_pdbx_database_status.recvd_initial_deposition_date" | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}')
    fi
    date_ts=$(date -d "$date_str" +%s)
    date_num=$(date -d "$date_str" +%d-%m-%y)
    status=""
    if [[ $date_ts -gt $(date -d "2018-04-30" +%s) ]]; then
        status="After"
    else
        status="Before"
    fi
    pdb_id=$(echo "$pdb_id" | tr '[:lower:]' '[:upper:]')
    echo -e "$pdb_id\t$date_str\t$date_num\t$status"
}

export -f getPDBDate

# Stage1: Extract Download mappable PDB structures and identify AF-untrained set
uniprot_pdb_cols=$uniprot_pdb_cols"_pdbinfo"
python parse_uniprot_retrieval.py  $UNIPROT_PDB_RAW "${BFVD_PREFIX}-${uniprot_pdb_cols}.tsv"

for file in $PDB_DIR/*.pdb; do
    # 1. Extract the date string (e.g., 17-MAY-20) from the first line
    raw_date=$(head -1 "$file" | grep -oE '[0-9]{2}-[A-Z]{3}-[0-9]{2}') # DD-MMM-YY
    filename=$(basename -- "$file" .pdb)
    
    # 2. Check if a date was found
    if [[ -n "$raw_date" ]]; then
        # 3. Convert to numeric format (e.g., 17-05-20)
        # %d = Day, %m = Month number, %y = Year (2 digits)
        file_ts=$(date -d "$raw_date" +%s)
        num_date=$(date -d "$raw_date" +%d-%m-%y)

        # 5. Compare
        if [[ $file_ts -gt $AF_CUTOFF ]]; then
            status="After"
        else
            status="Before"
        fi
        
        echo -e "$filename\t$raw_date\t$num_date\t$status" >> $PDB_STATUS
    else
        echo "$file doesn't have date information" >&2
        exit 1
    fi    
done

awk -F"\t" 'BEGIN {ORS=","} NR==FNR {d[$1]=1;next} ($1 in d==0) {print }' $PDB_STATUS  <(awk -F"\t" '{print $4}' "${BFVD_PREFIX}-${uniprot_pdb_cols}.tsv" | sort | uniq ) > pdb_na.list

for file in $PDB_DIR/*.cif; do
    # 1. Extract the date string (e.g., 2020-05-17) from the first line
    raw_date=$(grep "_pdbx_database_status.recvd_initial_deposition_date" $file | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}') # YEAR-MM-DD
    filename=$(basename -- "$file" .cif)
    
    # 2. Check if a date was found
    if [[ -n "$raw_date" ]]; then
        # 3. Convert to numeric format (e.g., 20-05-17)
        # %d = Day, %m = Month number, %y = Year (2 digits)
        file_ts=$(date -d "$raw_date" +%s)
        num_date=$(date -d "$raw_date" +%d-%m-%y)


        # 5. Compare
        if [[ $file_ts -gt $AF_CUTOFF ]]; then
            status="After"
        else
            status="Before"
        fi
        
        echo -e "$filename\t$raw_date\t$num_date\t$status" >> $PDB_STATUS
    else
        echo "$filename No_Date_Found"
        exit 1
    fi    
done

awk -F"\t" '
NR==FNR {af_train[$1]=$4;next}

{
    af_train[$4]=="" ? status="NAN" : status=af_train[$4]
    split($5,pdb_info,";")
    pdb_concat=$4";"pdb_info[1]";"pdb_info[4]";"status
    if ($1 in data) {data[$1]=data[$1]" "pdb_concat}
    else {   
        data[$1]=$1"\t"$2"\t"$3"\t"pdb_concat
    }
}

END {
    for (id in data) {
        print data[id]
    }
}
' $PDB_STATUS "${BFVD_PREFIX}-${uniprot_pdb_cols}.tsv" > "${BFVD_PREFIX}_concat-${uniprot_pdb_cols}.tsv"

# Select only UniProt accessions with AF-untrained PDB structures with highest coverage
awk -F"\t" '
BEGIN {OFS="\t"; print "UniProtID","Review","Length","PDBID","Chain","Coverage","AF_Training_Status"}
NR==FNR {af_train[$1]=$4;next}

{
    split($5,pdb_info,";")
    if ($1 in seen==0) {
        seen[$1]=$1"\t"$2"\t"$3"\t"$4
        af_check[$1]=af_train[$4]; chain[$1]=pdb_info[1]"\t"pdb_info[2]; highest_cov[$1]=pdb_info[4]
    }
    else if (af_check[$1]=="Before") {
        pass
    } else {
        if (af_train[$4]=="Before") {
            af_check[$1]=af_train[$4]; pass
        } else if (pdb_info[4]>highest_cov[$1]) {
            seen[$1]=$1"\t"$2"\t"$3"\t"$4
            af_check[$1]=af_train[$4]; chain[$1]=pdb_info[1]"\t"pdb_info[2]; highest_cov[$1]=pdb_info[4]
        }
    }
}

END {
    for (id in seen) {
        if (af_check[id]=="After") {
            print seen[id]"\t"chain[id]"\t"highest_cov[id]"\t"af_check[id]
        }
    }
}
' $PDB_STATUS "${BFVD_PREFIX}-${uniprot_pdb_cols}.tsv" > ${BFVD_PREFIX}-af_untrained_candidates.tsv

# Stage2: Filter by searching to PDB sequences before AF training date
fasta_file=$WKDIR/fasta/uniprot_2025_03_virus_1line.fasta
awk -F">|\t" '
NR==FNR {if ($0~/^>/) {header=$2; next} seqs[header]=seqs[header] $0; next}
{print ">"$1; p
    rint seqs[$1]
    if (id in to_delete==0) {
            
        }
}
' $fasta_file <(
    awk -F"\t" '$3<=1500 && $7 > 0.9 {print $1}' ${BFVD_PREFIX}-af_untrained_candidates.tsv
) > "${BFVD_PREFIX}-candidate.fasta"

$MMSEQS databases PDB "${BENCHMARK_DIR}/mmseqsdb/pdb" $SCRATCH
$MMSEQS easy-search "${BFVD_PREFIX}-candidate.fasta" "${BENCHMARK_DIR}/mmseqsdb/pdb" "${BFVD_PREFIX}-candidate_pdb_search" $SCRATCH \
    --format-output "query,target,fident,evalue,qcov,tcov,alnlen,qstart,qend,qlen,tstart,tend,tlen,bits"

# query, target, seqID, evalue, qcov, tcov, alnlen, qstart, qend, qlen, tstart, tend, tlen, bits
awk -F"\t" 'BEGIN {OFS="\t"} 
    $3 > 0.7 && $5 > 0.9 {print $2} # {print $1,$2,$3,$4,$5,$6,$10,$7}
' ${BFVD_PREFIX}-candidate_pdb_search | cut -d_ -f1 | sort | uniq > "${BFVD_PREFIX}-candidate_pdb_seq_matched.list"

# Fetch release date information if not included in PDB status file
awk -F"\t" 'NR==FNR {record[tolower($1)]=$0;next} !($1 in record) {print $1}' $PDB_STATUS "${BFVD_PREFIX}-candidate_pdb_seq_matched.list" |\
    xargs -I {} bash -c 'getPDBDate "{}"' > "${BENCHMARK_DIR}/pdb_date_status-matched_candidate_seq.tsv"

awk -F"\t" 'BEGIN {OFS="\t"} 
    $3 > 0.7 && $5 > 0.9 {print $1,$2,$3,$4,$5,$6,$10,$7}
' ${BFVD_PREFIX}-candidate_pdb_search > "${BFVD_PREFIX}-tmp.tsv"

awk -F"\t" 'BEGIN {OFS="\t"} FNR==1 { file_num++ }
file_num!=3 {status[$1]=$4; next}
file_num==3 {
    split($2,pdb_id_chain,"_"); pdb_id=toupper(pdb_id_chain[1]);
    print $0, status[pdb_id]
}' "${BENCHMARK_DIR}/pdb_date_status-matched_candidate_seq.tsv" "${BENCHMARK_DIR}/pdb_date_status.tsv" "${BFVD_PREFIX}-tmp.tsv" > "${BFVD_PREFIX}-candidate_pdb_seq_matched.tsv" 
rm "${BFVD_PREFIX}-tmp.tsv"

awk -F"\t" '
NR==FNR {if ($NF=="Before") {to_delete[$1]=1}; next}
{ 
    if ($0~/^>/) {
        split($0,h,">"); header=h[2];
        if (header in to_delete) {
            del=1
        } else {
            del=0
        }
    } else {
        if (del==0) {
            print ">"header; print $0
        }
    }
}
' ${BFVD_PREFIX}-candidate_pdb_seq_matched.tsv "${BFVD_PREFIX}-candidate.fasta" > "${BFVD_PREFIX}-candidate_filtered.fasta"

# Cluster candidate sequences to reduce the redundancy
$MMSEQS easy-cluster "${BFVD_PREFIX}-candidate_filtered.fasta" $BENCHMARK_DIR/candidate_filtered $SCRATCH \
    --min-seq-id 0.7 -c 0.8

awk -F"\t" '
    NR==FNR {bn[$1]=$1;next} $1 in bn {print $0}
' <(cut -d$'\t' -f1 ${BENCHMARK_DIR}/candidate_filtered_cluster.tsv | sort | uniq) "${BFVD_PREFIX}-af_untrained_candidates.tsv" > "${BENCHMARK_DIR}/candidate_filtered_info.tsv"