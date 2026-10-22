#!/bin/bash

OUTDIR=metadata

# File paths
SOURCE1=../uniref30
SOURCE2=../uniref30_pred_swap
REPS=uniref30_2302_db_virusdb_rep_processed.fasta
BIOLOGICAL=../virusDB/seq/uniref30_2302_db_virusdb_diff_processed.fasta
LOGAN1=../virusDB/logan/prediction_concat
LOGAN2=../virusDB/logan/prediction_concat_swap
awk -F">" '
    NR==FNR {if ($0~/>/) {swapped[$2]=$2};next}
    $0~/>/{
        split($2,name,"UniRef100_")
        {
            if ($2 in swapped) {print name[2]"\tdiff"}
            else {print name[2]"\tsame"}
        }
    }
    ' $BIOLOGICAL $REPS > tmp

awk -F"\t|/" 'NR==FNR {entry[$1]=$2;next}
    {split($NF,file,"_unrelaxed");path[file[1]]=$0;next}
    END {
        for (i in entry) {
            if (i in path) {
                print i"\t"path[i]"\t"entry[i]
            }
        }
    }
' tmp <(find $SOURCE1) <(find $SOURCE2 -name "*rank_001*.pdb") > $OUTDIR/entry_filepath_rep.tsv

awk -F"\t|/" 'NR==FNR {entry[$1]=$2;next}
    {split($NF,file,"_unrelaxed");path[file[1]]=$0;next}
    END {
        for (i in entry) {
            if (i in path) {
                print i"\t"path[i]"\t"entry[i]
            }
        }
    }
' tmp <(find $LOGAN1 -name "*rank_001*.pdb") <(find $LOGAN2 -name "*rank_001*.pdb") > $OUTDIR/logan-entry_filepath_rep.tsv
rm tmp

# length & unk count
awk 'BEGIN {OFS="\t"} $0~/^>/ {split($0,header,"UniRef100_");next} {unk=0;split($0,letter,""); for (i in letter) {if (letter[i]=="X") {unk++}}; print header[2],length($0),unk}' $REPS > $OUTDIR/entry_length_unk.tsv

 get scores
awk '
    {if ($3=="same") {

    } else {

    }
        print $0
    }
' $OUTDIR/entry_filepath_rep.tsv
getScore() {
    read -a info <<< $1
    plddt=$(jq -a '.plddt | add/length' ${info[1]} | xargs -I {} printf "%.3f" {})
    ptm=$(jq -a '.ptm' ${info[1]})
    echo "${info[0]}"$'\t'$plddt$'\t'$ptm
}
export -f getScore
sed -e 's/unrelaxed/scores/g' -e 's/pdb/json/g' -e 's/\/uniref30\//\/uniref30_pred\//' $OUTDIR/entry_filepath_rep.tsv | xargs -I {} -P 32 bash -c  'getScore "$@"' _ {} > $OUTDIR/bfvd_r1-model_plddt_ptm.tsv

sed -e 's/unrelaxed/scores/g' -e 's/pdb/json/g' $OUTDIR/logan-entry_filepath_rep.tsv | xargs -I {} -P 32 bash -c  'getScore "$@"' _ {} > $OUTDIR/logan_r1-model_plddt_ptm.tsv

for f in $(ls logan/msa); do name=${f%%.a3m}; n=$(grep -E "^>" logan/msa/$f| wc -l); echo $name$'\t'$n;done > $OUTDIR/logan-model_nmsa.tsv

awk -F"\t" 'BEGIN {OFS="\t"} NR==FNR {nmsa[$1]=$2;next} {print $1,nmsa[$1],$2,$3}' $OUTDIR/logan-model_nmsa.tsv $OUTDIR/logan_r1-model_plddt_ptm.tsv > $OUTDIR/logan-model_nmsa_plddt_ptm.tsv

awk -F"\t" 'NR==FNR {logan[$1]=$0;next} {if ($1 in logan) {print logan[$1]"\tlogan"} else {print $0"\tnologan"}} ' $OUTDIR/logan-entry_filepath_rep.tsv $OUTDIR/entry_filepath_rep.tsv > $OUTDIR/bfvd_logan-entry_filepath_rep.tsv
rm tmp

awk -F"\t" 'NR==FNR {score[$2]=$3"\t"$4;nmsa[$2]=$1;next} 
    {print $0"\t"nmsa[$1]"\t"score[$1]}
    ' <(paste <(sort -t, -k1 ../virusDB/result/uniref30_2302_virus_nmsa.csv| cut -d, -f 2) <(sort -t$'\t' -k1 $OUTDIR/bfvd_r1-model_plddt_ptm.tsv)) $OUTDIR/entry_length_unk.tsv > $OUTDIR/bfvd-model_length_unk_nmsa_plddt_ptm.tsv

awk -F"\t" 'BEGIN {OFS="\t"} NR==FNR {logan[$1]=$2"\t"$3"\t"$4;next} {if ($1 in logan) {print $1,$2,$3,logan[$1],"logan"} else {print $0,"nologan"}}' $OUTDIR/logan-model_nmsa_plddt_ptm.tsv $OUTDIR/bfvd-model_length_unk_nmsa_plddt_ptm.tsv > $OUTDIR/bfvd_logan-model_length_unk_nmsa_plddt_ptm_logan.tsv

rm $OUTDIR/bfvd_r1-model_plddt_ptm.tsv $OUTDIR/entry_length_unk.tsv

# Combine with cluster data
awk -F"\t" 'NR==FNR {split($1,r,"_unrelaxed");split($2,m,"_unrelaxed"); rep[m[1]]=r[1];cnt[r[1]]++;next} { if ($1 in rep==0) {rep[$1]="";cnt[""]=0}} {print $0"\t"rep[$1]"\t"cnt[rep[$1]]}' bfvd_cluster/clu_70_cluster.tsv $OUTDIR/bfvd-model_length_unk_nmsa_plddt_ptm.tsv > $OUTDIR/bfvd_clu70-model_length_unk_nmsa_plddt_ptm_rep_repcnt.tsv

awk -F"\t" 'NR==FNR {split($1,r,"_unrelaxed");split($2,m,"_unrelaxed"); rep[m[1]]=r[1];cnt[r[1]]++;next} { if ($1 in rep==0) {rep[$1]="";cnt[""]=0}} {print $0"\t"rep[$1]"\t"cnt[rep[$1]]}' bfvd_cluster/bfvd_logan_clu_70_cluster.tsv $OUTDIR/bfvd_logan-model_length_unk_nmsa_plddt_ptm_logan.tsv > $OUTDIR/bfvd_logan_clu70-model_length_unk_nmsa_plddt_ptm_logan_rep_repcnt.tsv
