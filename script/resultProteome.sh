#!/bin/bash

PROTEOMEDIR=/home/seamustard52/bfvd-analysis/proteome_cover
FASTADIR=$PROTEOMEDIR/fasta
PREDICTION=$PROTEOMEDIR/prediction
SEARCH=$PROTEOMEDIR/proteome_bfvd_logan.tsv
METADIR=/home/seamustard52/bfvd-analysis/metadata
RESULTDIR=/home/seamustard52/bfvd-analysis/result

out=proteome-entry_uniprotid_virus_proteome.tsv
#for file in $(ls $FASTADIR)
#do
#    #virus=${file%%.fasta}
#    grep ">" $FASTADIR/$file | grep -Eo "\|.*\||OS=.*OX=[0-9]+" | 
#        sed -e 's/OX=/\t/g' -e 's/OS=/\t/g' -e 's/sp|//g' -e 's/|//g' | 
#        awk -F"\t" 'FNR%2==1{printf $0"\t" } NR%2==0 {print $2"\t"$3}' >> tmp.out
#done
#
#awk -F"\t" 'NR==FNR {data[$1]=$0;next} {split($1,path,"/");split(path[length(path)],id,"_");print path[length(path)]"\t"data[id[1]]}' tmp.out <(ls $PREDICTION/*.done.txt | cut -d . -f 1) > $METADIR/$out
#rm tmp.out
#
#search_best=bfvd_logan_proteome_search_best-entry_uniprotid_virus_proteome_query_etarget_eval_tmtarget_qtm.tsv
#awk -F"\t" '
#    BEGIN {OFS="\t"} 
#    NR==FNR {split($1,name,"_unrelaxed"); 
#    if (name[1] in entry==0) {entry[name[1]]=$1;e_tar[name[1]]=$2;e_val[name[1]]=$3;tm_tar[name[1]]=$2;tm_val[name[1]]=$9;next} 
#    if ($3 < e_val[name[1]]) {e_val[name[1]]=$3;e_tar[name[1]]=$2} 
#    if ($9 > tm_val[name[1]]) {tm_val[name[1]]=$9;tm_tar[name[1]]=$2}; next} 
#    {print $0,entry[$1],e_tar[$1],e_val[$1],tm_tar[$1],tm_val[$1]} 
#    ' $SEARCH $METADIR/$out > $RESULTDIR/$search_best

### Pick by best e-value
#search_best=bfvd_logan_proteome_search_best-entry_uniprotid_virus_proteome_query_etarget_eval_qtm.tsv
#awk -F"\t" '
#    BEGIN {OFS="\t"} 
#    NR==FNR {split($1,name,"_unrelaxed"); 
#    if (name[1] in entry==0) {entry[name[1]]=$1;e_tar[name[1]]=$2;e_val[name[1]]=$3;tm_val[name[1]]=$9;next} 
#    if ($3 < e_val[name[1]]) {e_val[name[1]]=$3;e_tar[name[1]]=$2;tm_val[name[1]]=$9};next} 
#    {print $0,entry[$1],e_tar[$1],e_val[$1],tm_val[$1]} 
#    ' $SEARCH $METADIR/$out > $RESULTDIR/$search_best

### Check the e-value cutoff 0.1 and save the highest qTM
search_best=bfvd_logan_proteome_search_best_eval01-entry_uniprotid_virus_proteome_query_target_eval_qtm.tsv
awk -F"\t" '
    BEGIN {OFS="\t"} 
    NR==FNR {split($1,name,"_unrelaxed"); 
    if (name[1] in entry==0 && $3<0.1) {entry[name[1]]=$1;tar[name[1]]=$2;e_val[name[1]]=$3;tm_val[name[1]]=$9;next} 
    if (name[1] in entry && $3<0.1 && $9 > tm_val[name[1]]) {tm_val[name[1]]=$9;e_val[name[1]]=$3}; next} 
    {print $0,entry[$1],tar[$1],e_val[$1],tm_val[$1]} 
    ' $SEARCH $METADIR/$out > $RESULTDIR/$search_best


#getScore() {
#    read -a info <<< $1
#    plddt=$(jq -a '.plddt | add/length' ${info[1]} | xargs -I {} printf "%.3f" {})
#    ptm=$(jq -a '.ptm' ${info[1]})
#    echo "${info[0]}"$'\t'$plddt$'\t'$ptm
#}
#export -f getScore
#
#for i in $(ls $PREDICTION/*rank_001*.json)
#do
#    path=${i%%_scores*}
#    id=${path##*/} 
#    plddt=$(jq -a '.plddt | add/length' $i)
#    echo $id$'\t'$plddt
#done > $RESULTDIR/proteome-entry_avgplddt.tsv
