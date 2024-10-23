#!/bin/bash
BFVD=/home/seamustard52/bfvd-analysis
RESDIR=$BFVD/result
METADIR=$BFVD/metadata

### Result for BFVD cluster
CLUDIR=$BFVD/bfvd_cluster
#cluout1=bfvd_cluster-cov_nonsing_sing_total.csv
#echo "coverage,non-singleton,singleton,total" > $RESDIR/$cluout1
#for i in $(ls $CLUDIR/clu_*.tsv)
#do
#    cov=$(echo ${i##*clu_}| cut -d _ -f1)
#    result=$(awk -F"\t" '{rep[$1]++} END {for (r in rep) {if (rep[r]==1) {singleton++} else {nonsingleton++};total++}; print nonsingleton","singleton","total}' $i)
#    echo $cov","$result >> $RESDIR/$cluout1
#done

#cluout2=bfvd_logan_cluster-cov_nonsing_sing_total.csv
#echo "coverage,non-singleton,singleton,total" > $RESDIR/$cluout2
#for i in $(ls $CLUDIR/bfvd_logan_clu*.tsv)
#do
#    cov=$(echo ${i##*bfvd_logan_clu_}| cut -d _ -f1)
#    result=$(awk -F"\t" '{rep[$1]++} END {for (r in rep) {if (rep[r]==1) {singleton++} else {nonsingleton++};total++}; print nonsingleton","singleton","total}' $i)
#    echo $cov","$result >> $RESDIR/$cluout2
#done

### Result for Joint cluster
#JOINTCLU=$BFVD/virus_cluster/virus_clu70.tsv
JOINTCLU=$BFVD/virus_cluster/virus_logan_clu70.tsv
#awk -F"\t" ' BEGIN {OFS="\t"}
#    {if ($1 in rep==0) {rep[$1]=0;bfvd[$1]=0;nomburg[$1]=0}} 
#    {rep[$1]++ ;if ($2~/unrelaxed/) {bfvd[$1]++} else {nomburg[$1]++};next} 
#    END {for (r in rep) {print r,bfvd[r],nomburg[r]}}
#    ' $JOINTCLU > $RESDIR/virus_logan_clu_coverage-rep_bfvd_nomburg.tsv

### Result for Annotation
ANNOTDIR=/home/seamustard52/virusDB/application/Small_ACCs
#for i in $(ls $ANNOTDIR/gac6*/INHERIT_phages/*.fasta);
#do
#    dir=${i%/INHERIT_phages*}
#    file=${i##*/}
#    contig=${file%.fasta}
#
#    cds=$(sed -ne 's/^.*CDSs: //p' $dir/initial_polished/bakta/$contig.txt )
#    # not annotated with bakta
#    bakta_hypothetical=$(sed -ne 's/^.*hypotheticals: //p' $dir/initial_polished/bakta/$contig.txt )
#    afdb=$(grep -c 'hypothetical' $dir/INHERIT_phages/foldseek_hits/$contig"_besthits.tab")
#    bfvd_logan=$(grep -c 'hypothetical' /home/seamustard52/bfvd-analysis/phage_annotation/logan/$contig"_besthits.tab")
#    #AFDB + BFVD
#    afdb_bfvd=$(awk -F"\t" '$1~/hypothetical/ {split($1,file,".pdb");print file[1]}'  $dir/INHERIT_phages/foldseek_hits/$contig"_besthits.tab" /home/seamustard52/bfvd-analysis/phage_annotation/logan/$contig"_besthits.tab"| sort -u | wc -l)
#    nomburg=$(grep -c 'hypothetical' $dir/INHERIT_phages/foldseek_hits_nomburg/$contig"_besthits.tab")
#
#    echo "contig_"$contig$'\t'$cds$'\t'$bakta_hypothetical$'\t'$afdb$'\t'$bfvd_logan$'\t'$afdb_bfvd$'\t'$nomburg
#done > $RESDIR/bfvd_logan-total_num_hypoths.txt

### Major1: fragment analysis
#awk -F"\t" 'NR==FNR {fragment[$1]=$1;next} 
#    {split($1,entry,"_") ;
#     if ($NF==1) {singleton=1} else {singleton=0}
#     if (entry[1] in fragment) {fragmented=1} else {fragmented=0};
#     print $1"\t"entry[1]"\t"singleton"\t"fragmented}' $METADIR/BFVD_accessions_with_status_fragment.txt $METADIR/bfvd_logan_clu70-model_length_unk_nmsa_plddt_ptm_logan_rep_repcnt.tsv > $RESDIR/bfvd_logan_fragment-uniprotid_entry_singleton_fragment.tsv

#awk -F"\t" '
#    $3==0 && $4==0 {nosing_nofrag++;next} 
#    $3==1 && $4==0 {sing_nofrag++;next}
#    $3==0 && $4==1 {nosing_frag++;next}
#    $3==1 && $4==1 {sing_frag++;next}
#    END {print nosing_nofrag,nosing_frag,sing_nofrag,sing_frag}
#' $RESDIR/bfvd_logan_fragment-uniprotid_entry_singleton_fragment.tsv


#awk -F"\t" 'BEGIN {OFS="\t"} NR==FNR {fragment[$1]=$1;next} 
#    {if ($2 in fragment) {print $1,$2,$4,1} 
#    else {print $1,$2,$4,0}}
#' $METADIR/BFVD_accessions_with_status_fragment.txt <(awk -F"\t" '
#    NR==FNR {
#        if ($1 in entry==0) {entry[$1]=$1;bfvd[$1]=0;nomburg[$1]=0};
#        if ($2~/unrelaxed/) {bfvd[$1]++} else {nomburg[$1]++}; next
#    }
#    NR!=FNR { 
#        if ($2!~/unrelaxed/) {next};
#            split($2,entry,"_unrelaxed"); split(entry[1],id,"_")
#        if (bfvd[$1]>0 && nomburg[$1]>0) {print entry[1]"\t"id[1]"\t"$2"\tShared"}
#        else if (bfvd[$1]>0 && nomburg[$1]==0) {print entry[1]"\t"id[1]"\t"$2"\tBFVD"}
#        else if (bfvd[$1]==0 && nomburg[$1]>0) {print entry[1]"\t"id[1]"\t"$2"\tNomburg"}
#        }' $JOINTCLU $JOINTCLU ) > $RESDIR/bfvd_logan_fragment-uniprotid_entry_viralcluster_fragment.tsv

#awk -F"\t" '
#    $3=="Shared" && $4==0 {shared_nofrag++;next} 
#    $3=="BFVD" && $4==0 {uniq_nofrag++;next}
#    $3=="Shared" && $4==1 {shared_frag++;next}
#    $3=="BFVD" && $4==1 {uniq_frag++;next}
#    END {print shared_nofrag,shared_frag,uniq_nofrag,uniq_frag}
#' $RESDIR/bfvd_logan_fragment-uniprotid_entry_viralcluster_fragment.tsv

### DB Coverage result
#awk -F"\t" 'BEGIN {OFS="\t"} 
#    NR==FNR {split($1,entry,"_unrelaxed");result[entry[1]]=$1"\t"$2;number[entry[1]]=$4"\t"$3"\t"$5;next} 
#    { if ($1 in result) {print $1,result[$1],$2,number[$1],"PDB100"}
#      else {print $1,"","",$2,0,"","","PDB100"}
#    }
#' <(awk -F"\t" '!($1 in result){result[$1]=$1"\t"$2"\t"$3"\t"$4"\t"$6;qtm[$1]=$4;next} qtm[$1] < $4 {result[$1]=$1"\t"$2"\t"$3"\t"$4"\t"$6;qtm[$1]=$4} END {for (r in result) {print result[r]}}' db_cover/bfvd_logan_pdb100_aln) $METADIR/bfvd_logan-model_length_unk_nmsa_plddt_ptm_logan.tsv > $RESDIR/bfvd_logan_dbsearch-entry_query_target_qlen_qtm_eval_lddt_db.tsv
#
#awk -F"\t" 'BEGIN {OFS="\t"} 
#    NR==FNR {split($1,entry,"_unrelaxed");result[entry[1]]=$1"\t"$2;number[entry[1]]=$4"\t"$3"\t"$5;next} 
#    { if ($1 in result) {print $1,result[$1],$2,number[$1],"AFDB50"}
#      else {print $1,"","",$2,0,"","","AFDB50"}
#    }
#' <(awk -F"\t" '!($1 in result){result[$1]=$1"\t"$2"\t"$3"\t"$4"\t"$6;qtm[$1]=$4;next} qtm[$1] < $4 {result[$1]=$1"\t"$2"\t"$3"\t"$4"\t"$6;qtm[$1]=$4} END {for (r in result) {print result[r]}}' db_cover/bfvd_logan_afdb50_aln) $METADIR/bfvd_logan-model_length_unk_nmsa_plddt_ptm_logan.tsv >> $RESDIR/bfvd_logan_dbsearch-entry_query_target_qlen_qtm_eval_lddt_db.tsv
#
#awk -F"\t" 'BEGIN {OFS="\t"} NR==FNR {nmem[$1]=$NF;next} {print $1,nmem[$1],$5,$NF}' $METADIR/bfvd_logan_clu70-model_length_unk_nmsa_plddt_ptm_logan_rep_repcnt.tsv $RESDIR/bfvd_logan_dbsearch-entry_query_target_qlen_qtm_eval_lddt_db.tsv > $RESDIR/bfvd_logan_dbsearch-entry_nmem_qtm_db.tsv


# ETC: pLDDT comparsion between consensus vs biological structures
#awk -F"\t" 'NR==FNR{score[$1]=$0;next} $3=="diff" {print score[$1]}' <(awk -F"\t|," 'NR==FNR {new[$1]=$5;next} {old[$1]=$3;next} END {for (i in new) {print i"\t"old[i]"\t"new[i]}}' metadata/bfvd-model_length_unk_nmsa_plddt_ptm.tsv ../virusDB/result/uniref30_2302_virus_r1.csv) metadata/entry_filepath_rep.tsv > result/diff_oldplddt_newplddt.tsv
