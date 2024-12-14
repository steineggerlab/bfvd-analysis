#!/bin/bash

OUTDIR=/home/seamustard52/bfvd-analysis/3d-beacons-bfvd/metadata

METADATA=/home/seamustard52/bfvd-analysis/metadata/bfvd_logan-model_length_unk_nmsa_plddt_ptm_logan.tsv
TAXDATA=/home/seamustard52/bfvd-analysis/distribute/bfvd_taxid_rank_scientificname_lineage.tsv
LENDATA=/home/seamustard52/virusDB/distribute/uniref30_2302_virus/uniref30_2302_virus-rep_length.tsv
UNIPARC=/home/seamustard52/bfvd-analysis/3d-beacons-bfvd/uniparc.list

awk -F"\t" '
    BEGIN {OFS="\t"}
    NR==FNR {len[$1]=$2;next}
    {
        split($1,acc,"_");
        if (len[acc[1]]==$2) {
            last[acc[1]]=$2;
            start[$1]=1; end[$1]=$2;
        }
        else {
            if (acc[1] in last==0) {
                last[acc[1]]=$2
                start[$1]=1;end[$1]=$2
            } else {
                start[$1]=last[acc[1]]+1;end[$1]=last[acc[1]]+$2
                last[acc[1]]+=$2;
            }
        }
        print $1,acc[1],start[$1],end[$1],len[acc[1]],$5
    }
' $LENDATA $METADATA > $OUTDIR/tmp1

awk -F"\t" '
    BEGIN {OFS="\t"}
    NR==FNR {split($1,acc,"_"); 
    tax[acc[1]]=$2"\t"$4;next}
    {if ($2 in tax)
        {print $0,tax[$2]}
    else {
        print $0,"",""
    }
    }
' $TAXDATA $OUTDIR/tmp1 > $OUTDIR/tmp2

awk -F"\t" ' 
    BEGIN {OFS="\t"}
    NR==FNR {uniparc[$1]=$1;next}
    {
        if ($2 in uniparc) {
            print $0,"UNIPARC"
        } else {
            print $0,"UNIPROT"
        }
    }
' $UNIPARC $OUTDIR/tmp2 > $OUTDIR/bfvd_logan-entry_acc_start_end_len_plddt_taxid_organism_src.tsv

rm $OUTDIR/tmp*
#rm $OUTDIR/tmp1
#rm $OUTDIR/tmp2
#awk -F"\t" 'BEGIN {OFS="\t"} {split($1,entry,"_");print entry[1],$2,$4}' $TAXDATA | sort | uniq > $OUTDIR/entry_taxid_scientific.tsv
