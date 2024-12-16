#!/bin/bash

VIRUSKEYS=/home/mmirdit/tmp/bfvd/uniref30_2302_virus-members_keys.tsv
DISTRIBUTE=distribute
LCADIR=$DISTRIBUTE/lca
UNIREF30=/fast/databases/colabfold_db_all/uniref30_2302_db
#BFVD=/fast/databases/foldseek/bfvd_logan/bfvd
BFVD=/fast/databases/foldseek/bfvd/bfvd

## Update mapping file
awk -F"\t" 'NR==FNR {split($1,uniref100,"_");tax[uniref100[2]]=$2;next} {print $1"\t"tax[$2]}' $LCADIR/uniref100_2302_accession_mapping.tsv $UNIREF30"_seq_h_accessions.tsv" > $LCADIR/uniref30_2302_db_mapping
mv $UNIREF30"_mapping" $UNIREF30"_mapping_orig"
ln -s $LCADIR/uniref30_2302_db_mapping $UNIREF30"_mapping"
##lca
mmseqs createsubdb $VIRUSKEYS $UNIREF30"_aln" $LCADIR/uniref30_2302_db_virus_aln
mmseqs lca $UNIREF30 $LCADIR/uniref30_2302_db_virus_aln $LCADIR/uniref30_2302_db_virus_lca
mmseqs createtsv $UNIREF30 $LCADIR/uniref30_2302_db_virus_lca $LCADIR/uniref30_2302_virus_lca.tsv
sed -i 's/UniRef100_//g' $LCADIR/uniref30_2302_virus_lca.tsv

#tmp alignment for lineage
awk -F"\t" '{print $1"\t"$1}' $BFVD".lookup" > $LCADIR/tmp
mmseqs tsv2db $LCADIR/tmp $LCADIR/tmp_db
mmseqs addtaxonomy $BFVD $LCADIR/tmp_db $LCADIR/bfvd_taxlineage --tax-lineage 1 
mmseqs createtsv $BFVD $LCADIR/bfvd_taxlineage $LCADIR/bfvd_taxlineage.tsv
awk -F"\t" 'BEGIN {OFS="\t"} {print $1".pdb",$3,$4,$5,$6} ' $LCADIR/bfvd_taxlineage.tsv > $DISTRIBUTE/bfvd_taxid_rank_scientificname_lineage.tsv
