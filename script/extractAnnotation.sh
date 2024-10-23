#!/bin/bash

BFVD=/home/seamustard52/bfvd-analysis
RESDIR=$BFVD/result

ANNOTDIR=/home/seamustard52/virusDB/application/Small_ACCs

#ls $ANNOTDIR/gac6*/INHERIT_phages/

for i in $(ls $ANNOTDIR/gac6*/INHERIT_phages/*.fasta);
do
    dir=${i%/INHERIT_phages*}
    file=${i##*/}
    contig=${file%.fasta}

    cds=$(sed -ne 's/^.*CDSs: //p' $dir/initial_polished/bakta/$contig.txt )
    # not annotated with bakta
    bakta_hypothetical=$(sed -ne 's/^.*hypotheticals: //p' $dir/initial_polished/bakta/$contig.txt )
    afdb=$(grep -c 'hypothetical' $dir/INHERIT_phages/foldseek_hits/$contig"_besthits.tab")
    bfvd_logan=$(grep -c 'hypothetical' /home/seamustard52/bfvd-analysis/phage_annotation/logan/$contig"_besthits.tab")
    #AFDB + BFVD
    afdb_bfvd=$(awk -F"\t" '$1~/hypothetical/ {split($1,file,".pdb");print file[1]}'  $dir/INHERIT_phages/foldseek_hits/$contig"_besthits.tab" /home/seamustard52/bfvd-analysis/phage_annotation/logan/$contig"_besthits.tab"| sort -u | wc -l)
    nomburg=$(grep -c 'hypothetical' $dir/INHERIT_phages/foldseek_hits_nomburg/$contig"_besthits.tab")

    echo "contig_"$contig$'\t'$cds$'\t'$bakta_hypothetical$'\t'$afdb$'\t'$bfvd_logan$'\t'$afdb_bfvd$'\t'$nomburg
done

