#!/bin/bash

PROTEOMEDIR=/home/seamustard52/bfvd-analysis/proteome_cover

awk '$0~/^>/{split($0,head,"|");header=">"head[2];next} $0!~/^>/ {seq[header]=seq[header]$0;next} END {for (i in seq) {print i; print seq[i]} }' $PROTEOMEDIR/fasta/*.fasta > $PROTEOMEDIR/concat.fasta

awk '
    FNR%2==1{header[FNR]=$0;next} 
    {
        h=header[FNR-1]
        seq=$0
        if ((FNR%2==0))
        {
            print h
            print seq
        }

        if ((FNR%2==0) && (length($0)>1500))
        {
            nchop=(int(length(seq)/1501)+1);lenchop=int(length(seq)/nchop);
        
            nth=1;
            while (nth < nchop) {
                chopseq=substr(seq, lenchop*(nth-1) +1, lenchop)
                print h"_"nth;
                print chopseq;
                nth++
            }
            chopseq=substr(seq,lenchop*(nth-1)+1, length(seq))
            print h"_"nth;
            print chopseq;
        }

    }' $PROTEOMEDIR/concat.fasta > $PROTEOMEDIR/concat_processed.fasta

#rm $PROTEOMEDIR/concat.fasta
