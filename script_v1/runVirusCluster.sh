#!/bin/bash


nomburgPDB=../virusDB/all_pdb/nomburg
#bfvd=<FOLDSEEK_DB>/bfvd_logan/bfvd
bfvd=<FOLDSEEK_DB>/bfvd/bfvd
out=virus_cluster
nomburg=$out/db/nomburg
db=$out/db/virus_logan

foldseek concatdbs $bfvd $nomburg $db
foldseek concatdbs $bfvd"_h" $nomburg"_h" $db"_h"
foldseek concatdbs $bfvd"_ca" $nomburg"_ca" $db"_ca"
foldseek concatdbs $bfvd"_ss" $nomburg"_ss" $db"_ss"

foldseek cluster $db $out/virus_logan_clu70 $SCRATCH -c 0.7 --remove-tmp-files false
foldseek createtsv $db $db $out/virus_logan_clu70 $out/virus_logan_clu70.tsv

mv $SCRATCH $out/tmp_logan
