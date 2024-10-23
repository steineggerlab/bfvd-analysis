#!/bin/bash

#SBATCH -t 0
#SBATCH -J virusCluster
#SBATCH -p compute
#SBATCH -w hulk
#SBATCH -c 64
#SBATCH -o log/sbatch_virus_cluster_logan.out

nomburgPDB=/home/seamustard52/virusDB/all_pdb/nomburg
#bfvd=/fast/databases/foldseek/bfvd/bfvd
bfvd=/fast/databases/foldseek/bfvd_logan/bfvd
out=/home/seamustard52/bfvd-analysis/virus_cluster
nomburg=$out/db/nomburg
#db=$out/db/virus
db=$out/db/virus_logan

#foldseek createdb $nomburgPDB $nomburg

foldseek concatdbs $bfvd $nomburg $db
foldseek concatdbs $bfvd"_h" $nomburg"_h" $db"_h"
foldseek concatdbs $bfvd"_ca" $nomburg"_ca" $db"_ca"
foldseek concatdbs $bfvd"_ss" $nomburg"_ss" $db"_ss"

foldseek cluster $db $out/virus_logan_clu70 $SCRATCH -c 0.7 --remove-tmp-files false
foldseek createtsv $db $db $out/virus_logan_clu70 $out/virus_logan_clu70.tsv

mv $SCRATCH $out/tmp_logan
