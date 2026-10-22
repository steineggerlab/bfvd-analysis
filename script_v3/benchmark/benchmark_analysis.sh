#!/bin/bash
set -euo pipefail
# set -x
if [ -z ${SLURM_SUBMIT_DIR+x} ]; then
    cd $(dirname $0)
    source ../common/env.sh
else 
    cd $SLURM_SUBMIT_DIR
    source ./common/env.sh
fi

if false; then
benchmark/benchmark_prediction.sh
benchmark/benchmark_prediction.sh -t
benchmark/benchmark_prediction.sh -r 12
benchmark/benchmark_prediction.sh -t -r 12

benchmark/benchmark_prediction.sh -m logan
benchmark/benchmark_prediction.sh -m logan -t
benchmark/benchmark_prediction.sh -m logan -r 12
benchmark/benchmark_prediction.sh -m logan -r 12 -t
fi

if false; then
for dir in $BENCHMARK_DIR/prediction/*; do
    dirname=$(basename "${dir}")
    base=$(echo "${dirname}" | cut -d'_' -f1)
    recycle=$(echo "${dirname}" | cut -d'_' -f2)
    template=$(echo "${dirname}" | awk -F'_' '{print ($3 == "" ? "notemplate" : $3)}')
    for scorefile in $dir/*_scores*.json; do
       fn=$(basename "${scorefile}")
       entry=${fn%%_scores*}
       rank=$(echo "${fn}" | grep -oE 'rank_[0-9]+' | cut -d'_' -f2)
       model=$(echo "${fn}" | grep -oE 'model_[^_]+'| cut -d'_' -f2)
       avg_plddt=$(jq '.plddt | add / length' "${scorefile}")
       echo -e "${entry}\t${rank}\t${model}\t${avg_plddt}\t${base}\t${recycle}\t${template}" >> $METADIR/benchmark-entry_rank_model_plddt_base_recycle_template.tsv
   done
done

for dir in $BENCHMARK_DIR/prediction/*; do
    dirname=$(basename "${dir}")
    base=$(echo "${dirname}" | cut -d'_' -f1)
    recycle=$(echo "${dirname}" | cut -d'_' -f2)
    template=$(echo "${dirname}" | awk -F'_' '{print ($3 == "" ? "notemplate" : $3)}')

    mapping_file=$BENCHMARK_DIR/candidate_filtered_info.tsv
    pdbdir=$BENCHMARK_DIR/pdbs
    for pdbfile in $dir/*_unrelaxed_rank_*.pdb; do
        fn=$(basename "${pdbfile}")
        model=${fn%%_unrelaxed_rank*}
        rank=$(echo "${fn}" | grep -oE 'rank_[0-9]+' | cut -d'_' -f2)

        read -r pdbid chain < <(awk -F"\t" -v model="$model" '$1 == model {print $4, $5}' "$mapping_file")
        ref_pdb=""
        if [ -f "$pdbdir/${pdbid}.pdb" ]; then
            ref_pdb=$pdbdir/${pdbid}.pdb
        elif [ -f "$pdbdir/${pdbid}.cif" ]; then
            ref_pdb=$pdbdir/${pdbid}.cif
        else
            echo "Reference PDB ${pdbid} for ${model} not found" >> $BENCHMARK_DIR/missing_reference_pdbs.txt
            continue
            # exit 1
        fi
        foldseek easy-search $pdbfile $ref_pdb $BENCHMARK_DIR/foldseek_aln/"${model}_${rank}-${pdbid}_${chain}-${base}_${recycle}_${template}.m8" $BENCHMARK_DIR/foldseek_tmp \
            --format-output "query,target,qtmscore,ttmscore,alntmscore,evalue,lddt,qstart,qend,qlen,tstart,tend,tlen,alnlen" --exact-tmscore 1  --exhaustive-search 1
    done
done
fi
for file in $BENCHMARK_DIR/foldseek_aln/*.m8; do
    entry=$(basename "${file}" .m8)
    awk '
    BEGINFILE {
        split(FILENAME, a, "/");fn=a[length(a)];split(fn,b,"-");
        split(b[1], m, "_"); split(b[2], gt, "_"); split(b[3],options,"."); split(options[1], options, "_");
        model=m[1]; rank=m[2]; pdbid=gt[1]; chain=gt[2]; base=options[1]; recycle=options[2]; template=options[3];
        record=model"\t"rank"\t"base"\t"recycle"\t"template"\t"b[2]
    }
    { # Only retain the matched PDB in UniProt record
        if ($2==b[2] || $2==pdbid) {
            # transform $3,$4,$7 to numeric with 3 decimal places
            $3=sprintf("%.3f", $3); $4=sprintf("%.3f", $4); $7=sprintf("%.3f", $7);
            record=record"\t"$3"\t"$4"\t"$7"\t"$14
        }
    }
    ENDFILE {
        split(record, fields, "\t")
        if (length(fields) < 6) {
            record=record"\tNA\tNA\tNA\tNA"
        }
        print record
    }
    ' "$file" >> $METADIR/benchmark_aln2pdb-model_rank_base_recycle_template_pdb_qtmscore_ttmscore_lddt_alnlen.tsv
done

awk -F"\t" '
BEGINFILE {file_idx++}
{
    if (file_idx == 1) {rec[$1]=$2"\t"$3"\t"$4}
    else if (file_idx == 2) {
        $4=sprintf("%.3f", $4);
        scores[$1"\t"$2"\t"$5"\t"$6"\t"$7]=$1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"rec[$1]"\t"$4
    }
    else if (file_idx == 3) {
        if (scores[$1"\t"$2"\t"$3"\t"$4"\t"$5] != "") {
        print scores[$1"\t"$2"\t"$3"\t"$4"\t"$5]"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10
        }
    }
}
' $METADIR/benchmark-entry_nbasemsa_ntemplate_nloganmsa.tsv \
$METADIR/benchmark-entry_rank_model_plddt_base_recycle_template.tsv \
$METADIR/benchmark_aln2pdb-model_rank_base_recycle_template_pdb_qtmscore_ttmscore_lddt_alnlen.tsv > $METADIR/benchmark_combined-entry_rank_model_base_recycle_template_nbasemsa_ntemplate_nloganmsa_plddt_pdb_qtmscore_ttmscore_lddt_alnlen.tsv
