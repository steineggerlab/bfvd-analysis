#!/bin/bash

DISTRIBUTE=/home/seamustard52/bfvd-analysis/distribute
### Make sure to copy taxonomy and mapping file in foldseekdb before distribute
tar --owner=0 --group=0 -czvf $DISTRIBUTE/bfvd_foldseekdb.tar.gz -C /fast/databases/foldseek/bfvd_logan . 
tar --owner=0 --group=0 -czvf $DISTRIBUTE/bfvd_foldcompdb.tar.gz -C /fast/databases/foldcomp/bfvd_logan .
tar --owner=0 --group=0 -czvf $DISTRIBUTE/bfvd.tar.gz -C /home/seamustard52/bfvd-analysis/pdbs_bfvd_logan .

aws s3 cp $DISTRIBUTE/bfvd.tar.gz s3://bfvd --endpoint-url https://4dfc310765e5ee60c351741d16a92839.r2.cloudflarestorage.com
aws s3 cp $DISTRIBUTE/bfvd_foldseekdb.tar.gz s3://bfvd --endpoint-url https://4dfc310765e5ee60c351741d16a92839.r2.cloudflarestorage.com
aws s3 cp $DISTRIBUTE/bfvd_foldcompdb.tar.gz s3://bfvd --endpoint-url https://4dfc310765e5ee60c351741d16a92839.r2.cloudflarestorage.com

#tar --owner=0 --group=0 -czvf $DISTRIBUTE/bfvd_base_foldseekdb.tar.gz -C /fast/databases/foldseek/bfvd . 
#tar --owner=0 --group=0 -czvf $DISTRIBUTE/bfvd_base_foldcompdb.tar.gz -C /fast/databases/foldcomp/bfvd .
#tar --owner=0 --group=0 -czvf $DISTRIBUTE/bfvd_base.tar.gz -C /home/seamustard52/bfvd-analysis/pdbs .

#aws s3 cp $DISTRIBUTE/bfvd_base.tar.gz s3://bfvd --endpoint-url https://4dfc310765e5ee60c351741d16a92839.r2.cloudflarestorage.com
#aws s3 cp $DISTRIBUTE/bfvd_base_foldseekdb.tar.gz s3://bfvd --endpoint-url https://4dfc310765e5ee60c351741d16a92839.r2.cloudflarestorage.com
#aws s3 cp $DISTRIBUTE/bfvd_base_foldcompdb.tar.gz s3://bfvd --endpoint-url https://4dfc310765e5ee60c351741d16a92839.r2.cloudflarestorage.com
