#!/usr/bin/env bash
for i in GANDER_1.310memLib_cl=?_reads.fq; do
    echo $i;
    cmd=`printf '{"command": ["bbcorporate/nf_ONT_pipe", "-qs 1", "-r c9821d7", "--ref DSB127_701_NoStop_UpdatedHisReference.fa", "--gff3 DSB127_701_NoStop_UpdatedHisReference.gff3", "--s3dir s3://bb-aws-genomics-workflows-store/nf-pipeline-testing/nf_ONT_pipe/GANDER_1.310memLib", "--fqfile %s", "--intervals DSB127_701_NoStop_UpdatedHisMap:4218-5090", "--ploidy 2"]}' $i`
    #echo $cmd
    aws --profile birch batch submit-job \
    --job-name nfONTpipe_job_testing \
    --job-queue priority-gwfcore \
    --job-definition nextflow-nextflow-resources \
    --container-overrides "$cmd"

    echo Submitted $i
    #exit 1

done;
