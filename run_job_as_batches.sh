#!/usr/bin/env bash

    #cmd=`printf '{"command": ["bbcorporate/nf_ONT_pipe", "-r d973094", "--ref DSB127_701_NoStop_UpdatedHisReference.fa", "--gff3 DSB127_701_NoStop_UpdatedHisReference.gff3", "--s3dir s3://bb-aws-genomics-workflows-store/nf-pipeline-testing/nf_ONT_pipe/snip.310memberLibrary", "--fqfile %s", "--intervals DSB127_701_NoStop_UpdatedHisMap:4218-5090", "--ploidy 1"]}' $i`
    #echo $cmd
    aws --profile birch batch submit-job \
    --job-name nfONTpipe_job_testing \
    --job-queue priority-gwfcore \
    --job-definition nextflow-nextflow-resources \
    --container-overrides '{"command": ["bbcorporate/nf_ONT_pipe", "-r d973094", "--ref INS00182_DSB127_701_NoStop_UpdatedHisMapEngTag.fa", "--gff3 INS00182_DSB127_701_NoStop_UpdatedHisMapEngTag.gff3", "--gb INS00182_DSB127_701_NoStop_UpdatedHisMapEngTag.gb", "--s3dir s3://bb-aws-genomics-workflows-store/nf-pipeline-testing/nf_ONT_pipe/snip.310memberLibrary", "--fqfile snip.250k.reads.fq", "--intervals DSB127_701_NoStop_UpdatedHisMap:4218-5090", "--ploidy 1", "--expreadlen 960", "--slop 1"]}

    #echo Submitted $i
    
