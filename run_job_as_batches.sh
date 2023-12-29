#!/usr/bin/env bash
this_release=${1:-v0.1}
cont_pieces=`printf '{"command": ["bbcorporate/nf_ONT_pipe", "-r %s", "--ref INS00182_DSB127_701_NoStop_UpdatedHisMapEngTag.fa", "--gff3 INS00182_DSB127_701_NoStop_UpdatedHisMapEngTag.gff3", "--gb INS00182_DSB127_701_NoStop_UpdatedHisMapEngTag.gb", "--s3dir s3://bb-aws-genomics-workflows-store/nf-pipeline-testing/nf_ONT_pipe/snip.310memberLibrary", "--fqfile snip.250k.reads.fq", "--intervals DSB127_701_NoStop_UpdatedHisMap:4218-5090", "--ploidy 1", "--expreadlen 960", "--slop 1", "--depth 10"]}' $this_release`

echo $cont_pieces

#exit 0

    aws --profile birch batch submit-job \
    --job-name nfONTpipe_job_testing \
    --job-queue priority-gwfcore \
    --job-definition nextflow-nextflow-resources \
    --container-overrides 5'{"command": ["bbcorporate/nf_ONT_pipe", "-r Latest", "--ref INS00182_DSB127_701_NoStop_UpdatedHisMapEngTag.fa", "--gff3 INS00182_DSB127_701_NoStop_UpdatedHisMapEngTag.gff3", "--gb INS00182_DSB127_701_NoStop_UpdatedHisMapEngTag.gb", "--s3dir s3://bb-aws-genomics-workflows-store/nf-pipeline-testing/nf_ONT_pipe/snip.310memberLibrary", "--fqfile snip.250k.reads.fq", "--intervals DSB127_701_NoStop_UpdatedHisMap:4218-5090", "--ploidy 1", "--expreadlen 960", "--slop 1", "--depth 10"]}'

    #echo Submitted $i
    
