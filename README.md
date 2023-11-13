# nf_ONT_pipe
NF docker AWS pipe for ONT reads at BB


### Notes for launching a test run
**Oct 19 2023**  
Test pipeline with 

```
nf nf_containerAddGATKMyself.groovy --intervals pBB0212:5208-6083 --ref pBB0212.withCDSsf.noVERSION.resetOri.fa --ploidy 7 --snpeffdb rotated --gatk_vcf out.vcf
```

Need these input files in `.` (from plasmidsaurus seq run EXP00074/Orders_Q8D_raw_reads)
```
Orders_Q8D_1_Pool_1.fastq
Orders_Q8D_2_Pool_2.fastq
Orders_Q8D_3_Pool_3.fastq
Orders_Q8D_4_Pool_4.fastq
```

These files are in `s3://bb-aws-genomics-workflows-store/rawdata/nanopore_data/misc/`. 

The plasmid reference fasta file that is suitable for testing is in this dir `s3://bb-aws-genomics-workflows-store/nf-pipeline-testing/`

Also need this gff3 annotation file (replacing gb file format preferred by snpeff)  
`pBB0212.withCDSsf.noVERSION.resetOri.testing.gff3`  



# Running the pipeline, test files

S3 dir is `s3://bb-aws-genomics-workflows-store/nf-pipeline-testing/nf_ONT_pipe/`

Test input files are already set up there.

NF doesn't understand AWS profiles, or $AWS_PROFILE. You need to set up these ENV VARS  in the shell you will launch the pipeline from:  `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY`

You might have to login to ECR with somethign like this. Might not need the `--profile` part if you have the ENV VARS set up.

```
aws --profile XXXXX ecr get-login-password --region us-west-1 | docker login --username AWS --password-stdin 454262641088.dkr.ecr.us-west-1.amazonaws.com`
```

Images are either public (biocontainers/quay.io) or in company ECR

```
REPOSITORY                                                     TAG                               IMAGE ID       CREATED         SIZE
454262641088.dkr.ecr.us-west-1.amazonaws.com/nf_ont_pipe_ecr   ngmlrSamBcftools                  39be7c433321   10 days ago     2.65GB
quay.io/biocontainers/gatk4                                    4.4.0.0--py36hdfd78af_0           b22efbf7d870   8 months ago    1.06GB

```



run this

```
nf run -qs 1 main.nf --ref pBB0212.withCDSsf.noVERSION.resetOri.fa --gff3 pBB0212.withCDSsf.noVERSION.resetOri.testing.gff3 --s3dir s3://bb-aws-genomics-workflows-store/nf-pipeline-testing/nf_ONT_pipe --intervals pBB0212:5208-6083 --ploidy 2

```

