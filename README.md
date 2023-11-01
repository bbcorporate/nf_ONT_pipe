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
