
// BB ONT processing NF pipeline script

// this version:-
// 	works with a container from docker hub instead of local
// 	incorporates gatk4/createsequencedictionary image from nf-core
//      and snpeff image from biocontainers
//      I'm running my own process {} for these steps
//      Moving from snpeff to bcftools csq
//      needs input gff3 not gb


//   In this script, these vars are defined implicitly
// launchDir   The directory where the workflow is run (requires version 20.04.0 or later).
// projectDir  The directory where the main script is located (previously baseDir)




// include modules we need
//include { GATK4_CREATESEQUENCEDICTIONARY } from './modules/nf-core/gatk4/createsequencedictionary/main'
//include { GATK4_HAPLOTYPECALLER } from './modules/nf-core/gatk4/haplotypecaller/main'  





// use this as a path, so needs to be absolute
// params.ref = "/Users/simonp/serenity_projects/birchbio/testingNF+Docker2/pBB0212.withCDSsf.noVERSION.resetOri.fa"
// replace with implicit variable
// launchDir is the directory where the workflow is run
// projectDir is the directory where the main script is located
// process {publishDir...} is the directory where the output is published

// pass in --s3dir as base dir for S3 files

/*
// trying to remove trailing slash so S3 doesn't make a nested empty named dir
// doesn't work though
if ( params.s3dir.endsWith('/') ) {
    params.abc = 'boo'
    println("params.s3dir ends with /")
    println("up to the end is ${params.s3dir.substring(0, params.s3dir.length() -1 )}")
    params.s3dir = params.s3dir.substring(0, params.s3dir.length() -1 )
}
println("params.s3dir is: ${params.s3dir}")

println("params.abc is: ${params.abc}")
*/

params.path_ref = "${params.s3dir}/${params.ref}"
params.path_gff3 = "${params.s3dir}/${params.gff3}"
params.enc = "${params.path_ref}-enc.2.ngm"
params.ht  = "${params.path_ref}-ht-13-2.2.ngm"
//params.fai = "${params.path_ref}.fai"    // *** this adds .fai to fa filename

println("This is the ONT-UMI pipeline")
println("===================================")

/*
println("projectDir (main.script.nf in here) is: $projectDir")
println("launchDir (nf run ... in here) is: $launchDir")
println("params.path_ref is: $params.path_ref")
 */

//for debugging
//System.exit(0)

// steps with my container
///////////////////////////////

process bcftools_csq {
    container '454262641088.dkr.ecr.us-west-1.amazonaws.com/nf_ont_pipe_ecr:ngmlrSamBcftools'
    publishDir = "${params.s3dir}"
    
    input:
    path this_ref
    path this_gff3
    path gatk_vcf_out
    
    output:
    path '*.csq.ann.vcf'

    script:
    """
bcftools csq -p a  -f ${this_ref} -g ${this_gff3}  --verbose 2 -o out.csq.ann.vcf ${gatk_vcf_out} 
"""
}

process ngmlr {
    container '454262641088.dkr.ecr.us-west-1.amazonaws.com/nf_ont_pipe_ecr:ngmlrSamBcftools'
    
    input:
    path this_fq
    path this_ref
    path this_enc  // ngm index file
    path this_ht   // ngm index file
    
    output:
    path "out.ngmlr.rg.sam"
    
    script:
    """
ngmlr -x ont -r ${this_ref} -q ${this_fq} -o out.ngmlr.rg.sam  --rg-id Orders_Q8D_1_Pool_1  --rg-sm sample1  --rg-lb library1  --rg-pl ONT -t 8
"""

}

process samtools_post_process {
    container '454262641088.dkr.ecr.us-west-1.amazonaws.com/nf_ont_pipe_ecr:ngmlrSamBcftools'
	
    input:
    path this_sam
    
    output:
    path "out.ngmlr.rg.sort.bam", emit: sorted_bam_out
    path "out.ngmlr.rg.sort.bam.bai", emit: sorted_bam_out_bai
    
    script:
    """
samtools view -b $this_sam > out.ngmlr.rg.bam
samtools sort out.ngmlr.rg.bam > out.ngmlr.rg.sort.bam
samtools index out.ngmlr.rg.sort.bam
"""
}

process index_reference {
    debug true
    container '454262641088.dkr.ecr.us-west-1.amazonaws.com/nf_ont_pipe_ecr:ngmlrSamBcftools'
    publishDir = "${params.s3dir}"
    //stageOutMode = 'copy'
    
    input:
    path this_ref
    
    output:
    path "${this_ref}.fai"
    script:
    """
samtools faidx $this_ref
"""
}


//GATK steps with biocontainer
////////////////////////////////

process create_seq_dict {
    // need to build container with java and gatk4
    // use the gatk4 container, just not the nf-core modules, which
    // are so slow it's not even a joke
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    publishDir = "$launchDir"
    input:
    path this_ref
    output:
    path '*.dict', emit: dict
    script:
	"""
gatk CreateSequenceDictionary -R $this_ref 
"""
}

process haplotype_caller {
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'
    publishDir = "$launchDir"
    input:
    path sorted_bam
    path sorted_bam_bai
    val intervals
    path path_ref
    path path_ref_index
    path path_ref_dict
    output:
    path 'gatk.*.vcf'

    script:
    
	""" date && 
gatk HaplotypeCaller -I $sorted_bam -O gatk.out.vcf -R $path_ref  -ploidy ${params.ploidy} -L $intervals  --disable-read-filter WellformedReadFilter --do-not-run-physical-phasing true && date
"""
    
}

// testing processes

process foo {
    debug true
    
    input:
    val x   // not value!!!
    //tuple x


    //println("data is: ${x}")
    
    script:
    """
echo "echo data is: '${x}'"
"""
}

process zoo {
    debug true
    publishDir = "$launchDir"
    input:
    path x
    output:
    path 'test*.txt'
    script:
	"""
for i in ${x}; do head -4 \$i | cut -c 1-10 > test.\$i.txt ; done  
"""
}
    
// main workflow
//////////////////////

workflow {
    // testing batches
    // this works
    // Channel.of(1..23).map { "this_is_batch ${it}" }.view()
    // Channel.of(1..23)| buffer(size: 4, remainder: true) | map { "this is the batch ${it}" } | foo
    // Channel.fromPath("testChunk*.fq") | ngmlr(...) | samtools_post_process
    //Channel.of(1..23)| buffer(size: 4, remainder: true) | foo
    // this works to make buffers, but have to split the buffers further inside
    // the process zoo
    //Channel.fromPath("seq*") | buffer(size: 4, remainder: true) | zoo

    
    /* THIS WAS THE PREVIOUS WORKFLOW */

    def fq_glob = params.s3dir + '/Orders_Q8D_*_1.fastq'
    println("fq_glob is: $fq_glob")
    def fq_files = Channel.fromPath( fq_glob )
    //println("fq_files is: $fq_files")
    //System.exit(0)

    // can't do this
    //def fasta_ref_channel = Channel.fromPath(params.ref)
    //index_reference(fasta_ref_channel)
    //System.exit(0)

    index_reference(params.path_ref)
    create_seq_dict(params.path_ref)

    //System.exit(0)

    
    ngmlr(fq_files,params.path_ref, params.enc, params.ht)  | samtools_post_process

/*
//    println("** Overwriting intervals **")
//    params.intervals = "pBB0212:5400-5450"
*/
    haplotype_caller( samtools_post_process.out[0] , samtools_post_process.out[1], params.intervals, params.path_ref, index_reference.out, create_seq_dict.out )
    bcftools_csq( params.path_ref, params.path_gff3 , haplotype_caller.out )

    // below is testing just bcftools csq piece
    //bcftools_csq( params.path_ref, params.path_gff3 , "${launchDir}/gatk.out.vcf" )
    //    bcftools_csq( ${params.gatk_vcf_out} )
    // println("Using manually diploidified VCF file")
    // bcftools_csq( "${workflow.launchDir}/out.diploidified.unphased.vcf" )

}

