
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

params.extrabinparams = ''

params.path_ref = "${params.s3dir}/${params.ref}"
params.path_gff3 = "${params.s3dir}/${params.gff3}"
params.enc = "${params.path_ref}-enc.2.ngm"
params.ht  = "${params.path_ref}-ht-13-2.2.ngm"
//params.fai = "${params.path_ref}.fai"    // *** this adds .fai to fa filename
println("====================================")
println("=   This is the ONT-UMI pipeline   =")
println("====================================")


println("projectDir (main.script.nf in here) is: $projectDir")
println("launchDir (nf run ... in here) is: $launchDir")
println("params.path_ref is: $params.path_ref")
 

//for debugging
//System.exit(0)

// steps with my container
///////////////////////////////

process bin_reads_by_umi {
    debug true
    cpus 1
    memory '6 GB'
    //label = [ 'process_medium', 'error_retry' ]
    container '454262641088.dkr.ecr.us-west-1.amazonaws.com/nf_ont_pipe_ecr:binReadsComplete3'
    publishDir = "${params.s3dir}"
    println("header in process bin_reads")
    
    input:
    path this_fq
    path this_gb
    
    output:
    path "*cl=*reads.fq"

    script:
    """
bin_reads_by_umi.py -d ${params.depth} -gb ${params.gb}  -fq ${this_fq}  -expectedreadlength ${params.expreadlen} ${params.extrabinparams} -slop ${params.slop} 
"""
}

process bcftools_csq {
    container '454262641088.dkr.ecr.us-west-1.amazonaws.com/nf_ont_pipe_ecr:binReadsComplete3'
    publishDir = "${params.s3dir}"
    cpus 1
    memory '1 GB'
    //label = [ 'process_ultralow', 'error_retry' ]      
    
    input:
    path this_ref
    path this_gff3
    path gatk_vcf_out
    
    output:
    path "${gatk_vcf_out}.csq.ann.vcf"

    script:
    """
bcftools csq -p a  -f ${this_ref} -g ${this_gff3}  --verbose 2 -o ${gatk_vcf_out}.csq.ann.vcf ${gatk_vcf_out} 
"""
}

process ngmlr {
    cpus 8
    memory '6 GB'
    container '454262641088.dkr.ecr.us-west-1.amazonaws.com/nf_ont_pipe_ecr:binReadsComplete3'
    //label = [ 'process_medium', 'error_retry' ]
    
    input:
    path this_fq
    path this_ref
    path this_enc  // ngm index file
    path this_ht   // ngm index file
    
    output:
    path "${this_fq}.ngmlr.rg.sam"
    
    script:
    """
ngmlr -x ont -r ${this_ref} -q ${this_fq} -o ${this_fq}.ngmlr.rg.sam  --rg-id Orders_Q8D_1_Pool_1  --rg-sm sample1  --rg-lb library1  --rg-pl ONT -t 8
"""

}

process samtools_post_process {
    container '454262641088.dkr.ecr.us-west-1.amazonaws.com/nf_ont_pipe_ecr:binReadsComplete3'
    label = [ 'process_ultralow', 'error_retry' ]    
	
    input:
    path this_sam
    
    output:
    path "${this_sam}.rg.sort.bam", emit: sorted_bam_out
    path "${this_sam}.rg.sort.bam.bai", emit: sorted_bam_out_bai
    
    script:
    """
samtools view -b $this_sam > out.ngmlr.rg.bam
samtools sort out.ngmlr.rg.bam > ${this_sam}.rg.sort.bam
samtools index ${this_sam}.rg.sort.bam
"""
}


process index_reference {
    // outputs reference fa file 
    // AND
    // index .fai file
    debug true
    container '454262641088.dkr.ecr.us-west-1.amazonaws.com/nf_ont_pipe_ecr:binReadsComplete3'
    publishDir = "${params.s3dir}"
    label = [ 'process_ultralow', 'error_retry' ]
    //stageOutMode = 'copy'
    
    input:
    path this_ref
    
    output:
    //path this_ref
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
    //label = [ 'process_medium', 'error_retry' ]    
    cpus 1
    memory '1 GB'
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
    //label = [ 'process_medium', 'error_retry' ]
    cpus 1
    memory '1GB'

    input:
    path sorted_bam
    path sorted_bam_bai
    val intervals
    path path_ref
    path path_ref_index
    path path_ref_dict

    output:
    path "${sorted_bam}.vcf"
    //path path_ref

    script:
    
	"""
gatk --java-options "-Xmx4g" HaplotypeCaller -I $sorted_bam -O ${sorted_bam}.vcf -R $path_ref  -ploidy ${params.ploidy} -L $intervals  --disable-read-filter WellformedReadFilter --do-not-run-physical-phasing true  --seconds-between-progress-updates 10
"""
    
}

// testing processes

process foo {
    debug true
    label = [ 'process_ultralow', 'error_retry' ]
    
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
    label = [ 'process_ultralow', 'error_retry' ]
    
    input:
    path x
    output:
    path 'test*.txt'
    script:
	"""
for i in ${x}; do head -4 \$i | cut -c 1-10 > test.\$i.txt ; done  
"""
}
    
    process list_file{
        debug true
        input:
        path a_file
        output:
        path "out.txt"
        script:
        """
        echo HELLO
        ls -ltrh ${a_file} > out.txt
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

    
    // working on linking bin_reads
    def gb_path = params.s3dir + '/' + params.gb
    def fq_path = params.s3dir + '/' +  params.fqfile
    def fq_chann = Channel.fromPath(fq_path)
    println("s3dir is ${params.s3dir}")
    println("fqfile is ${params.fqfile}")
    println("fq_path is: $fq_path")
    println("depth is: $params.depth")
    
    //list_file(fq_chann)
    
/*   NOT USING THIS FOR NOW.
CONNECT OUTPUT BACK TO NGMLR WHEN BATCHES ARE FIXED UP
    bin_reads_by_umi(fq_chann, gb_path)
    //println("after bin_reads call BHGF")
*/

    index_reference(params.path_ref)
    create_seq_dict(params.path_ref)
    

    //def binned_fq_chann = Channel.fromPath(bin_reads_by_umi.out)
    
    //ngmlr(bin_reads_by_umi.out, params.path_ref, params.enc, params.ht)  | samtools_post_process
//CONNECT TO PRE-EXISTING OUTPUT FILES FROM BIN_UMI - SEE ABOVE
    def existing_cluster_fq = Channel.fromPath("$launchDir/GANDER_1_cl=?_*reads.fq")
    ngmlr(existing_cluster_fq, params.path_ref, params.enc, params.ht)  | samtools_post_process
    haplotype_caller( samtools_post_process.out[0] , samtools_post_process.out[1], params.intervals, params.path_ref, index_reference.out, create_seq_dict.out )
    //haplotype_caller( samtools_post_process.out[0] , samtools_post_process.out[1], index_reference.out[0], index_reference.out[1], create_seq_dict.out, params.intervals )
    bcftools_csq( params.path_ref, params.path_gff3 , haplotype_caller.out )





    // below is testing just bcftools csq piece
    //bcftools_csq( params.path_ref, params.path_gff3 , "${launchDir}/gatk.out.vcf" )
    //    bcftools_csq( ${params.gatk_vcf_out} )
    // println("Using manually diploidified VCF file")
    // bcftools_csq( "${workflow.launchDir}/out.diploidified.unphased.vcf" )

}

