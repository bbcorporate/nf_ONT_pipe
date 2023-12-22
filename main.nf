
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

def buffer_size = 2

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
    debug true
    cpus 1
    memory '1 GB'
    //label = [ 'process_ultralow', 'error_retry' ]      
    
    input:
    path this_ref
    path this_gff3
    path gatk_vcf_out
    
    output:
    path "*.csq.ann.vcf"

    stub:
    println("In bcf_csq gatk_vcf_out is ${gatk_vcf_out}")
    """
    for vcf_name in ${gatk_vcf_out.join(' ')}; do
        echo STUB
        touch \$vcf_name.csq.ann.vcf
    done

"""

    script:
    """
for vcf_name in ${gatk_vcf_out.join(' ')}; do
    bcftools csq -p a  -f ${this_ref} -g ${this_gff3}  --verbose 2 -o \$vcf_name.csq.ann.vcf \$vcf_name
done
"""
}
process show_thing { 
     cpus 1
    memory '10 MB'
    debug true
    
    input:
    val p1

    script:
    println("thing is ${p1}")
    """
    date
    """
}
process show_six {
    cpus 1
    memory '10 MB'
    debug true
    
    input:
    val p1
    val p2
    val p3
    val p4
    val p5
    val p6

    script:
    println("p1 is ${p1}")
    println("p2 is ${p2}")
    println("p3 is ${p3}")
    println("p4 is ${p4}")
    println("p5 is ${p5}")
    println("p6 is ${p6}")
    """
    """


}
process show_this {
    cpus 1 
    memory '10 MB'
    input:
    val f1
    val f2
    val f3
    val f4

    script:
    println("f1 is $f1")
    println("f2 is $f2")
    println("f3 is $f3")
    println("f4 is $f4")
    """
    echo ${f1}
    """
}
process make_ngmlr_filenames {
    cpus 1
    memory '10 MB'
    input:
    path this_fq

    output:
    path this_fq
    val params.path_ref
    val params.enc
    val params.ht

    script:
    '''
    date
    ''' 
}
process make_bcftools_filenames {
    cpus 1
    memory '10 MB'
    input:
    path this_vcf
    output:
    val params.path_ref
    val params.path_gff3
    path this_vcf

    script:
    println("bcf filenames this_vcf are ${this_vcf}")
    """
    date 
    """
 }

process make_gatk_filenames {
    cpus 1
    memory '10 MB'

    input:
    path this_bam  // was path this_bam -- it's a list
    path this_bai

    output:
    path this_bam
    path this_bai // was added_bai  // can't add suffix to each member of a list like ${this_bam}.bai
    val params.intervals
    val params.path_ref
    val "${params.path_ref}.fai"
    val params.path_ref.replace('.fa','.dict') //seriously

/*
    exec:   // just run groovy code
    //added_bai = this_bam.collect(it + '.bai')
    //println(this_bam.getClass())
    //println("Working in ${workflow.workDir}")
    println("IN make_gatk_filename{}")
    //println("IN MAKE_GATK_FLNM: task workdir is ${task.workDir}")
    */
    
    script:
    added_bai = this_bam.collect {"${params.s3dir}/${it}.bai"}   // was .join(' ')
    """
    date
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
    //was -r ${this_ref} 
    // then ${params.path_ref}  // but container won't find the path, even though I can tab complete and then find it.
    // -r ../../../${params.ref}  // this works, but still need to handle the enc and ht files
    """
ngmlr -x ont -r ${this_ref} -q ${this_fq} -o ${this_fq}.ngmlr.rg.sam  --rg-id Orders_Q8D_1_Pool_1  --rg-sm sample1  --rg-lb library1  --rg-pl ONT -t 8
"""

}

process ngmlr_samtools_batch {
    cpus 9
    memory '6 GB'
    container '454262641088.dkr.ecr.us-west-1.amazonaws.com/nf_ont_pipe_ecr:binReadsComplete3'
    debug true
    //publishDir = "${params.s3dir}"
    //label = [ 'process_medium', 'error_retry' ]
    
    input:
    path this_fq_batch
    path this_ref
    path this_enc  // ngm index file
    path this_ht   // ngm index file
    
    output:
    path "*.rg.sort.bam"
    path "*.rg.sort.bam.bai"  // make this in process make_gatk_filenames
    
    script:
    //was -r ${this_ref} 
    // then ${params.path_ref}  // but container won't find the path, even though I can tab complete and then find it.
    // -r ../../../${params.ref}  // this works, but still need to handle the enc and ht files
    """
    for fq_name in ${this_fq_batch.join(' ')}; do
        #echo IN \$PWD NGMLR_BATCH: this fq file is \$fq_name
        ngmlr -x ont -q \$fq_name  -r ${this_ref} -o \$fq_name.ngmlr.rg.sam  --rg-id Orders_Q8D_1_Pool_1  --rg-sm sample1  --rg-lb library1  --rg-pl ONT -t 9
        samtools view -b \$fq_name.ngmlr.rg.sam > \$fq_name.ngmlr.rg.bam
        samtools sort \$fq_name.ngmlr.rg.bam > \$fq_name.rg.sort.bam
        samtools index \$fq_name.rg.sort.bam
        #echo wrote \$fq_name.rg.sort.bam and .bam.bai
    done
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
process haplotype_caller_batch  {
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
    path "*.vcf"
    //path path_ref

    stub:
    """
    for bam_name in ${sorted_bam.join(' ')}; do
        echo STUB
        touch \$bam_name.vcf
    done
"""
    script:
    
	"""
    for bam_name in ${sorted_bam.join(' ')}; do
        echo IN GATK_BATCH \$PWD : this bam file is \$bam_name
        gatk --java-options "-Xmx4g" HaplotypeCaller -I \$bam_name -O \$bam_name.vcf -R $path_ref  -ploidy ${params.ploidy} -L $intervals  --disable-read-filter WellformedReadFilter --do-not-run-physical-phasing true
    done
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
    // Channel.of(1..23)| buffer(size: buffer_size, remainder: true) | map { "this is the batch ${it}" } | foo
    // Channel.fromPath("testChunk*.fq") | ngmlr(...) | samtools_post_process
    //Channel.of(1..23)| buffer(size: buffer_size, remainder: true) | foo
    // this works to make buffers, but have to split the buffers further inside
    // the process zoo
    //Channel.fromPath("seq*") | buffer(size: buffer_size, remainder: true) | zoo

    
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
    def existing_cluster_fq = Channel.fromPath("${params.s3dir}/GANDER_1_cl=?_*reads.snip.fq")
    
    
    // this takes 7 input fq files, splits into buffers of 3 and adds the filenames needed to run ngmlr and prints out the filenames
    // the first filename is a list
    //existing_cluster_fq | buffer (size: buffer_size, remainder: true ) | make_ngmlr_filenames | show_this

    //make a second buffer for the sorted bam files from samtools post processing
    // we're making a list of lists
    //  [[/Users/.../work/1b/76543f09caf1240c438aebd7ecfbff/GANDER_1_cl=1_735reads.fq.rg.sort.bam, /Users/.../work/1b/76543f09caf1240c438aebd7ecfbff/GANDER_1_cl=4_298reads.fq.rg.sort.bam, 
    //  /Users/.../1b/76543f09caf1240c438aebd7ecfbff/GANDER_1_cl=8_2reads.fq.rg.sort.bam], /Users/simon... ]]
    // => I don't need the second buffer()
    //existing_cluster_fq | buffer (size: buffer_size, remainder: true ) | make_ngmlr_filenames | ngmlr_samtools_batch | buffer(size: buffer_size, remainder: true) | show_thing

    // this prints out six things 1,2 are lists, 3,4,5,6 are scalar things
    //existing_cluster_fq | buffer (size: buffer_size, remainder: true ) | make_ngmlr_filenames | ngmlr_samtools_batch | make_gatk_filenames | show_six

    //makes vcf files
    //existing_cluster_fq | buffer (size: buffer_size, remainder: true ) | make_ngmlr_filenames | ngmlr_samtools_batch | make_gatk_filenames | haplotype_caller_batch

    // goes to annotated vcf files
    existing_cluster_fq | buffer (size: buffer_size, remainder: true ) | make_ngmlr_filenames | ngmlr_samtools_batch | make_gatk_filenames | haplotype_caller_batch | make_bcftools_filenames | bcftools_csq
    
    // try simple connection to this
    //bcftools_csq( params.path_ref, params.path_gff3 , haplotype_caller_batch.out )


    //def vcf_files =Channel.fromPath(haplotype_caller_batch.out)
    // try something like this 
    // vcf_files | buffer (size: buffer_size, remainder: true) | 

    /*
    // took out params.path_ref, 
    // container can't find it, even though it's there!!!
    // need to refer to it as ../../../theNameOfThe.fa
    //ngmlr(existing_cluster_fq, params.path_ref, params.enc, params.ht)  | samtools_post_process
    haplotype_caller( samtools_post_process.out[0] , samtools_post_process.out[1], params.intervals, params.path_ref, index_reference.out, create_seq_dict.out )
    //haplotype_caller( samtools_post_process.out[0] , samtools_post_process.out[1], params.intervals, params.path_ref, index_reference.out, create_seq_dict.out  )
    bcftools_csq( params.path_ref, params.path_gff3 , haplotype_caller.out )

*/


    // below is testing just bcftools csq piece
    //bcftools_csq( params.path_ref, params.path_gff3 , "${launchDir}/gatk.out.vcf" )
    //    bcftools_csq( ${params.gatk_vcf_out} )
    // println("Using manually diploidified VCF file")
    // bcftools_csq( "${workflow.launchDir}/out.diploidified.unphased.vcf" )

}

