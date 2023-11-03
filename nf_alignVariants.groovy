
// BB ONT processing NF pipeline script

// this version:-
// 	works with a container from docker hub instead of local
// 	incorporates gatk4/createsequencedictionary image from nf-core
//      and snpeff image from biocontainers
//      I'm running my own process {} for these steps


//   In this script, these vars are defined implicitly
// launchDir   The directory where the workflow is run (requires version 20.04.0 or later).
// projectDir  The directory where the main script is located (previously baseDir)




// include modules we need
//include { GATK4_CREATESEQUENCEDICTIONARY } from './modules/nf-core/gatk4/createsequencedictionary/main'
//include { GATK4_HAPLOTYPECALLER } from './modules/nf-core/gatk4/haplotypecaller/main'  


// use this as a path, so needs to be absolute
// params.ref = "/Users/simonp/serenity_projects/birchbio/testingNF+Docker2/pBB0212.withCDSsf.noVERSION.resetOri.fa"
// replace with implicit variable
params.ref_fasta = "${projectDir}/${params.ref}"
params.enc = "${params.ref_fasta}-enc.2.ngm"
params.ht  = "${params.ref_fasta}-ht-13-2.2.ngm"
params.fai = "${params.ref_fasta}.fai"    // *** this adds .fai to fa filename

println("This is the  Birch ONT-UMI pipeline")
println("===================================")

println("projectDir is: $projectDir")
println("params.ref_fasta is: $params.ref_fasta")

//exit for debugging
//System.exit(0)

process simple_head {
	container 'procho/ont-pipe:ngmlrOverflowFixSamtoolsCmdLine'
	input:
	path this_file

	output:
	path "${this_file}.head.txt"
	
	script:
"""
head ${this_file} > ${this_file}.head.txt
"""
}


process build_snpeff_db {
    publishDir = "$projectDir"
//    container 'procho/ont-pipe:ngmlrOverflowFixSamtoolsCmdLine'
    input:
    path gb_file

    output:
    path "snpeff_db"
    path "snpeff.nfpipe.config"
    script:
    """
buildSnpEffDB.py -gb ${gb_file} -config snpeff.nfpipe.config -version ${params.snpeffdb}  -desc '${params.desc}' -data snpeff_db/data -continueanyway
"""
}


process snpeff {
    container 'quay.io/biocontainers/snpeff:5.1--hdfd78af_2'
    publishDir = "$projectDir"
    input:
    path gatk_vcf
    path snpeff_config
    output:
    path '*.ann.vcf', emit: ann_vcf_out
    path '*.stats.csv'
    script:
    // I made a sym link in project dir to snpeff dir in /opt/snpEff
    """
snpEff ann -dataDir snpeff_db/data -c ${snpeff_config}   ${params.snpeffdb}  ${gatk_vcf}  > gatk.ann.vcf
"""
}
    
process ngmlr {
    //stageInMode 'copy'
    container 'procho/ont-pipe:ngmlrOverflowFixSamtoolsCmdLine'
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
    container 'procho/ont-pipe:ngmlrOverflowFixSamtoolsCmdLine'
	
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

process create_seq_dict {
    // need to build container with java and gatk4
    // use the gatk4 container, just not the nf-core modules, which
    // are so slow it's not even a joke
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'

    publishDir = "$projectDir"
    input:
    path this_ref
    output:
    path '*.dict', emit: dict
    script:
	"""
gatk CreateSequenceDictionary -R $this_ref 
"""
}


process index_reference {
    publishDir = "$projectDir"
    input:
    path this_ref
    output:
    path "${this_ref}.fai"
    script:
    """
samtools faidx $this_ref
"""
}

process haplotype_caller {
    container 'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0'
    publishDir = "$projectDir"
    input:
    path sorted_bam
    path sorted_bam_bai
    val intervals
    path ref_fasta
    path ref_fasta_index
    path ref_fasta_dict
    output:
    path '*.vcf'

    script:
	"""
gatk HaplotypeCaller -I $sorted_bam -O gatk.out.vcf -R $ref_fasta -L $intervals -ploidy ${params.ploidy} --do-not-run-physical-phasing true --disable-read-filter WellformedReadFilter
"""
}

workflow {
    def fq_files = Channel.fromPath('Orders_Q8D_[1]*.fastq')

    //prepare snpeff db
    
    //build_snpeff_db( "$projectDir/${params.gb}")


    /*
    index_reference(params.ref_fasta)
    create_seq_dict(params.ref_fasta)
    ngmlr(fq_files,params.ref_fasta, params.enc, params.ht)  | samtools_post_process
    //something like this
    haplotype_caller( samtools_post_process.out[0] , samtools_post_process.out[1], params.intervals, params.ref_fasta, index_reference.out, create_seq_dict.out ) | snpeff
     */
    
    // snpeff
    snpeff( "${projectDir}/out.vcf" , build_snpeff_db.out[1])
}
