
nextflow.enable.dsl = 2

include { FASTQC } from './modules/fastqc.nf'
include { FASTP } from './modules/fastp'
include { HISAT2_INDEX_REFERENCE ; HISAT2_ALIGN } from './modules/hisat2.nf'
include { SAMTOOLS } from './modules/samtools'
include { CUFFLINKS } from './modules/cufflinks'

log.info """\
         RNAseq differential analysis using NextFlow 
         =============================
         outdir: ${params.outdir}
         basedir : ${params.basedir}
         """
         .stripIndent()
 
params.outdir = 'results'


workflow {
    
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true ) 
    FASTQC( read_pairs_ch )
    FASTP( read_pairs_ch )
    HISAT2_INDEX_REFERENCE( params.reference_genome )
    HISAT2_ALIGN(FASTP.out.sample_trimmed, HISAT2_INDEX_REFERENCE.out)
    SAMTOOLS(HISAT2_ALIGN.out.sample_sam)
    CUFFLINKS(SAMTOOLS.out.sample_bam, params.reference_annotation)
}

