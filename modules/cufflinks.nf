process CUFFLINKS {
    label 'cufflinks'
    publishDir params.outdir
     
    input:
    tuple val(sample_name), path(sorted_bam)
    path(annotation)
    
    output:
    path('transcripts.gtf'), emit: cufflinks_gtf 

    promise([RETURN(NOT(EMPTY_FILE('transcripts.gtf'))), "if [ \$(grep -Eo 'FPKM \"[^\"]*\"' transcripts.gtf | sort | uniq -c | wc -l) -eq 1 ]; then exit 1; fi", COMMAND_LOGGED_NO_ERROR(), INPUTS_NOT_CHANGED()])
    
    script:
    """
    cufflinks -G ${annotation} ${sorted_bam}  --num-threads ${params.threads}
    """
}
