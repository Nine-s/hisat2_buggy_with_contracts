process SAMTOOLS {
    label 'samtools'
    publishDir params.outdir
    
    input:
    tuple val(sample_name), path(sam_file)
    
    output:
    tuple val(sample_name), path("${sample_name}.sorted.bam"), emit: sample_bam 

    promise([FOR_ALL('f', ITER("*.sorted.bam"), { f -> IF_THEN(EMPTY_FILE(f), "exit 1")}), FOR_ALL("x", ITER("*.sam"), {RETURN(GREATER_THAN(NUM("\$(du \$(basename \$x).sam | sed 's/\t.*//g')"), NUM("\$(du \$(basename \$x).sorted.bam | sed 's/\t.*//g')")))}), COMMAND_LOGGED_NO_ERROR(), INPUTS_NOT_CHANGED()])
    
    script:
    """
    samtools view -bS ${sample_name}.sam -@ ${params.threads} | samtools sort -o ${sample_name}.sorted.bam -T tmp  -@ ${params.threads} 
    """
    
}
