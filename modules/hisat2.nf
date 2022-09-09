process HISAT2_INDEX_REFERENCE {
    label 'hisat2'
    publishDir params.outdir
    
    input:
    path(reference)

    output:
    tuple path(reference), path("${reference.baseName}*.ht2")
    
    require(["""exit 0
import sys
import os

for file in [f for f in os.listdir() if f.endswith(".fa")]:
        with open(file) as f:
                for line in f:
                        if line[0] not in [">", "A", "C", "T", "G", "U", "N", ";"]:
                                sys.exit(1)"""])
    promise([COMMAND_LOGGED_NO_ERROR(), INPUTS_NOT_CHANGED()])

    script:
    """
    hisat2-build ${reference} ${reference.baseName} -p ${params.threads} 
    """
}

process HISAT2_ALIGN {
    label 'hisat2'
    publishDir params.outdir
    
    input:
    tuple val(sample_name), path(reads)
    tuple path(reference), path(index)

    output:
    path "${sample_name}_summary.log", emit: log
    tuple val(sample_name), path("${sample_name}.sam"), emit: sample_sam 

    promise([
        FOR_ALL("f", ITER("*.sam"), { f -> IF_THEN(OR(EMPTY_FILE(f), NOT(COND("head -n 1 " + f + " | grep -E '^@'"))), "exit 1")}), 
        FOR_ALL("f", ITER("*.log"), { f -> IF_THEN(EMPTY_FILE(f), "exit 1")}), 
        FOR_ALL("f", 
            ITER("*_summary.log"), 
            { f ->
                RETURN(
                    GREATER_THAN(
                        NUM("\$(grep Overall " + f + " | grep -Eo '[0-9]*\\.[0-9]*%' | sed 's/.[0-9][0-9]%//g')"), 
                        NUM(5)
                    )
                )
            }
        ), COMMAND_LOGGED_NO_ERROR(), INPUTS_NOT_CHANGED()
    ])

    script:
    """
    hisat2 -x ${reference.baseName} -1 ${reads[0]} -2 ${reads[1]} --new-summary --summary-file ${sample_name}_summary.log -S ${sample_name}.sam --thread ${params.threads} 
    sort ${sample_name}.sam > ${sample_name}.sam
    """
}
