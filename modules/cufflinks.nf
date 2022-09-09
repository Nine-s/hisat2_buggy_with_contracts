process CUFFLINKS {
    label 'cufflinks'
    publishDir params.outdir
     
    input:
    tuple val(sample_name), path(sorted_bam)
    path(annotation)
    
    output:
    path('transcripts.gtf'), emit: cufflinks_gtf 

    promise([RETURN(NOT(EMPTY_FILE('transcripts.gtf'))), """exit 0
import sys

bag = dict()
with open("transcripts.gtf") as t:
    for line in t:
        line = line.split()
        if "FPKM" not in line:
            continue
        index = line.index("FPKM") + 1
        if line[index] not in bag:
            bag[line[index]] = 0
        bag[line[index]] += 1
if len(bag) < 2:
    sys.exit(1)""", COMMAND_LOGGED_NO_ERROR(), INPUTS_NOT_CHANGED()])
    
    script:
    """
    cufflinks -G ${annotation} ${sorted_bam}  --num-threads ${params.threads}
    """
}
