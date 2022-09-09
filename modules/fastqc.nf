process FASTQC {
    label 'fastqc'
    tag "fastqc $sample_id"
    publishDir params.outdir

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*_fastqc.zip", emit: zip)

    require(["""exit 0
import sys
import os

counter = 0
length = 0
for file in os.listdir():
    if file.endswith(".fastq"):
        with open(file) as f:
            for line in f:
                if counter == 0:
                    if line[0] != "@" or line[1].isspace():
                        sys.exit(1)
                elif counter == 1:
#                   if any(map(lambda x: x not in ["G", "A", "T", "C", "N"], line[:-1])):
#                       sys.exit(1)
                    length = len(line)
                elif counter == 2:
                    if line[0] != "+":
                        sys.exit(1)
                elif counter == 3:
                    if len(line) != length:
                        sys.exit(1)
                counter = (counter + 1) % 4
if counter != 0:
    sys.exit(1)"""])
    promise([FOR_ALL("f", ITER("*_fastqc.zip"), {f -> IF_THEN(EMPTY_FILE(f), "exit 1")}), 
        CONTRACT(FOR_ALL("f", ITER("./*.zip"), { zipfile -> "unzip " + zipfile + "; " + IF_THEN(EQUAL(COUNT_PATTERN("\$(basename " + zipfile + ")" + "/summary.txt", "FAIL"), NUM(0)), "exit 1")}), "debug"), COMMAND_LOGGED_NO_ERROR(), INPUTS_NOT_CHANGED()])

    script:
    """
    fastqc ${reads[0]} ${reads[1]} --thread ${params.threads}
    """
}
