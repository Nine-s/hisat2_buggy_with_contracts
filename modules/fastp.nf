process FASTP {
    label 'fastp'
    publishDir params.outdir

    input:
    //tuple val(name), path(reads_1), path(reads_2)
    tuple val(name), path(reads)

    output:
    tuple val(name), path("${name}*.trimmed.fastq"), emit: sample_trimmed
    path "${name}_fastp.json", emit: json_report

    require(["""#!/usr/bin/env python3
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
    promise([FOR_ALL("f", ITER("*_fastp.json"), {f -> IF_THEN(EMPTY_FILE(f), "exit 1")}), FOR_ALL("f", ITER("*.trimmed.fastq"), {f -> IF_THEN(EMPTY_FILE(f), "exit 1")}), """#!/usr/bin/env python3
import json
import sys
import os

file = [f for f in os.listdir() if f.endswith("_fastp.json")][0]
with open(file) as log:
        j_dict = json.load(log)
        summary = j_dict["summary"]
        reads_before = summary["before_filtering"]["total_reads"]
        reads_after = summary["after_filtering"]["total_reads"]
        sys.exit(1 if (reads_before - reads_after) / reads_before > 0.95 else 0)""", """#!/usr/bin/env python3
import sys
import os

counter = 0
length = 0
for file in os.listdir():
    if file.endswith(".trimmed.fastq"):
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
    sys.exit(1)""", COMMAND_LOGGED_NO_ERROR(), INPUTS_NOT_CHANGED()])

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${name}.R1.trimmed.fastq -O ${name}.R2.trimmed.fastq --thread ${task.cpus} --json ${name}_fastp.json --thread ${params.threads}
    """
}
