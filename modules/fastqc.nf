process FASTQC {
    label 'fastqc'
    tag "fastqc $sample_id"
    publishDir params.outdir

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*_fastqc.zip", emit: zip)

    require([
	FOR_ALL("f", ITER("*.fastq"), 
		{ f -> 
			IF_THEN(
				NOT(
					EQUAL(
						NUM("\$(( \$(wc -l $f | cut -d' ' -f1)/4*4 ))"), 
						NUM("\$(wc -l $f | cut -d' ' -f1)")
					)
				), 
				"exit 1"
			)
		}
	) 
])
    promise([FOR_ALL("f", ITER("*_fastqc.zip"), {f -> IF_THEN(EMPTY_FILE(f), "exit 1")}), 
        CONTRACT(FOR_ALL("f", ITER("./*.zip"), { zipfile -> "unzip " + zipfile + "; " + IF_THEN(EQUAL(COUNT_PATTERN("\$(basename " + zipfile + ")" + "/summary.txt", "FAIL"), NUM(0)), "exit 1")}), "debug"), COMMAND_LOGGED_NO_ERROR(), INPUTS_NOT_CHANGED()])

    script:
    """
    fastqc ${reads[0]} ${reads[1]} --thread ${params.threads}
    """
}
