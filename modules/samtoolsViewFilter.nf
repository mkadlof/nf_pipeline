process samtoolsViewFilter {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(inputBam)

    output:
    tuple val(sampleId), path('output_filtered.bam')

    script:
    """
    samtools view -@ ${params.threads} -b ${inputBam} -F 2820 -T 30 > output_filtered.bam
    """
}