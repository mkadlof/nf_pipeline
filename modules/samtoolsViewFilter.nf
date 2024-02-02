process samtoolsViewFilter {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(inputBam), path(inputBamBai)

    output:
    tuple val(sampleId), path('output_filtered.bam'), path('output_filtered.bam.bai')

    script:
    """
    samtools view -@ ${params.threads} -b ${inputBam} -F 2820 > output_filtered.bam
    samtools index output_filtered.bam
    """
}