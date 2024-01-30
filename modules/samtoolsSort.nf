process samtoolsSort {
    input:
    tuple val(sampleId), path(bam_file)

    output:
    tuple val(sampleId), path('sorted_reads.bam')

    script:
    """
    samtools sort -@ ${params.threads} ${bam_file} -o sorted_reads.bam
    """
}