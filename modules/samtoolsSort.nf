process samtoolsSort {
    input:
    tuple val(sampleId), path(bam_file)

    output:
    tuple val(sampleId), path('sorted_reads.bam'), path('sorted_reads.bam.bai')

    script:
    """
    samtools sort -@ ${params.threads} ${bam_file} -o sorted_reads.bam
    samtools index sorted_reads.bam
    """
}