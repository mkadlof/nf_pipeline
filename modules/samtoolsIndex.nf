process samtoolsIndex {
    publishDir "results/${sampleId}", mode: 'symlink', pattern: 'output_filtered.bam.bai'

    input:
    tuple val(sampleId), path(bam_file)

    output:
    path "${bam_file}.bai"

    script:
    """
    samtools index ${bam_file}
    """
}