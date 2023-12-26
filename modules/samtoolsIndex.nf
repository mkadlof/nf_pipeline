process samtoolsIndex {
    input:
    path bam_file

    output:
    path "${bam_file}.bai"

    script:
    """
    samtools index ${bam_file}
    """
}