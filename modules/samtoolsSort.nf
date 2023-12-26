process samtoolsSort {
    input:
    path bam_file

    output:
    path 'sorted_reads.bam'

    script:
    """
    samtools sort ${bam_file} -o sorted_reads.bam
    """
}