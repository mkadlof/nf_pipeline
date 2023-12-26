process samtoolsViewFilter {
    input:
    path inputBam

    output:
    path 'output_filtered.bam'

    script:
    """
    samtools view -b ${inputBam} -F 2820 -T 30 > output_filtered.bam
    """
}