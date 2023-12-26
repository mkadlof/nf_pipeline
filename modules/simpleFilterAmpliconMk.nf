process simpleFilterAmpliconMk {
    input:
    path inputBam
    path index

    output:
    path 'output_sorted_downsampled.bam'

    script:
    """
    simple_filter_amplicon_mk_illumina.py --cycles 30 --mode paired ${inputBam} output_sorted_downsampled.bam
    """
}