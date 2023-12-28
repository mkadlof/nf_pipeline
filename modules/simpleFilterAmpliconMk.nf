process simpleFilterAmpliconMk {
    input:
    tuple val(sampleId1), path(inputBam)
    path(indexBai)

    output:
    tuple val(sampleId1), path('output_sorted_downsampled.bam')

    script:
    """
    simple_filter_amplicon_mk_illumina.py --cycles 30 --mode paired ${inputBam} output_sorted_downsampled.bam
    """
}