process simpleFilterAmpliconMk {
    input:
    tuple val(sampleId), path(inputBam), path(indexBai)

    output:
    tuple val(sampleId), path('output_sorted_downsampled.bam'), path('output_sorted_downsampled.bam.bai')

    script:
    """
    simple_filter_amplicon_mk_illumina.py --cycles ${params.cycles} --mode single ${inputBam} output_sorted_downsampled.bam
    """
}