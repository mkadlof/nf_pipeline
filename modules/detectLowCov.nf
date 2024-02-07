process detectLowCov {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path('input.bam'), path('input.bam.bai')

    output:
    tuple val(sampleId), path('low_coverage.bed')

    script:
    """
    bedtools genomecov -bga -ibam input.bam | awk '\$4 < ${params.coverage_threshold}' | bedtools merge > low_coverage.bed
    """
}