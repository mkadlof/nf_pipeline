process detectLowCoverage {
    publishDir "results/${sampleId1}", mode: 'symlink'

    input:
    tuple val(sampleId1), path('input.bam'), path('input.bam.bai')

    output:
    tuple val(sampleId1), path('low_coverage.bed')

    script:
    """
    bedtools genomecov -bga -ibam input.bam | awk '\$4 < ${params.coverage_threshold}' | bedtools merge > low_coverage.bed
    """
}