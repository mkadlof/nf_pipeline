process applyLowCovToRef {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    path('reference.fasta')
    tuple val(sampleId), path('low_coverage.bed')

    output:
    tuple val(sampleId), path("reference_masked.fasta")

    script:
    """
    bedtools maskfasta -fi reference.fasta -bed low_coverage.bed -fo reference_masked.fasta
    """
}