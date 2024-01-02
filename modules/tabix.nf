process tabix {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(variants_vcf, name: 'variants.vcf.gz')

    output:
    path "variants.vcf.gz.tbi"

    script:
    """
    tabix -p vcf variants.vcf.gz
    """
}









