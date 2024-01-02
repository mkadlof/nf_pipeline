process vcf2fasta {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    path reference_genome_fasta, name: 'reference_genome.fasta'
    tuple val(sampleId), path(variants_vcf, name: 'variants.vcf.gz')
    path variants_vcf_tbi, name: 'variants.vcf.gz.tbi'

    output:
    path "output.fasta"

    script:
    """
    bcftools consensus -f reference_genome.fasta variants.vcf.gz -o output.fasta
    """
}









