process vcf2fasta {

    input:
    path reference_genome_fasta, name: 'reference_genome.fasta'
    tuple val(sampleId1), path(variants_vcf, name: 'variants.vcf.gz')
    tuple val(sampleId2), path(variants_vcf_tbi, name: 'variants.vcf.gz.tbi')

    output:
    tuple val(sampleId1), path("output.fasta")

    script:
    """
    bcftools consensus -f reference_genome.fasta variants.vcf.gz -o output-tmp.fasta
    awk '/^>/{print ">${sampleId1}"; next} 1' output-tmp.fasta > output.fasta
    """
}









