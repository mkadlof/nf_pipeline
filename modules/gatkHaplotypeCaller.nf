process gatkHaplotypeCaller {
    publishDir "results/${sampleId1}", mode: 'symlink'

    input:
    path reference_genome_fasta, name: 'reference_genome.fasta'
    path reference_genome_fai, name: 'reference_genome.fasta.fai'
    path reference_genome_dict, name: 'reference_genome.dict'
    tuple val(sampleId1), path(inputBam)
    tuple val(sampleId2), path(inputBamBai)

    output:
    tuple val(sampleId1), path('gatkHaplotypeCaller.vcf.gz')

    script:
    """
    ${params.gatkPath} HaplotypeCaller --java-options "-Xmx12g" \
     -R ${reference_genome_fasta} \
     -I ${inputBam} \
     -O gatkHaplotypeCaller.vcf.gz \
     -ploidy 1
    """
}