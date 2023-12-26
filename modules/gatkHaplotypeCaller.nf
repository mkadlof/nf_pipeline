process gatkHaplotypeCaller {
    input:
    path reference_genome_fasta, name: 'reference_genome.fasta'
    path reference_genome_fai, name: 'reference_genome.fasta.fai'
    path reference_genome_dict, name: 'reference_genome.dict'
    path inputBam
    output:
    path 'gatkHaplotypeCaller.vcf.gz'
    path 'gatkHaplotypeCaller.bam'

    script:
    """
    ${params.gatkPath} HaplotypeCaller --java-options "-Xmx12g" \
     -R ${reference_genome_fasta} \
     -I ${inputBam} \
     -O gatkHaplotypeCaller.vcf.gz \
     -bamout gatkHaplotypeCaller.bam \
     -ploidy 1
    """
}