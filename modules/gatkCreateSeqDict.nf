process gatkCreateSeqDict {
    input:
    path reference_genome_fasta, name: 'reference_genome.fasta'

    output:
    path "reference_genome.dict"

    script:
    """
    ${params.gatkPath} CreateSequenceDictionary \
     -R ${reference_genome_fasta}
    """
}