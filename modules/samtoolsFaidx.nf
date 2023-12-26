process samtoolsFaidx{
    input:
    path reference_genome
    output:
    path "${reference_genome}.fai"
    
    script:
    """
    samtools faidx ${reference_genome}
    """
}