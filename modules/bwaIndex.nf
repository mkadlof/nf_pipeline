process bwaIndex {
    input:
    path reference_genome
    output:
    tuple path("${reference_genome}.amb"), path("${reference_genome}.ann"), path("${reference_genome}.bwt"), path("${reference_genome}.pac"), path("${reference_genome}.sa")

    script:
    """
    bwa index ${reference_genome}
    """
}