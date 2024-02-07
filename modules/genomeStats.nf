process genomeStats {
    input:
    path reference_genome

    output:
    path 'genomeStats.json.gz'

    script:
    """
    genomeStats.py ${reference_genome}
    """
}