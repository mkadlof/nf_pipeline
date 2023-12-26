process bamStats {
    input:
    path inputBam
    path indexBai
    output:
    path 'bamStats.json.gz'

    script:
    """
    bamStats.py ${inputBam}
    """
}