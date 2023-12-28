process bamStats {
    input:
    tuple val(sampleID), path(inputBam)
    path indexBai

    output:
    tuple val(sampleID), path('bamStats.json.gz')

    script:
    """
    bamStats.py ${inputBam} --sample ${sampleID}
    """
}