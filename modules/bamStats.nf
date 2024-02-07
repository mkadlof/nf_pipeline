process bamStats {
    input:
    tuple val(sampleId), path(inputBam), path(inputBamBai)

    output:
    tuple val(sampleId), path('bamStats.json.gz')

    script:
    """
    bamStats.py ${inputBam} --sample ${sampleId}
    """
}