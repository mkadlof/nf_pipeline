process combineJsonFiles_x4 {
    publishDir "results/${sampleId1}", mode: 'symlink'

    input:
    // from genomeStats
    path(json1, name: '1.json.gz')
    // from bamStats
    tuple val(sampleId1), path(json2, name: '2.json.gz')
    tuple val(sampleId2), path(json3, name: '3.json.gz')
    tuple val(sampleId3), path(json4, name: '4.json.gz')

    output:
    path 'stats.json.gz'

    script:
    """
    combineJsonFiles.py --output stats.json.gz ${json1} ${json2} ${json3} ${json4}
    """
}