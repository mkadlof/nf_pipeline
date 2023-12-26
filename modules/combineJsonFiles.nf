process combineJsonFiles_x2 {
    input:
    path json1, name: '1.json.gz'
    path json2, name: '2.json.gz'
    output:
    path 'stats.json.gz'

    script:
    """
    combineJsonFiles.py --output stats.json.gz ${json1} ${json2}
    """
}

process combineJsonFiles_x3 {
    input:
    path json1, name: '1.json.gz'
    path json2, name: '2.json.gz'
    path json3, name: '3.json.gz'
    output:
    path 'stats.json.gz'

    script:
    """
    combineJsonFiles.py --output stats.json.gz ${json1} ${json2} ${json3}
    """
}

process combineJsonFiles_x4 {
    input:
    path json1, name: '1.json.gz'
    path json2, name: '2.json.gz'
    path json3, name: '3.json.gz'
    path json4, name: '4.json.gz'
    output:
    path 'stats.json.gz'

    script:
    """
    combineJsonFiles.py --output stats.json.gz ${json1} ${json2} ${json3} ${json4}
    """
}