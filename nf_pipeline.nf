
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

process bwaMapping {
    input:
    path reference_genome
    // index
    tuple path("${reference_genome}.amb"), path("${reference_genome}.ann"), path("${reference_genome}.bwt"), path("${reference_genome}.pac"), path("${reference_genome}.sa")
    // reads
    tuple val(sampleId), path(reads)

    output:
    path 'mapped_reads.bam'

    script:
    """
    bwa mem -t `nproc` ${reference_genome} ${reads[0]} ${reads[1]} | samtools view -Sb - > mapped_reads.bam
    """
}




workflow {
    reference_genome = Channel.fromPath(params.reference_genome)
    reads = Channel.fromFilePairs(params.reads)
    genomeStats(reference_genome)
    index = bwaIndex(reference_genome)
    bwaMapping(reference_genome, index, reads)
}
