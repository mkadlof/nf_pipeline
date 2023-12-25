
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
    bwa mem -t `nproc` ${reference_genome} ${reads[0]} ${reads[1]} | samtools view -b - > mapped_reads.bam
    """
}

process samtoolsViewFilter {
    input:
    path inputBam

    output:
    path 'output_filtered.bam'

    script:
    """
    samtools view -b ${inputBam} -F 2820 -T 30 > output_filtered.bam
    """
}

process samtoolsSort {
    input:
    path bam_file

    output:
    path 'sorted_reads.bam'

    script:
    """
    samtools sort ${bam_file} -o sorted_reads.bam
    """
}

process samtoolsIndex {
    input:
    path bam_file

    output:
    path "${bam_file}.bai"

    script:
    """
    samtools index ${bam_file}
    """
}

process simpleFilterAmpliconMk {
    input:
    path inputBam
    path index

    output:
    path 'output_sorted_downsampled.bam'

    script:
    """
    simple_filter_amplicon_mk_illumina.py --cycles 30 --mode paired ${inputBam} output_sorted_downsampled.bam
    """
}


workflow {
    reference_genome = Channel.fromPath(params.reference_genome)
    reads = Channel.fromFilePairs(params.reads)
    genomeStats(reference_genome)
    index = bwaIndex(reference_genome)
    bwaMapping(reference_genome, index, reads)
    samtoolsViewFilter(bwaMapping.out)
    samtoolsSort(samtoolsViewFilter.out)
    samtoolsIndex(samtoolsSort.out)
    simpleFilterAmpliconMk(samtoolsSort.out, samtoolsIndex.out)
}
