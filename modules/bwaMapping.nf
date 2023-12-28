process bwaMapping {
    input:
    path reference_genome
    // index
    tuple path("${reference_genome}.amb"), path("${reference_genome}.ann"), path("${reference_genome}.bwt"), path("${reference_genome}.pac"), path("${reference_genome}.sa")
    // reads
    tuple val(sampleId), path(reads)

    output:
    tuple val(sampleId), path('mapped_reads.bam')

    script:
    """
    bwa mem -t `nproc` ${reference_genome} ${reads[0]} ${reads[1]} | samtools view -b - > mapped_reads.bam
    """
}