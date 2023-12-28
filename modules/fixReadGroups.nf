process fixReadGroups {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(bam)

    output:
    tuple val(sampleId), path("${bam.baseName}.rg.bam")
    tuple val(sampleId), path("${bam.baseName}.rg.bai")

    script:
    """
    java -jar ${params.picardPath} AddOrReplaceReadGroups \
     I=${bam} \
     O=${bam.baseName}.rg.bam \
     SORT_ORDER=coordinate \
     RGID=foo \
     RGLB=bar \
     RGPL=illumina \
     RGSM=${sampleId} \
     RGPU=unit1 \
     CREATE_INDEX=True
    """
}