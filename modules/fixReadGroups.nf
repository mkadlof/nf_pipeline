process fixReadGroups {
    input:
    path bam
    output:
    path "${bam.baseName}.with_rg.bam"

    script:
    """
    java -jar ${params.picardPath} AddOrReplaceReadGroups \
     I=${bam} \
     O=${bam.baseName}.with_rg.bam \
     SORT_ORDER=coordinate \
     RGID=foo \
     RGLB=bar \
     RGPL=illumina \
     RGSM=Sample1 \
     RGPU=unit1 \
     CREATE_INDEX=True
    """
}