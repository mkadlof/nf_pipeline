process iVar {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(bam), path(bai)
    path(primers)

    output:
    tuple val(sampleId), path("ivar_output.bam"), path('ivar_output.bam.bai')

    script:
    """
    ivar trim -i ${bam} \
              -b ${primers} \
              -e \
              -p ivar_output_tmp.bam
    samtools sort -@ ${params.threads} ivar_output_tmp.bam -o ivar_output.bam
    samtools index ivar_output.bam
    """
}