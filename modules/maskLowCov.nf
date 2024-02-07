process maskLowCov {
    publishDir "results/${sampleId1}", mode: 'symlink'

    input:
    tuple val(sampleId1), path('reference_masked.fasta'),  path('input.fasta')

    output:
    tuple val(sampleId1), path('output.fasta')

    script:
    """
    cat reference_masked.fasta input.fasta > unaligned.fasta
    mafft --auto unaligned.fasta > aligned.fasta
    parseMaskedAlignedFasta.py aligned.fasta output.fasta
    """
}