process Scispeak {
    
    publishDir "${params.outdir}/scispeak/${sample_id}", mode: 'symlink'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("scispeak_${sample_id}_R[12].fastq.gz"), emit: reads_ch
    path("scispeak_${sample_id}_log.json")

    script:
    """
    scispeak \
        -i ${reads[0]} \
        -I ${reads[1]} \
        -o scispeak_${sample_id} \
        -t ${params.n_threads} \
        --log;
    """

}
