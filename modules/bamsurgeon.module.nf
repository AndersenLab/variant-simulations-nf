/*
    Bamsurgeon spike-ins
*/

process bamsurgeon_spike_snps {

    tag { "${varset}:${real_or_simulated}:${var_type}" }

    container "andersenlab/bamsurgeon"

    label 'md'

    input:
        tuple val(varset), \
              val(real_or_simulated), \
              val(var_type), \
              path("varfile.tsv"), \
              path("in.bam"), \
              path("in.bam.bai")

    // output:
    //     path("${varset_n}.bam")

    """
    addsnv.py \\
        --reference ${params.reference} \\
        --varfile varfile.tsv \\
        --tmpdir . \\
        --procs ${task.cpus} \\
        --maxdepth 2000 \\
        --mindepth 1 \\
        -m 1.0 \\
        --bamfile in.bam \\
        --aligner mem \
        --picardjar /usr/local/bin/picard.jar \
        --outbam ${varset}_${var_type}_${real_or_simulated}.bam
    """

}


process bamsurgeon_spike_indels {

    tag { "${varset}:${real_or_simulated}:${var_type}" }
    
    container "andersenlab/bamsurgeon"

    label 'md'

    input:
        tuple val(varset), \
              val(real_or_simulated), \
              val(var_type), \
              path("varfile.tsv"), \
              path("in.bam"), \
              path("in.bam.bai")

    // output:
    //     path("${varset_n}.bam")

    """
    addindel.py \\
        --varfile varfile.tsv \\
        --reference ${params.reference} \\
        --tmpdir . \\
        --procs ${task.cpus} \\
        --maxdepth 2000 \\
        --mindepth 1 \\
        -m 1.0 \\
        --bamfile in.bam \\
        --aligner mem \
        --picardjar /usr/local/bin/picard.jar \
        --outbam ${varset}_${var_type}_${real_or_simulated}.bam
    """

}