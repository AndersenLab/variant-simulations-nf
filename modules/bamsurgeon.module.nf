/*
    Bamsurgeon spike-ins
*/

params.picard_path = "/opt/conda/share/picard-2.22.3-0/picard.jar"

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
        --tagreads \\
        --force \\
        -m 1.0 \\
        --bamfile in.bam \\
        --aligner mem \
        --picardjar ${params.picard_path} \
        --outbam ${varset}_${var_type}_${real_or_simulated}.bam

    makevcf.py

        samtools quickcheck ${varset}_${var_type}_${real_or_simulated}.bam
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
        --tagreads \\
        --force \\
        -m 1.0 \\
        --bamfile in.bam \\
        --aligner mem \
        --picardjar ${params.picard_path} \
        --outbam ${varset}_${var_type}_${real_or_simulated}.bam
        
        # Check that bam aligned properly
        samtools quickcheck ${varset}_${var_type}_${real_or_simulated}.bam
    """

}