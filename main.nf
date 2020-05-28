#!/usr/bin/env nextflow 
/*
    Andersen Lab Variant simulations pipeline
    Authors:
    - Daniel Cook <danielecook@gmail.com>
*/
nextflow.preview.dsl=2

// Params
reference = file(params.reference, checkIfExists: true)
bamfile = file(params.bamfile, checkIfExists: true)

def log_summary() {

    out =  '''
    Variant Simulations                                           
    '''
}

log_summary()


process gen_varset_real {

    tag { "${varset}:${simulation_type}:${var_type}"}

    conda "bedtools=2.29.2"

    input:
        tuple val(varset), \
              val(simulation_type), \
              val(var_type), \
              path(bam)

    output:
        tuple val(varset_n), path("${varset_n}.tsv")

    """
        bedtools
    """

}

process generate_snv_bams {

    conda "bamsurgeon="

    input:
        tuple val(varset_n), path("varset.tsv"), path("in.bam")

    output:
        path("${varset_n}.bam")

    """
    python bamsurgeon/addsnv.py \\
                        --reference ${reference} \
                        --tmpdir . \\
                        --procs ${task.cpus} \\
                        --maxdepth 2000 \\
                        --mindepth 1 \\
                        -m 1.0 \\
                        -v varset.tsv \\
                        --bamfile in.bam \\
                        --aligner mem \
                        --picardjar /lscr2/andersenlab/dec211/variant-caller-analysis/tools/picard.jar \
                        --outbam ${varset_n}.bam
    """
}



process generate_indel_bams {

    container "lethalfang/bamsurgeon"

    input:
        tuple val(varset_n), path("varset.tsv"), path("in.bam")

    output:
        path("${varset_n}.bam")

    """
    python bamsurgeon/addindel.py \\
                        --reference ${reference} \
                        --tmpdir . \\
                        --procs ${task.cpus} \\
                        --maxdepth 2000 \\
                        --mindepth 1 \\
                        -m 1.0 \\
                        -v varset.tsv \\
                        --bamfile in.bam \\
                        --aligner mem \
                        --picardjar /lscr2/andersenlab/dec211/variant-caller-analysis/tools/picard.jar \
                        --outbam ${varset_n}.bam
    """

}


workflow {

    varsets = Channel.from(1..10)

    varsets.combine(["real", "simulated"])
           .combine(["snp", "indel"])
           .combine([bamfile])
           .branch {
               real: it[1] == "real"
               simulated: it[1] == "simulated"
           }.set { varsets_in }
    
    varsets_in.real | gen_varset_real


}