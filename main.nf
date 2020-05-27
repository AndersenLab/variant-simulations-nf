#!/usr/bin/env nextflow 
/*
    Andersen Lab Variant simulations pipeline
    Authors:
    - Daniel Cook <danielecook@gmail.com>
*/
nextflow.preview.dsl=2
params.reference = "_assets/genomes/c_elegans/PRJNA13758/WS276/c_elegans.PRJNA13758/WS276.genome.fa.gz"
params.bamfile = "_assets/bam/N2.bam"

// Params
reference = file(params.reference, checkIfExists: true)
bamfile = file(params.bamfile, checkIfExists: true)

def log_summary() {

    out =  '''
    Variant Simulations                                           
    '''

}

log_summary()


process generate_varsets {

    input:
        val(varset)

    output:
        tuple val(varset_n), path("${varset_n}.tsv")

    """
    
    """

}

process generate_snv_bams {

    container "lethalfang/bamsurgeon"

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


}