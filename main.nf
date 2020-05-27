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
        val(var)

    output:
        path("varset.tsv")

    """
    
    """

}

process generate_sites {

    """
    python bamsurgeon/addsnv.py --reference {reference} \
                        --tmpdir ${TMPDIR} \
                        --procs ${task.cpus} \
                        --maxdepth 2000 \
                        --mindepth 1 \
                        -m 1.0 \
                        -v ${varset} \
                        --bamfile {config[bam_location]} \
                        --aligner mem \
                        --picardjar /lscr2/andersenlab/dec211/variant-caller-analysis/tools/picard.jar \
                        --outbam {output.spiked_bam}
    """

}


process bamsurgeon {

    container "lethalfang/bamsurgeon"


}

workflow {

    Channel.from


}