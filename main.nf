#!/usr/bin/env nextflow 
/*
    Andersen Lab Variant simulations pipeline
    Authors:
    - Daniel Cook <danielecook@gmail.com>
*/
nextflow.preview.dsl=2

assert System.getenv("NXF_VER") == "20.04.1"

// Params
params.output = "output"
reference = file(params.reference, checkIfExists: true)
bam_file = file(params.bam_file, checkIfExists: true)

def log_summary() {

    out =  """
    Variant Simulations     

    bam_file=${bam_file}
    vcf_file=${params.vcf_file}
    reference=${reference}
    n_variants=${params.n_variants}                                      
    """
}

log_summary()

include { bamsurgeon_spike_snps;
          bamsurgeon_spike_indels; } from './modules/bamsurgeon.module.nf'
include { gen_varset_real;
          gen_varset_simulated;
          process_varset; } from './modules/varsets.module.nf'

workflow {

    varsets = Channel.from(1..params.n_varsets)

    // branch real vs. simulated
    varsets.combine(["real", "simulated"])
           .combine(["snps", "indels"])
           .combine([bam_file])
           .branch {
               real: it[1] == "real"
               simulated: it[1] == "simulated"
           }.set { varsets_in }
    
    // generate real vs simulated test sets
    varsets_in.real | gen_varset_real
    varsets_in.simulated | gen_varset_simulated

    gen_varset_real.out.to_process_varset
                   .mix(gen_varset_simulated.out) | process_varset

    // Branch snps and indels to bamsurgeon
    process_varset.out.to_bamsurgeon
                  .combine([file(params.bam_file, checkIfExists: true)])
                  .combine([file(params.bam_file + ".bai", checkIfExists: true)])
                  .branch {
        snps: it[2] == "snps"
        indels: it[2] == "indels"
    }.set { bamsurgeon_in }

    bamsurgeon_in.snps | bamsurgeon_spike_snps
    bamsurgeon_in.indels | bamsurgeon_spike_indels

    // Combine variant datasets for evaluation
    process_varset.out.to_combine.collectFile(name: "${params.output}/varset_summary.tsv",
                                             keepHeader: true,
                                             skip: 1)

}