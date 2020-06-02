#!/usr/bin/env nextflow 
/*
    Andersen Lab Variant simulations pipeline
    Authors:
    - Daniel Cook <danielecook@gmail.com>
*/
nextflow.preview.dsl=2

assert System.getenv("NXF_VER") == "20.04.1"
assert params.n_variants % 2 == 0 // must spike an even number of variants

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

include { gen_varset_real;
          gen_varset_simulated;
          process_varset; 
          mix_varsets;
          resample_varset; } from './modules/varsets.module.nf'
include { bamsurgeon_spike_snps;
          bamsurgeon_spike_indels; 
          process_truth_sets; } from './modules/bamsurgeon.module.nf'

workflow {

    real_or_simulated = Channel.from(["real", "simulated"])

    // branch real vs. simulated
    real_or_simulated.combine(["snps", "indels"])
                     .combine(1..params.n_varsets)
                     .combine([bam_file])
                     .branch {
               real: it[0] == "real"
               simulated: it[0] == "simulated"
           }.set { varsets_in }
    
    // generate real vs simulated test sets
    varsets_in.real | gen_varset_real
    varsets_in.simulated | gen_varset_simulated

    // Group by generation type and variant type and shuffle
    gen_varset_real.out.to_mix_varset
                   .concat(gen_varset_simulated.out.to_mix_varset).groupTuple(by: [0,1]) | mix_varsets

    // Resample variants
    gen_varset_real.out.to_resample_varset
                   .concat(gen_varset_simulated.out.to_resample_varset)
                   .combine(mix_varsets.out, by: [0,1]).view() | resample_varset | process_varset

    // Branch snps and indels to bamsurgeon
    process_varset.out.to_bamsurgeon
                  .combine([file(params.bam_file, checkIfExists: true)])
                  .combine([file(params.bam_file + ".bai", checkIfExists: true)])
                  .branch {
        snps: it[1] == "snps"
        indels: it[1] == "indels"
    }.set { bamsurgeon_in }

    bamsurgeon_in.snps | bamsurgeon_spike_snps
    bamsurgeon_in.indels | bamsurgeon_spike_indels

    // Combine variant datasets for evaluation
    process_varset.out.to_combine.collectFile(name: "${params.output}/varset_summary.tsv",
                                              keepHeader: true,
                                              skip: 1)

    // Combine Truth Sets
    bamsurgeon_spike_snps.out.truth_vcf
                         .concat(bamsurgeon_spike_indels.out.truth_vcf) | process_truth_sets
                         
    process_truth_sets.out.collectFile(name: "${params.output}/truth_summary.tsv",
                                      keepHeader: true,
                                      skip: 1)

}