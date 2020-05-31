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
    varset_size=${params.varset_size}                                      
    """
}

log_summary()


process gen_varset_real {

    tag { "${varset}:${real_or_simulated}:${var_type}" }

    //conda "bedtools=2.29.2 bcftools=1.10 bedops=2.4.39"
    publishDir "${params.output}/varsets/vcf", mode: 'copy', pattern: "*vcf*"

    input:
        tuple val(varset), \
              val(real_or_simulated), \
              val(var_type), \
              path(bam)

    output:
        tuple val(varset), \
              val(real_or_simulated), \
              val(var_type), \
              path("out.tsv"), emit: "to_process_varset"
        tuple path("${varset}_${var_type}.vcf.gz"), \
              path("${varset}_${var_type}.vcf.gz.csi")

    """
            # Generate simulated variant set
            sc sample -n ${params.varset_size} --types=${var_type} ${params.vcf_file} | \
            bcftools sort -O z > ${varset}_${var_type}.vcf.gz
            bcftools index ${varset}_${var_type}.vcf.gz
            
            # Convert to 
            if [[ "${var_type}" == "snps" ]]; then
                bcftools view ${varset}_${var_type}.vcf.gz | \
                bcftools query -f "%CHROM\t%POS\t%POS\t1.0\t%ALT{0}\t%ALT{0}\n" > out.tsv
            elif [[ "${var_type}" == "indels" ]]; then
                {
                    # Insertions ~ VCF -> BED; 0-based; But VCF starts indels at -1
                    bcftools view ${varset}_${var_type}.vcf.gz | \\
                    vcf2bed --insertions - | \\
                    cut -f 1,2,7 | \
                    awk '{ print \$1, \$2 + 1, \$2 + 2, "INS", substr(\$3, 2), \$3; }'; 
                    
                    # Deletions; Similar to above but 
                    bcftools view ${varset}_${var_type}.vcf.gz | \\
                    vcf2bed --deletions - | \\
                    cut -f 1,2,6 | \\
                    awk '{ print \$1, \$2 + 1, \$2 + length(substr(\$3, 2)), "DEL", substr(\$3, 2), \$3; }';
                } > out.tsv
            fi;
    """

}

process gen_varset_simulated {

    tag { "${varset}:${real_or_simulated}:${var_type}" }

    input:
        tuple val(varset), \
              val(real_or_simulated), \
              val(var_type), \
              path(bam)

    output:
        tuple val(varset), \
              val(real_or_simulated), \
              val(var_type), \
              path("out.tsv")

    """
        if [[ "${var_type}" == "snps" ]]; then

            # Gen 2x vars and filter for non-ref positions
            # This is an inefficient but easy way to do it (-;
            paste  <(sc rand -n ${params.varset_size*2} ${params.reference}) \\
                   <(sc rand -n ${params.varset_size*2} ${params.reference} | cut -f 4) | \
            awk '\$4 != \$5 { print \$1, \$2, \$3, "1.0", \$4, \$4 }' | \\
            head -n ${params.varset_size} > out.tsv

        elif [[ "${var_type}" == "indels" ]]; then
            # Insertions ~ Generate random seq from elsewhere in the genome to insert
            {
                paste  <(sc rand -n ${params.varset_size/2} --dist=2-50 ${params.reference} | cut -f 1,2,4) \\
                       <(sc rand -n ${params.varset_size/2} --dist=2-50 ${params.reference} | cut -f 4) | \\
                    awk '{ print \$1, \$2, \$2+1, "INS", \$4, substr(\$4,1), \$4 }';
                # Deletions
                sc rand -n ${params.varset_size/2} --dist=2-50 ${params.reference} | \
                awk '{ print \$1, \$2, \$3 + length(substr(\$4,2)) + 1, "DEL", substr(\$4,2), \$4; }';
            } | sort -k 1,1 -k 2,2n - > out.tsv
        fi;
    """
}


process process_varset {

    tag { "${varset}:${real_or_simulated}:${var_type}" }

    publishDir "${params.output}/varsets", mode: 'copy', pattern: "*.tsv"

    input:
        tuple val(varset), \
              val(real_or_simulated), \
              val(var_type), \
              path("varset.tsv")

    output:
        tuple val(varset), \
              val(real_or_simulated), \
              val(var_type), 
              path("${varset}_${var_type}_${real_or_simulated}.tsv"), emit: "to_bamsurgeon"
        path "varset_for_analysis.combine", emit: "to_combine"

    """
        sort -k 1,1 -k 2,2n varset.tsv | \\
        tr ' ' '\t' | \\
        awk -v OFS='\t' '{ print \$0, "${varset}", "${real_or_simulated}", "${var_type}" }' > ${varset}_${var_type}_${real_or_simulated}.tsv

        # Add a header to varset for downstream analysis
        {
            echo -e "CHROM\tPOS\tALT\tstart\tstop\tsnp_ins_del\tspike_sequence\tvarset\treal_or_simulated\tvar_type";
                awk -v OFS="\t" '{ if ("${var_type}" == "snps") { col4="SNP" } else { col4=\$4 };
                       print \$1, \$2, \$6, \$2, \$3, col4, \$5, "${varset}", "${real_or_simulated}", "${var_type}" }' varset.tsv
        } > varset_for_analysis.combine

    """
}


process bamsurgeon_spike_snps {

    container "erikwaskiewicz/bamsurgeon"

    input:
        tuple val(varset), \
              val(real_or_simulated), \
              val(var_type), 
              path("varfile.tsv")

    // output:
    //     path("${varset_n}.bam")

    """
    addsnv.py \\
        --reference ${reference} \\
        --varfile varfile.tsv \\
        --tmpdir . \\
        --procs ${task.cpus} \\
        --maxdepth 2000 \\
        --mindepth 1 \\
        -m 1.0 \\
        -v varset.tsv \\
        --bamfile in.bam \\
        --aligner mem \
        --picardjar /usr/local/bin/picard.jar \
        --outbam ${varset}_${var_type}_${real_or_simulated}.bam
    """

}


process bamsurgeon_spike_indels {

    container "erikwaskiewicz/bamsurgeon"

    input:
        tuple val(varset), \
              val(real_or_simulated), \
              val(var_type), 
              path("varfile.tsv")

    // output:
    //     path("${varset_n}.bam")

    """
    addindel.py \\
        --varfile varfile.tsv \\
        --reference ${reference} \\
        --tmpdir . \\
        --procs ${task.cpus} \\
        --maxdepth 2000 \\
        --mindepth 1 \\
        -m 1.0 \\
        -v varset.tsv \\
        --bamfile in.bam \\
        --aligner mem \
        --picardjar /usr/local/bin/picard.jar \
        --outbam ${varset}_${var_type}_${real_or_simulated}.bam
    """

}


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
    process_varset.out.to_bamsurgeon.branch {
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