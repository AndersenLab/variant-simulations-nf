/*
    Varset generation and processing
*/

/*
    Varsets are output as:

    SNP ~ 1-based:
    chrom
    start
    end
    vaf
    allele (alt)

    INDEL ~ 0-based:
    chrom
    start
    end
    vaf
    INS/DEL
    allele

*/

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
            sc sample -n ${params.n_variants*2} --types=${var_type} ${params.vcf_file} | \
            bcftools filter --include '%FILTER="PASS"' | \
            bcftools sort -O z > ${varset}_${var_type}.vcf.gz
            bcftools index ${varset}_${var_type}.vcf.gz
            
            # Convert to 
            if [[ "${var_type}" == "snps" ]]; then
                bcftools view ${varset}_${var_type}.vcf.gz | \
                bcftools query -f "%CHROM\t%POS\t%POS\t1.0\t%ALT{0}\n" > out.tsv
            elif [[ "${var_type}" == "indels" ]]; then
                {
                    # Insertions ~ VCF -> BED; 0-based; But VCF starts indels at -1
                    bcftools view ${varset}_${var_type}.vcf.gz | \\
                    vcf2bed --insertions - | \\
                    cut -f 1,2,7 | \
                    awk '{ print \$1, \$2 + 1, \$2 + 2, "1.0", "INS", substr(\$3, 2); }'; 
                    
                    # Deletions; Similar to above but 
                    bcftools view ${varset}_${var_type}.vcf.gz | \\
                    vcf2bed --deletions - | \\
                    cut -f 1,2,6 | \\
                    awk '{ print \$1, \$2 + 1, \$2 + length(\$3), "1.0", "DEL", substr(\$3, 2); }';
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
            paste  <(sc rand -n ${params.n_variants*2} ${params.reference}) \\
                   <(sc rand -n ${params.n_variants*2} ${params.reference} | cut -f 4) | \
            awk '\$4 != \$5 { print \$1, \$2, \$3, "1.0", \$5 }' | \\
            head -n ${params.n_variants} > out.tsv

        elif [[ "${var_type}" == "indels" ]]; then
            # Insertions ~ Generate random seq from elsewhere in the genome to insert
            {
                paste  <(sc rand -n ${params.n_variants/2} --dist=2-30 ${params.reference} | cut -f 1,2,4) \\
                       <(sc rand -n ${params.n_variants/2} --dist=2-30 ${params.reference} | cut -f 4) | \\
                    awk '{ print \$1, \$2, \$2+1, "1.0", "INS", \$4 }';
                # Deletions
                sc rand -n ${params.n_variants/2} --dist=2-30 ${params.reference} | \
                awk '{ print \$1, \$2, \$2 + length(\$4), "1.0", "DEL", \$4; }';
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
            echo -e "CHROM\tPOS\tVAR\tstart\tstop\tsnp_ins_del\tspike_sequence\tvarset\treal_or_simulated\tvar_type";
                awk -v OFS="\t" '{ if ("${var_type}" == "snps") { vargroup="SNP"; allele=\$5 } else { vargroup=\$5; allele=\$6 };
                       print \$1, \$2, allele, \$2, \$3, \$4, vargroup, "${varset}", "${real_or_simulated}", "${var_type}" }' varset.tsv
        } > varset_for_analysis.combine

    """
}