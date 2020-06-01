/*
    Bamsurgeon spike-ins
*/

params.picard_path = "/opt/conda/share/picard-2.22.3-0/picard.jar"

process bamsurgeon_spike_snps {

    publishDir "${params.output}/bam", mode: 'copy', pattern: "*bam*"
    publishDir "${params.output}/truth", mode: 'copy', pattern: "*vcf*"

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

    output:
        tuple path("${varset}_${var_type}_${real_or_simulated}.bam"), \
              path("${varset}_${var_type}_${real_or_simulated}.bam.bai")
        tuple path("${varset}_${var_type}_${real_or_simulated}.truth.vcf.gz"), \
              path("${varset}_${var_type}_${real_or_simulated}.truth.vcf.gz.csi")

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
        --outbam out.bam

    # Fix up the truth set
    bcftools reheader --fai ${params.reference}.fai out.addsnv.varfile.vcf | \\
    bcftools sort -O z > ${varset}_${var_type}_${real_or_simulated}.truth.vcf.gz
    bcftools index ${varset}_${var_type}_${real_or_simulated}.truth.vcf.gz

    samtools sort --threads ${task.cpus} -O BAM out.bam > ${varset}_${var_type}_${real_or_simulated}.bam
    samtools index -@ ${task.cpus} ${varset}_${var_type}_${real_or_simulated}.bam

    # Check that bam aligned properly
    samtools quickcheck ${varset}_${var_type}_${real_or_simulated}.bam
    """

}


process bamsurgeon_spike_indels {

    publishDir "${params.output}/bam", mode: 'copy', pattern: "*bam*"
    publishDir "${params.output}/truth", mode: 'copy', pattern: "*vcf*"

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

    output:
        tuple path("${varset}_${var_type}_${real_or_simulated}.bam"), \
              path("${varset}_${var_type}_${real_or_simulated}.bam.bai")
        tuple path("${varset}_${var_type}_${real_or_simulated}.truth.vcf.gz"), \
              path("${varset}_${var_type}_${real_or_simulated}.truth.vcf.gz.csi")

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
        --outbam out.bam

    # Fix up the truth set
    bcftools reheader --fai ${params.reference}.fai out.addindel.varfile.vcf | \\
    bcftools sort -O z > ${varset}_${var_type}_${real_or_simulated}.truth.vcf.gz
    bcftools index ${varset}_${var_type}_${real_or_simulated}.truth.vcf.gz

    samtools sort --threads ${task.cpus} -O BAM out.bam > ${varset}_${var_type}_${real_or_simulated}.bam
    samtools index -@ ${task.cpus} ${varset}_${var_type}_${real_or_simulated}.bam
    
    # Check that bam aligned properly
    samtools quickcheck ${varset}_${var_type}_${real_or_simulated}.bam
    samtools index ${varset}_${var_type}_${real_or_simulated}.bam
    """

}