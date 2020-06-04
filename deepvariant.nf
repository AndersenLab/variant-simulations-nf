/*
    Simple workflow to call variants using deep  variant
*/
nextflow.preview.dsl=2

assert System.getenv("NXF_VER") == "20.04.1"

process deepvariant {

    tag { "${fname}" }
    container "google/deepvariant:0.10.0"
    publishDir "callset/deepvariant", mode: 'copy'
    label 'xl'
    scratch false
    maxRetries 3
    errorStrategy 'retry'

    input:
        tuple val(fname), path(bam), path(bai)

    output:
        tuple path("${fname}.deepvariant.vcf.gz"), path("${fname}.deepvariant.vcf.gz.tbi")
        tuple path("${fname}.deepvariant.g.vcf.gz"), path("${fname}.deepvariant.g.vcf.gz.tbi"), emit: gvcf

    """
    /opt/deepvariant/bin/run_deepvariant \\
        --model_type=WGS \\
        --ref=${params.reference} \\
        --reads=${bam} \\
        --output_vcf=${fname}.deepvariant.vcf.gz \\
        --output_gvcf=${fname}.deepvariant.g.vcf.gz \\
        --num_shards=${task.cpus}
    """
}

workflow {
    bams = Channel.fromPath("output/bam/*.bam")
    bams.map { [it.getSimpleName(), it, it + ".bai"] } | deepvariant
}
