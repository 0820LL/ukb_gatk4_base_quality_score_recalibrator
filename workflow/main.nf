// Declare syntax version
nextflow.enable.dsl=2

process GATK4_BASERECALIBRATOR {

    container = "${projectDir}/../singularity-images/depot.galaxyproject.org-singularity-gatk4-4.3.0.0--py36hdfd78af_0.img"

    input:
    path cram
    path fasta
    path bed_file
    path known_vcf
    path known_indels
    path standard_indels
    path fai
    path dict
    path tbi

    output:
    path "*.table"

    script:
    """
    gatk --java-options "-Xmx4g" BaseRecalibrator \\
        --input $cram \\
        --output ${params.prefix}.table \\
        --reference $fasta \\
        --intervals $bed_file \\
        --known-sites $known_vcf \\
        --known-sites $known_indels \\
        --known-sites $standard_indels \\
        --tmp-dir .
    cp ${params.prefix}.table ${launchDir}/${params.outdir}/
    """
}

workflow{
    bam             = Channel.of(params.cram)
    fasta           = Channel.of(params.fasta)
    bed_file        = Channel.of(params.bed_file)
    known_vcf       = Channel.of(params.known_vcf)
    known_indels    = Channel.of(params.known_indels)
    standard_indels = Channel.of(params.standard_indels)
    fasta_dict      = Channel.of(params.fasta_dict)
    fasta_fai       = Channel.of(params.fasta_fai)
    tbi             = Channel.of(params.known_vcf_tbi, params.known_indels_tbi, params.standard_indels_tbi).collect()
    GATK4_BASERECALIBRATOR(bam, fasta, bed_file, known_vcf, known_indels, standard_indels, fasta_fai, fasta_dict, tbi)
}

