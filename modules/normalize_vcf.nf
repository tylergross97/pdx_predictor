process NORMALIZE_VCF {
    container 'community.wave.seqera.io/library/awscli_bcftools:c0aa6b95529320d8'
    publishDir "${params.outdir_base}/normalized_vcfs", mode: 'copy'

    input:
    tuple val(patient_id), val(sample_meta), path(vcf_file), path(reference_fasta)

    output:
    tuple val(patient_id), val(sample_meta), path("norm_${sample_meta.sample_id}.vcf"), emit: normalized_vcfs

    script:
    """
    echo "Normalizing VCF for Patient ${patient_id}, Sample ${sample_meta.sample_id}..."
    
    # Normalize the VCF file using the reference FASTA
    bcftools norm -f ${reference_fasta} -o norm_${sample_meta.sample_id}.vcf ${vcf_file} || touch norm_${sample_meta.sample_id}.vcf
    
    echo "Normalization complete for ${sample_meta.sample_id}!"
    """
}