process NORMALIZE_VCF {
    container 'community.wave.seqera.io/library/bcftools_awscli:e8024de90248c3f6'
    publishDir "${params.outdir_base}/normalized_vcfs", mode: 'copy'
    memory '8.GB'
    cpus 2
    time '2.h'

    input:
    tuple val(patient_id), val(sample_meta), path(vcf_file), path(reference_fasta)

    output:
    tuple val(patient_id), val(sample_meta), path("norm_${sample_meta.sample_id}.vcf"), emit: normalized_vcfs

    script:
    """
    echo "Normalizing VCF for Patient ${patient_id}, Sample ${sample_meta.sample_id}..."
    echo "Using reference: ${reference_fasta} (size: \$(ls -lh ${reference_fasta} | awk '{print \$5}'))"
    
    # Check if reference file exists and is not empty
    if [[ ! -s "${reference_fasta}" ]]; then
        echo "ERROR: Reference file ${reference_fasta} is empty or missing"
        touch norm_${sample_meta.sample_id}.vcf
        exit 1
    fi
    
    # Check if VCF file exists and is not empty
    if [[ ! -s "${vcf_file}" ]]; then
        echo "Warning: VCF file ${vcf_file} is empty or missing. Creating empty normalized VCF."
        touch norm_${sample_meta.sample_id}.vcf
        exit 0
    fi
    
    # Normalize the VCF file using the reference FASTA
    echo "Starting bcftools normalization..."
    bcftools norm -f ${reference_fasta} -o norm_${sample_meta.sample_id}.vcf ${vcf_file}
    
    # Check if normalization was successful
    if [[ ! -f "norm_${sample_meta.sample_id}.vcf" ]]; then
        echo "Warning: Normalization failed. Creating empty VCF file."
        touch norm_${sample_meta.sample_id}.vcf
    else
        echo "Normalization successful. Output size: \$(ls -lh norm_${sample_meta.sample_id}.vcf | awk '{print \$5}')"
    fi
    
    echo "Normalization complete for ${sample_meta.sample_id}!"
    """
}