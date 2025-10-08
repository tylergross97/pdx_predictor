process STAGE_REFERENCE {
    container 'community.wave.seqera.io/library/bcftools_awscli:e8024de90248c3f6'
    storeDir "${params.outdir_base}/reference"  // Permanent cache
    memory '2.GB'
    cpus 1
    
    output:
    path "hg19.fa", emit: reference
    
    script:
    """
    echo "Downloading reference genome (this may take a few minutes)..."
    aws s3 cp ${params.reference_s3_path} hg19.fa
    echo "Reference genome download complete: \$(ls -lh hg19.fa)"
    """
}