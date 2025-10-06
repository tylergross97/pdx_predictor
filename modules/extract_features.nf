process EXTRACT_FEATURES {
    container 'community.wave.seqera.io/library/awscli_pandas_python:da17666c6f28b039'
    publishDir "${params.outdir_base}/originator_vcf_features", mode: 'copy'

    input:
    tuple val(meta), path(vcf_file)

    output:
    tuple val(meta), path("features_${meta.patient_id}_${meta.sample_id}.csv"), emit: vcf_features

    script:
    """
    echo "Extracting features for ORIGINATOR VCF: Patient ${meta.patient_id}, Sample ${meta.sample_id}..."
    
    # Use pandas to parse the VCF file
    python3 <<EOF
import pandas as pd

# Read the VCF file, skipping header lines
with open("${vcf_file}") as f:
    lines = [line for line in f if not line.startswith('#')]

# Parse the VCF file into a DataFrame
columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
data = [line.strip().split('\\t')[:8] for line in lines]
df = pd.DataFrame(data, columns=columns)

# Save the selected columns to a CSV file
df[['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER']].to_csv("features_${meta.patient_id}_${meta.sample_id}.csv", index=False)
EOF

    echo "Feature extraction complete for ${meta.patient_id}_${meta.sample_id}!"
    """
}