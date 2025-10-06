process MERGE_FEATURES {
    container 'community.wave.seqera.io/library/awscli_pandas_python:da17666c6f28b039'
    publishDir "${params.outdir_base}/originator_vcf_features", mode: 'copy'

    input:
    val metas
    path feature_files

    output:
    path "all_originator_features.csv", emit: consolidated_features

    script:
    """
    echo "Merging all feature files into a single CSV..."

    # Create a simple metadata mapping file
    cat > metadata_mapping.txt << 'EOL'
${metas.withIndex().collect { meta, idx -> "${feature_files[idx]}\t${meta.patient_id}\t${meta.sample_id}" }.join('\n')}
EOL

    python3 <<EOF
import pandas as pd

# Initialize an empty DataFrame
merged_df = pd.DataFrame()

# Read metadata mapping
with open('metadata_mapping.txt', 'r') as f:
    for line in f:
        line = line.strip()
        if line:
            parts = line.split('\t')
            feature_file = parts[0]
            patient_id = parts[1]
            sample_id = parts[2]
            
            # Read the input feature file
            df = pd.read_csv(feature_file)
            
            # Add patient metadata as columns
            df['Patient_ID'] = patient_id
            df['Sample_ID'] = sample_id
            
            # Append to the merged DataFrame
            merged_df = pd.concat([merged_df, df], ignore_index=True)

# Save the consolidated file
merged_df.to_csv("all_originator_features.csv", index=False)
EOF

    echo "Consolidated features saved to all_originator_features.csv!"
    """
}