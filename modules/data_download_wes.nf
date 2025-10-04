 process DATA_DOWNLOAD_WES {
    container 'python:3.9-slim'
    publishDir "${params.outdir_base}/vcfs", mode: 'copy'

    input:
    path samplesheet

    output:
    path "*.vcf", emit: vcf_files

    script:
    """
    # Install required Python packages
    pip install pandas requests

    # Python script to process the samplesheet and download VCF files
    python - <<EOF
import pandas as pd
import requests
from urllib.parse import urljoin

# Load the samplesheet
samplesheet = pd.read_csv("${samplesheet}")

# Debug: Print the cleaned samplesheet
print("Cleaned samplesheet:")
print(samplesheet)

# Get a list of Patient IDs with an ORIGINATOR sample
originator_patients = samplesheet[samplesheet['Sample ID'] == 'ORIGINATOR']['Patient ID'].unique()
print("Originator Patient IDs:")
print(originator_patients)

# Filter rows where Patient ID is in the list and PDM Type is PDX
filtered_rows = samplesheet[
    (samplesheet['Patient ID'].isin(originator_patients)) &
    (samplesheet['PDM Type'].str.strip().str.lower() == 'pdx')
]

# Debug: Print rows matching both filters
print("Rows matching Patient ID and PDM Type filter:")
print(filtered_rows)

# Extract VCF URLs and prepend the base URL
base_url = "https://pdmdb.cancer.gov"
vcf_urls = filtered_rows['VCF'].str.strip().apply(lambda x: urljoin(base_url, x)).tolist()

# Debug: Print VCF URLs
print("VCF URLs with Base URL:")
print(vcf_urls)

# Download each VCF file
for url in vcf_urls:
    if url:
        print(f"Downloading: {url}")
        response = requests.get(url)
        filename = url.split("/")[-1]
        with open(filename, "wb") as f:
            f.write(response.content)
EOF
    """
}