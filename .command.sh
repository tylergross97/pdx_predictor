#!/bin/bash -ue
# Update package list and install required tools
apt-get update -qq
apt-get install -y -qq procps gawk coreutils grep

echo "Processing divergence analysis for Patient 184742: ORIGINATOR vs XEY"

# Extract variant lines (skip headers) and sort for comparison
echo "Extracting variants from originator VCF..."
grep -v '^#' 184742~138-R~ORIGINATOR~v2.0.4.51.0~WES.vcf | cut -f1-5 | sort > orig_variants.txt || touch orig_variants.txt

echo "Extracting variants from passage VCF..."
grep -v '^#' 184742~138-R~XEY~v2.0.4.51.0~WES.vcf | cut -f1-5 | sort > pass_variants.txt || touch pass_variants.txt

# Count total variants
ORIG_COUNT=$(wc -l < orig_variants.txt)
PASS_COUNT=$(wc -l < pass_variants.txt)

echo "Originator variants: $ORIG_COUNT"
echo "Passage variants: $PASS_COUNT"

# Find shared variants (intersection)
echo "Computing variant intersections..."
comm -12 orig_variants.txt pass_variants.txt > shared_variants.txt
SHARED_COUNT=$(wc -l < shared_variants.txt)

# Find unique variants
comm -23 orig_variants.txt pass_variants.txt > orig_only.txt
comm -13 orig_variants.txt pass_variants.txt > pass_only.txt

ORIG_ONLY=$(wc -l < orig_only.txt)
PASS_ONLY=$(wc -l < pass_only.txt)

echo "Shared variants: $SHARED_COUNT"
echo "Originator-only variants: $ORIG_ONLY"
echo "Passage-only variants: $PASS_ONLY"

# Calculate metrics and create CSV using gawk
gawk -v orig=$ORIG_COUNT -v pass=$PASS_COUNT -v shared=$SHARED_COUNT \
    -v orig_only=$ORIG_ONLY -v pass_only=$PASS_ONLY \
    -v patient_id="184742" \
    -v orig_sample="ORIGINATOR" \
    -v pass_sample="XEY" \
    -v passage="P0" \
    -v orig_specimen="138-R" \
    -v pass_specimen="138-R" \
    -v onco_code="CCRCC" \
    'BEGIN {
        total_union = orig + pass - shared
        jaccard_sim = (total_union > 0) ? shared / total_union : 0
        jaccard_dist = 1 - jaccard_sim
        gain_rate = (orig > 0) ? pass_only / orig : 0
        loss_rate = (orig > 0) ? orig_only / orig : 0
        rel_div = (total_union > 0) ? (orig_only + pass_only) / total_union : 0
        burden_change = (orig > 0) ? (pass - orig) / orig : 0

        # Print header
        print "patient_id,originator_sample,passage_sample,passage_number,originator_specimen,passage_specimen,onco_tree_code,total_originator_variants,total_passage_variants,shared_variants,originator_only_variants,passage_only_variants,total_union_variants,jaccard_similarity,jaccard_distance,variant_gain_rate,variant_loss_rate,relative_divergence,mutation_burden_change"

        # Print data
        printf "%s,%s,%s,%s,%s,%s,%s,%d,%d,%d,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n", 
               patient_id, orig_sample, pass_sample, passage, orig_specimen, pass_specimen, onco_code,
               orig, pass, shared, orig_only, pass_only, total_union, 
               jaccard_sim, jaccard_dist, gain_rate, loss_rate, rel_div, burden_change
    }' > 184742_P0_divergence_metrics.csv

echo "Divergence analysis complete!"
echo "Jaccard distance: $(gawk 'NR==2 {print $15}' FS=',' 184742_P0_divergence_metrics.csv)"

# Clean up intermediate files
rm -f orig_variants.txt pass_variants.txt shared_variants.txt orig_only.txt pass_only.txt
