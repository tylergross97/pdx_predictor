process LINE_BASED_DIVERGENCE {
    container 'community.wave.seqera.io/library/pip_awscli:8f8767791189d33e'
    publishDir "${params.outdir_base}/genomic_divergence", mode: 'copy'

    input:
    tuple val(patient_id), val(originator_meta), path(originator_vcf), val(passage_meta), path(passage_vcf)

    output:
    tuple val(patient_id), path("${patient_id}_${passage_meta.passage}_divergence_metrics.csv"), emit: divergence_metrics

    script:
    """
    echo "Processing divergence analysis for Patient ${patient_id}: ${originator_meta.sample_id} vs ${passage_meta.sample_id}"
    
    # Extract variant lines (skip headers) and sort for comparison
    echo "Extracting variants from originator VCF..."
    grep -v '^#' ${originator_vcf} | cut -f1-5 | sort > orig_variants.txt || touch orig_variants.txt
    
    echo "Extracting variants from passage VCF..."
    grep -v '^#' ${passage_vcf} | cut -f1-5 | sort > pass_variants.txt || touch pass_variants.txt
    
    # Count total variants
    ORIG_COUNT=\$(wc -l < orig_variants.txt)
    PASS_COUNT=\$(wc -l < pass_variants.txt)
    
    echo "Originator variants: \$ORIG_COUNT"
    echo "Passage variants: \$PASS_COUNT"
    
    # Find shared variants (intersection)
    echo "Computing variant intersections..."
    comm -12 orig_variants.txt pass_variants.txt > shared_variants.txt
    SHARED_COUNT=\$(wc -l < shared_variants.txt)
    
    # Find unique variants
    comm -23 orig_variants.txt pass_variants.txt > orig_only.txt
    comm -13 orig_variants.txt pass_variants.txt > pass_only.txt
    
    ORIG_ONLY=\$(wc -l < orig_only.txt)
    PASS_ONLY=\$(wc -l < pass_only.txt)
    
    echo "Shared variants: \$SHARED_COUNT"
    echo "Originator-only variants: \$ORIG_ONLY"
    echo "Passage-only variants: \$PASS_ONLY"
    
    # Calculate metrics and create CSV using gawk
    awk -v orig=\$ORIG_COUNT -v pass=\$PASS_COUNT -v shared=\$SHARED_COUNT \\
        -v orig_only=\$ORIG_ONLY -v pass_only=\$PASS_ONLY \\
        -v patient_id="${patient_id}" \\
        -v orig_sample="${originator_meta.sample_id}" \\
        -v pass_sample="${passage_meta.sample_id}" \\
        -v passage="${passage_meta.passage}" \\
        -v orig_specimen="${originator_meta.specimen_id}" \\
        -v pass_specimen="${passage_meta.specimen_id}" \\
        -v onco_code="${passage_meta.onco_tree_code}" \\
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
            printf "%s,%s,%s,%s,%s,%s,%s,%d,%d,%d,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\\n", 
                   patient_id, orig_sample, pass_sample, passage, orig_specimen, pass_specimen, onco_code,
                   orig, pass, shared, orig_only, pass_only, total_union, 
                   jaccard_sim, jaccard_dist, gain_rate, loss_rate, rel_div, burden_change
        }' > ${patient_id}_${passage_meta.passage}_divergence_metrics.csv
    
    echo "Divergence analysis complete!"
    echo "Jaccard distance: \$(awk 'NR==2 {print \$15}' FS=',' ${patient_id}_${passage_meta.passage}_divergence_metrics.csv)"
    
    # Clean up intermediate files
    rm -f orig_variants.txt pass_variants.txt shared_variants.txt orig_only.txt pass_only.txt
    """
}