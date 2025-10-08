include { LINE_BASED_DIVERGENCE } from '../modules/line_based_divergence'
include { NORMALIZE_VCF } from '../modules/normalize_vcf'

workflow LINE_BASED_DIVERGENCE_ANALYSIS {
    take:
    vcf_with_reference_ch  // Channel: [meta, vcf_file, reference_fasta]

    main:
    // Prepare input for VCF normalization
    normalize_input = vcf_with_reference_ch
        .map { meta, vcf, reference ->
            [meta.patient_id, meta, vcf, reference]
        }

    // Run VCF normalization - now each process gets the same cached reference
    NORMALIZE_VCF(normalize_input)

    // Transform normalized output back to the expected format for divergence analysis
    normalized_vcfs = NORMALIZE_VCF.out.normalized_vcfs
        .map { patient_id, sample_meta, normalized_vcf ->
            [sample_meta, normalized_vcf]
        }

    // Group samples by patient_id and create originator-passage pairs
    grouped_samples = normalized_vcfs
        .map { meta, vcf ->
            [meta.patient_id, meta, vcf]
        }
        .groupTuple(by: 0)
        .map { patient_id, metas, vcfs ->
            // Find originator sample
            def originator_idx = metas.findIndexOf { it.passage == 'o' }
            def passage_samples = []
            
            // Collect all passage samples
            metas.eachWithIndex { meta, idx ->
                if (meta.passage != 'o') {
                    passage_samples.add([meta, vcfs[idx]])
                }
            }
            
            // Create comparisons if we have both originator and passage samples
            if (originator_idx >= 0 && passage_samples.size() > 0) {
                def originator_meta = metas[originator_idx]
                def originator_vcf = vcfs[originator_idx]
                
                // Create a comparison for each passage sample
                return passage_samples.collect { passage_meta, passage_vcf ->
                    [patient_id, originator_meta, originator_vcf, passage_meta, passage_vcf]
                }
            } else {
                log.warn "Patient ${patient_id}: Missing originator (${originator_idx >= 0 ? 'found' : 'not found'}) or passage samples (${passage_samples.size()} found)"
                return []
            }
        }
        .flatten()
        .collate(5)  // Group into tuples of 5 elements
        .filter { it.size() == 5 }  // Only keep complete tuples

    // Debug: Show what comparisons will be made
    grouped_samples.view { patient_id, orig_meta, orig_vcf, pass_meta, pass_vcf ->
        "Creating comparison: Patient ${patient_id} - ${orig_meta.sample_id} (${orig_meta.passage}) vs ${pass_meta.sample_id} (${pass_meta.passage})"
    }

    // Run line-based divergence analysis on normalized VCFs
    LINE_BASED_DIVERGENCE(grouped_samples)

    emit:
    divergence_metrics = LINE_BASED_DIVERGENCE.out.divergence_metrics
    normalized_vcfs = NORMALIZE_VCF.out.normalized_vcfs
}