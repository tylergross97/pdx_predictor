#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { DATA_DOWNLOAD_WES } from './modules/data_download_wes.nf'
include { STAGE_REFERENCE } from './modules/stage_reference.nf'
include { LINE_BASED_DIVERGENCE_ANALYSIS } from './subworkflows/line_based_divergence_analysis.nf'
include { EXTRACT_FEATURES } from './modules/extract_features.nf'
include { MERGE_FEATURES } from './modules/merge_features.nf'

workflow {
    ch_wes_samplesheet = Channel.fromPath(params.wes_samplesheet)
        .splitCsv(header: true, quote: '"')
        .filter { row ->
            def pdmType = row['PDM Type']?.toString()?.strip()?.toLowerCase() ?: ''
            return pdmType == 'pdx' || pdmType == 'patient/originator specimen'
        }
        .map { row ->
            def determinePassage = { sampleId ->
                if (sampleId == 'ORIGINATOR') {
                    return 'o'
                }
                def cleanSampleId = sampleId.replaceAll(/(_RG-|RG-)/, '')
                def charCount = cleanSampleId.length()
                
                switch (charCount) {
                    case 3: return 'P0'
                    case 6: return 'P1'
                    case 9: return 'P2'
                    case 12: return 'P3'
                    case 15: return 'P4'
                    case 18: return 'P5'
                    case 21: return 'P6'
                    case 24: return 'P7'
                    case 27: return 'P8'
                    case 30: return 'P9'
                    default: return 'Unknown'
                }
            }
            
            def meta = [
                patient_id: row['Patient ID'],
                specimen_id: row['Specimen ID'], 
                sample_id: row['Sample ID'],
                passage: determinePassage(row['Sample ID']),
                onco_tree_code: row['OncoTreeCode']
            ]
            
            return [meta, row]
        }
    
    // Start reference download immediately (in parallel with VCF downloads)
    STAGE_REFERENCE()
    
    // Download VCF files
    DATA_DOWNLOAD_WES(ch_wes_samplesheet)
    
    // Wait for both reference and VCF files to be ready
    // This creates a channel where each VCF gets paired with the same reference
    vcf_with_reference = DATA_DOWNLOAD_WES.out.vcf_files
        .combine(STAGE_REFERENCE.out.reference)
        .map { meta, vcf_file, reference ->
            [meta, vcf_file, reference]
        }
    
    // Debug: Show when reference is ready
    STAGE_REFERENCE.out.reference.view { "Reference genome ready: ${it}" }
    
    // Run line-based genomic divergence analysis
    LINE_BASED_DIVERGENCE_ANALYSIS(vcf_with_reference)
    
    // Display results
    LINE_BASED_DIVERGENCE_ANALYSIS.out.divergence_metrics
        .view { patient_id, metrics_file ->
            "Genomic divergence metrics computed for patient ${patient_id}: ${metrics_file}"
        }
    
    // Optional: Collect all metrics into a summary
    LINE_BASED_DIVERGENCE_ANALYSIS.out.divergence_metrics
        .map { patient_id, metrics_file -> metrics_file }
        .collectFile(name: 'all_divergence_metrics.csv', keepHeader: true, storeDir: "${params.outdir_base}/summary")
        .view { "Summary metrics file created: ${it}" }

    // Extract features from ORIGINATOR VCF files
    EXTRACT_FEATURES(DATA_DOWNLOAD_WES.out.vcf_files.filter { meta, vcf_file -> meta.passage == 'o' })

    // Use toList() to preserve tuple structure
    collected_features = EXTRACT_FEATURES.out.toList()

    // Debug: view what we're collecting
    collected_features.view { "Collected features: $it" }

    // Separate metadata and files for the process
    metas = collected_features.map { tuples -> tuples.collect { tuple -> tuple[0] } }
    files = collected_features.map { tuples -> tuples.collect { tuple -> tuple[1] } }

    // Pass separated metadata and files to MERGE_FEATURES
    MERGE_FEATURES(metas, files)

    // Display the consolidated features file
    MERGE_FEATURES.out.consolidated_features
        .view { "Consolidated features file created: ${it}" }
}