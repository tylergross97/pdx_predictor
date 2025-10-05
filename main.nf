#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { DATA_DOWNLOAD_WES } from './modules/data_download_wes.nf'
include { LINE_BASED_DIVERGENCE_ANALYSIS } from './subworkflows/line_based_divergence_analysis.nf'

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
    
    // Download VCF files
    DATA_DOWNLOAD_WES(ch_wes_samplesheet)
    
    // Run line-based genomic divergence analysis
    LINE_BASED_DIVERGENCE_ANALYSIS(DATA_DOWNLOAD_WES.out.vcf_files)
    
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
}