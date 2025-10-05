#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { DATA_DOWNLOAD_WES } from './modules/data_download_wes.nf'

workflow {
    ch_wes_samplesheet = Channel.fromPath(params.wes_samplesheet)
        .splitCsv(header: true, quote: '"')  // Handle quoted fields properly
        .filter { row ->
            // Only process PDX samples or Patient/Originator Specimen samples
            def pdmType = row['PDM Type']?.toString()?.strip()?.toLowerCase() ?: ''
            return pdmType == 'pdx' || pdmType == 'patient/originator specimen'
        }
        .map { row ->
            // Function to determine passage based on Sample ID
            def determinePassage = { sampleId ->
                if (sampleId == 'ORIGINATOR') {
                    return 'o'  // ORIGINATOR gets 'o'
                }
                
                // Remove '_RG-' or 'RG-' prefixes/suffixes for character counting
                def cleanSampleId = sampleId.replaceAll(/(_RG-|RG-)/, '')
                
                // Count all characters in the cleaned sample ID
                def charCount = cleanSampleId.length()
                
                // Map character count to passage (multiples of 3)
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
            
            // Create meta map with the requested fields
            def meta = [
                patient_id: row['Patient ID'],
                specimen_id: row['Specimen ID'], 
                sample_id: row['Sample ID'],
                passage: determinePassage(row['Sample ID']),
                onco_tree_code: row['OncoTreeCode']
            ]
            
            return [meta, row]
        }
    
    DATA_DOWNLOAD_WES(ch_wes_samplesheet)
}