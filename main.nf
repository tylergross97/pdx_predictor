#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { DATA_DOWNLOAD_WES } from './modules/data_download_wes.nf'

workflow {
    ch_wes_samplesheet = Channel.fromPath(params.wes_samplesheet)
    DATA_DOWNLOAD_WES(ch_wes_samplesheet)
}