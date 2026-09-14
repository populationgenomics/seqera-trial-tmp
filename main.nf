#!/usr/bin/env nextflow

include { BUILD_SAMPLESHEET  } from './modules/local/build_samplesheet/main'
include { COMBINER_DATAPROC  } from './modules/local/run_combiner/main'


workflow {
    main:
    	if (!params.ar_guid) {
			println "Required --ar_guid argument not provided, jobs would not be labelled correctly"
			exit 1
		}

        if (!params.cohorts) {
            println "Required --cohorts argument not provided"
            exit 1
        }
        if (!params.metamist_project) {
            println "Required --metamist_project argument not provided"
            exit 1
        }

        BUILD_SAMPLESHEET()

        MC_HASH = BUILD_SAMPLESHEET.out.samplesheet
            .map { file -> file.baseName }

        // kick off a cheeky little combine
        COMBINER_DATAPROC(
            BUILD_SAMPLESHEET.out.samplesheet,
            MC_HASH,
        )

    publish:
        samplesheet = BUILD_SAMPLESHEET.out.samplesheet
        dp_output = COMBINER_DATAPROC.out.result
}

output {
    samplesheet {
    }
    dp_output {
        path { mc_hash, result -> "dataproc_run/${mc_hash}_outputs" }
    }
}