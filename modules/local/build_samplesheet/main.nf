process BUILD_SAMPLESHEET {
    machineType 'e2-micro'

    container "${params.dataproc_container}"

    output:
        path("*.tsv"), emit: samplesheet

    script:
        def cohort_list = params.cohorts.replace(',', ' ')
        """
        fetch_cohort_samplesheet.py \\
            --project ${params.metamist_project} \\
            --cohorts ${cohort_list} \\
            ${params.expected_sg_count ? "--expected-sg-count ${params.expected_sg_count}" : ''}
        """
}