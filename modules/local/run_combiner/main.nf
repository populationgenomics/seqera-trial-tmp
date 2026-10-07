process COMBINER_DATAPROC {
  container "${params.dataproc_container}"
  machineType 'e2-standard-2'

  input:
  path samplesheet
  val mc_hash

  output:
  tuple val(mc_hash), val(params.output_vds ?: "${params.outdir}/combiner/${mc_hash}.vds"), emit: result

  script:
  def vds_output = params.output_vds ?: "${params.outdir}/combiner/${mc_hash}.vds"
  def cohort_list = params.cohorts.replace(',', ' ')

  """
        # locate this script, staged in from the bin
        HAIL_SCRIPT="\$(which combiner_from_samplesheet.py)"

        submit_combiner_dataproc.py \\
            --cluster-name ${params.ar_guid} \\
            --region australia-southeast1 \\
            --project ${params.gcp_project} \\
            --script \${HAIL_SCRIPT} \\
            --staging-bucket ${params.dataproc_staging_bucket} \\
            --temp-bucket ${params.dataproc_temp_bucket} \\
            --autoscaling-policy ${params.dataproc_autoscaling_policy} \\
            --master-type ${params.dataproc_master_type} \\
            --worker-type ${params.dataproc_worker_type} \\
            --num-workers ${params.dataproc_num_workers} \\
            --num-secondary-workers ${params.dataproc_num_secondary_workers} \\
            --max-age ${params.dataproc_max_age} \\
            --max-idle ${params.dataproc_max_idle} \\
            --boot-disk-size ${params.dataproc_boot_disk_size_gb} \\
            --num-local-ssds ${params.dataproc_num_local_ssds} \\
            --service-account ${params.dataproc_service_account} \\
            --branch-factor ${params.branch_factor} \\
            --target-records ${params.target_records} \\
            --gvcf-batch-size ${params.gvcf_batch_size} \\
            ${params.preemptible_secondary_workers ? '--preemptible-secondary-workers' : ''} \\
            ${params.dataproc_internal_ip_only ? '--internal-ip-only' : ''} \\
            ${params.force_new ? '--force' : ''} \\
            --samplesheet ${samplesheet} \\
            --vds_output ${vds_output} \\
            --intervals ${params.intervals} \\
            --labels ar-guid=${params.ar_guid} compute-category="dataproc" dataset=${params.metamist_project}

        python3 -m cpg_utils.metamist_registration \\
             --project ${params.metamist_project} \\
             --output ${vds_output} \\
             --type vds \\
             --cohorts ${cohort_list} \\
             --meta stage=COMBINER_DATAPROC sequencing_type=genome ar-guid=${params.ar_guid}
        """
}

