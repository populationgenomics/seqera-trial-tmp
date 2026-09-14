#!/usr/bin/env python3

"""
Dataproc lifecycle orchestrator for running Hail jobs from Nextflow.
"""

import argparse
import json

from cpg_utils import to_path
from cpg_utils.dataproc_runner import (
    HailDataprocCluster,
    install_sigterm_cleanup,
    parse_label_kvs,
)


def parse_args():
    parser = argparse.ArgumentParser(description='Create a Dataproc cluster, run a Hail job, and tear down.')
    parser.add_argument('--cluster-name', required=True, help='Prefix for cluster name')
    parser.add_argument('--region', required=True, help='GCP region')
    parser.add_argument('--project', required=True, help='GCP project ID')
    parser.add_argument('--script', required=True, help='Local path to Hail Python script')
    parser.add_argument('--samplesheet', required=True, help='Path to sample sheet TSV')
    parser.add_argument('--staging-bucket', required=True, help='GCS bucket for staging')
    parser.add_argument('--temp-bucket', required=True, help='GCS bucket for Hail temp files')
    parser.add_argument('--vds_output', required=True, help='VDS path to write')
    parser.add_argument(
        '--intervals',
        required=True,
        help='Regions BED path (must be a gs:// URI, read directly by the Hail job).',
    )
    parser.add_argument(
        '--autoscaling-policy',
        default=None,
        help='Existing autoscaling policy ID or full resource URI.',
    )
    parser.add_argument(
        '--labels',
        nargs='*',
        default=[],
        help='Cluster labels as key=value pairs for billing/cost tracking.',
    )

    # Cluster sizing / lifecycle, forwarded to HailDataprocCluster.
    parser.add_argument('--master-type', default='n1-highmem-4', help='Master machine type.')
    parser.add_argument('--worker-type', default='n1-highmem-4', help='Worker machine type.')
    parser.add_argument('--num-workers', type=int, default=2, help='Number of primary workers.')
    parser.add_argument('--num-secondary-workers', type=int, default=2, help='Initial secondary workers.')
    parser.add_argument('--max-age', type=int, default=86400, help='Cluster auto-delete TTL (seconds).')
    parser.add_argument('--max-idle', type=int, default=900, help='Cluster idle-delete TTL (seconds).')
    parser.add_argument('--boot-disk-size', type=int, default=100, help='Per-node boot disk size (GB).')
    parser.add_argument(
        '--num-local-ssds', type=int, default=0, help='Local NVMe SSDs per node (375 GB each) for shuffle scratch.'
    )
    parser.add_argument(
        '--service-account', default=None, help='Service account for the cluster VMs (defaults to the compute SA).'
    )
    parser.add_argument(
        '--preemptible-secondary-workers',
        action='store_true',
        help='Use preemptible (spot) secondary workers.',
    )

    # Combiner tuning, forwarded to the Hail script.
    parser.add_argument('--branch-factor', type=int, default=50, help='Combiner merge-tree branch factor.')
    parser.add_argument('--target-records', type=int, default=200000, help='Target records per partition.')
    parser.add_argument('--gvcf-batch-size', type=int, default=50, help='gVCFs combined per batch.')
    parser.add_argument('--force', action='store_true', help='Force a fresh combine, ignoring any saved plan.')
    return parser.parse_args()


def main():
    args = parse_args()

    labels = parse_label_kvs(args.labels) if args.labels else None

    cluster = HailDataprocCluster(
        project=args.project,
        region=args.region,
        autoscaling_policy=args.autoscaling_policy,
        cluster_name_prefix=args.cluster_name,
        staging_bucket=args.staging_bucket,
        temp_bucket=args.temp_bucket,
        labels=labels,
        master_type=args.master_type,
        worker_type=args.worker_type,
        num_workers=args.num_workers,
        num_secondary_workers=args.num_secondary_workers,
        max_age_seconds=args.max_age,
        max_idle_seconds=args.max_idle,
        boot_disk_size_gb=args.boot_disk_size,
        preemptible_workers=args.preemptible_secondary_workers,
        num_local_ssds=args.num_local_ssds,
        service_account=args.service_account,
    )

    # Tear the cluster down promptly if Nextflow/Seqera cancels the task (SIGTERM).
    install_sigterm_cleanup(cluster)

    with cluster:
        script_uri = cluster.upload(args.script)
        samplesheet = cluster.upload(args.samplesheet)

        if not args.intervals.startswith('gs://'):
            raise ValueError('Supplied intervals must point to a GCP path')

        intervals_uri = args.intervals

        hail_args = [
            '--samplesheet',
            to_path(samplesheet).name,
            '--out',
            args.vds_output,
            '--intervals',
            intervals_uri,
            '--temp_dir',
            args.temp_bucket,
            '--branch',
            str(args.branch_factor),
            '--target_records',
            str(args.target_records),
            '--gvcf_batch',
            str(args.gvcf_batch_size),
        ]
        if args.force:
            hail_args.append('--force')

        submitted_job = cluster.run_job(
            script_uri,
            args=hail_args,
            file_uris=[samplesheet],
        )

        # collect this information in an output file
        output = {
            'cluster_name': cluster.name,
            'job_id': submitted_job.job_uuid,
            'region': args.region,
            'project': args.project,
            'hail_version': cluster.hail_version,
            'autoscaling_policy_uri': cluster.autoscaling_policy_uri,
            'labels': cluster.labels,
        }
    with open('output.json', 'w') as f:
        json.dump(output, f, indent=2)


if __name__ == '__main__':
    main()
