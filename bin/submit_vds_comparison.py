#!/usr/bin/env python3

"""
Dataproc lifecycle orchestrator for running Hail jobs from Nextflow.
"""

import argparse

from cpg_utils.dataproc_runner import HailDataprocCluster, parse_label_kvs


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create a Dataproc cluster, run a Hail job, and tear down."
    )
    parser.add_argument("--cluster-name", required=True, help="Prefix for cluster name")
    parser.add_argument("--region", required=True, help="GCP region")
    parser.add_argument("--project", required=True, help="GCP project ID")
    parser.add_argument("--script", required=True, help="Local path to Hail Python script")
    parser.add_argument("--staging-bucket", required=True, help="GCS bucket for staging")
    parser.add_argument("--temp-bucket", required=True, help="GCS bucket for Hail temp files")
    parser.add_argument("--checkpoints", required=True, help="GCS bucket for checkpointed MTs.")
    parser.add_argument("--vds1", required=True, help="First VDS in comparison")
    parser.add_argument("--vds2", required=True, help="Second VDS in comparison")
    parser.add_argument(
        "--autoscaling-policy",
        default=None,
        help="Existing autoscaling policy ID or full resource URI.",
    )
    parser.add_argument(
        "--labels",
        nargs="*",
        default=[],
        help="Cluster labels as key=value pairs for billing/cost tracking.",
    )
    parser.add_argument(
        "--internal-ip-only",
        action="store_true",
        help="Launch cluster VMs with internal IPs only (no external IPs).",
    )
    return parser.parse_args()

def main():
    args = parse_args()

    labels = parse_label_kvs(args.labels) if args.labels else None

    with HailDataprocCluster(
            project=args.project,
            region=args.region,
            autoscaling_policy=args.autoscaling_policy,
            cluster_name_prefix=args.cluster_name,
            staging_bucket=args.staging_bucket,
            temp_bucket=args.temp_bucket,
            labels=labels,
            num_secondary_workers=2,
            preemptible_workers=False,
            internal_ip_only=args.internal_ip_only,
    ) as cluster:
        script_uri = cluster.upload(args.script)
        _submitted_job = cluster.run_job(
            script_uri,
            args=[
                '--vds1', args.vds1,
                '--vds2', args.vds2,
                '--checkpoints', args.checkpoints,
            ],
        )


if __name__ == "__main__":
    main()
