#!/usr/bin/env python3

"""
Quick combiner-in-dataproc-in-nextflow test run
"""

import argparse
import gzip
import os
from csv import DictReader

import hail as hl
from hail.vds import new_combiner

from cpg_utils import to_path


def get_arguments() -> argparse.Namespace:
    """Parse command-line arguments"""

    parser = argparse.ArgumentParser()
    parser.add_argument('--samplesheet', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--temp_dir', required=True)
    parser.add_argument('--intervals', required=True)
    parser.add_argument('--force', action='store_true')
    parser.add_argument('--branch', type=int, default=50)
    parser.add_argument('--gvcf_batch', type=int, default=50)
    parser.add_argument('--target_records', type=int, default=200000)
    parser.add_argument('--prev_vds', default=None)
    return parser.parse_args()


def parse_samplesheet(ss_path: str) -> dict[str, str]:
    """Parse sample sheet file"""
    sample_to_gvcf: dict[str, str] = {}
    with open(ss_path) as fh:
        reader = DictReader(fh, delimiter='\t')
        for row in reader:
            sample_to_gvcf[row['sg_id']] = row['gvcf']
    return sample_to_gvcf


def read_bed_file_as_intervals(bed_path: str) -> list[hl.Interval]:
    """Manually interpret an input BED file (plain-text or gzipped) as a series of Intervals."""
    path = to_path(bed_path)
    intervals: list[hl.Interval] = []
    with gzip.open(path.open('rb'), 'rt') if str(path).endswith('.gz') else path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue

            chrom, start, end = stripped.split()[:3]

            start_locus = hl.Locus(chrom, int(start) + 1, reference_genome='GRCh38')
            end_locus = hl.Locus(chrom, int(end), reference_genome='GRCh38')
            intervals.append(hl.Interval(start_locus, end_locus, includes_start=True, includes_end=True))
    return intervals


def main() -> None:
    hl.init(default_reference='GRCh38')

    args = get_arguments()

    if to_path(args.out).exists():
        print(f'Output {args.out} already exists, doing nothing.')
        return

    # get the multicohort hash, derived from the samplesheet name. Saves passing two things...
    name_hash = to_path(args.samplesheet).stem

    combiner_temp = os.path.join(args.temp_dir, 'combiner_temp', name_hash)

    # get the mapping of samples and files relevant to this analysis
    sg_to_gvcf = parse_samplesheet(args.samplesheet)

    # sort by SG ID so names and paths stay aligned
    gvcf_sample_names = sorted(sg_to_gvcf)
    gvcfs = [sg_to_gvcf[sgid] for sgid in gvcf_sample_names]
    # for now, skip prior VDS integration, or considering removal from an existing VDS

    plan = os.path.join(args.temp_dir, 'combiner_plan', f'{name_hash}.json')

    intervals = read_bed_file_as_intervals(args.intervals)

    # new_combiner resumes from save_path itself, reapplying the tuning args over the saved plan
    combiner = new_combiner(
        output_path=args.out,
        save_path=plan,
        gvcf_paths=gvcfs,
        gvcf_sample_names=gvcf_sample_names,
        gvcf_external_header=gvcfs[0],
        reference_genome='GRCh38',
        temp_path=combiner_temp,
        intervals=intervals,
        force=args.force,
        branch_factor=args.branch,
        target_records=args.target_records,
        gvcf_batch_size=args.gvcf_batch,
    )
    combiner.run()


if __name__ == '__main__':
    main()
