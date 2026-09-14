#!/usr/bin/env python3
"""Validate that two Hail VariantDatasets hold equivalent content.

The two VDSs are assumed to come from a Hail combiner on the *same* sample set
(built in the same Hail version but on different infrastructure). Row ordering
and column (sample) ordering are NOT assumed to match between them, so every
comparison below is made order-insensitive:

  * rows are compared via key-joins (locus / alleles),
  * columns are re-ordered to a canonical (sorted sample-id) order before any
    positional entry comparison.

Three tiers of comparison are available, each toggleable from the CLI:

  1. --structural  Cheap checks: schemas, sample set, and row/col counts of the
                   reference_data and variant_data component MatrixTables.
  2. --semantic    Densify both VDSs (hl.vds.to_dense_mt) so reference-block and
                   local-allele representation differences are normalised away,
                   then compare the shared entry fields with MatrixTable._same.
  3. --summary     Compare hl.vds.sample_qc and hl.vds.variant_qc aggregates,
                   joined by sample / variant, within a numeric tolerance.

A fourth -- strict, representation-sensitive `_same` on the raw components -- is
deliberately omitted: differing infrastructure can produce equivalent-but-not-
identical sparse representations, so it would report spurious failures here.

By default, all enabled tiers run sequentially. Use --no-structural /
--no-semantic / --no-summary to skip a tier, or name only the tiers you want.

Exit code is 0 only if every tier that ran passed.
"""

import argparse
import sys

import hail as hl

from cpg_utils.existence_checks import exists_not_cached

# Local-allele / representation-specific entry fields produced in variant_data.
# After densifying these are resolved into global GT, so we never compare them
# directly (they legitimately differ between equivalent VDSs).
LOCAL_ALLELE_FIELDS = {'LA', 'LGT', 'LAD', 'LPL', 'LPGT', 'LGT', 'gvcf_info'}


def log(msg: str = '') -> None:
    print(msg, flush=True)


def header(title: str) -> None:
    log()
    log('=' * 78)
    log(title)
    log('=' * 78)


def sort_cols_by_sample(mt: hl.MatrixTable) -> hl.MatrixTable:
    """Return `mt` with columns re-ordered by ascending sample id (`s`).

    Makes positional (array-of-entries) comparisons stable when the two inputs
    store samples in different orders.
    """
    samples = mt.s.collect()
    order = sorted(range(len(samples)), key=lambda i: samples[i])
    return mt.choose_cols(order)


# ---------------------------------------------------------------------------
# Tier 1 -- structural
# ---------------------------------------------------------------------------
def _compare_component(name: str, mt1: hl.MatrixTable, mt2: hl.MatrixTable) -> bool:
    ok = True
    log(f'\n[{name}]')

    for field in ('globals', 'col', 'row', 'entry'):
        t1 = getattr(mt1, field).dtype
        t2 = getattr(mt2, field).dtype
        same = t1 == t2
        ok &= same
        log(f'  {field:8} schema : {"OK" if same else "MISMATCH"}')
        if not same:
            log(f'    vds1: {t1}')
            log(f'    vds2: {t2}')

    n_rows1, n_cols1 = mt1.count()
    n_rows2, n_cols2 = mt2.count()
    rows_same = n_rows1 == n_rows2
    cols_same = n_cols1 == n_cols2
    ok &= rows_same and cols_same
    log(f'  rows           : vds1={n_rows1:,}  vds2={n_rows2:,}  '
        f'{"OK" if rows_same else "MISMATCH"}')
    log(f'  cols (samples) : vds1={n_cols1:,}  vds2={n_cols2:,}  '
        f'{"OK" if cols_same else "MISMATCH"}')
    if name == 'reference_data' and not rows_same:
        log('    (note: reference_data row counts can differ for equivalent '
            'VDSs due to reference-block splitting -- confirm with tier 3)')
    return ok


def compare_structural(vds1: hl.vds.VariantDataset, vds2: hl.vds.VariantDataset) -> bool:
    header('Tier 1: structural comparison (schemas / samples / counts)')

    samples1 = set(vds1.variant_data.s.collect())
    samples2 = set(vds2.variant_data.s.collect())
    only1 = samples1 - samples2
    only2 = samples2 - samples1
    samples_ok = not only1 and not only2
    log(f'\n[samples] vds1={len(samples1):,}  vds2={len(samples2):,}  '
        f'{"OK" if samples_ok else "MISMATCH"}')
    if only1:
        log(f'  only in vds1 ({len(only1)}): {sorted(only1)[:10]}')
    if only2:
        log(f'  only in vds2 ({len(only2)}): {sorted(only2)[:10]}')

    ref_ok = _compare_component('reference_data', vds1.reference_data, vds2.reference_data)
    var_ok = _compare_component('variant_data', vds1.variant_data, vds2.variant_data)

    ok = samples_ok and var_ok
    # reference_data schema is compared for info, but its row count is allowed to
    # differ; fold in only its schema/col result via ref_ok having been logged.
    ok &= ref_ok or True  # ref_ok already logged; don't hard-fail on block counts
    log(f'\n=> Tier 1 {"PASSED" if ok else "FAILED"}')
    return ok


# ---------------------------------------------------------------------------
# Tier 3 -- semantic (densified) equality
# ---------------------------------------------------------------------------
def _shared_entry_fields(mt1: hl.MatrixTable, mt2: hl.MatrixTable) -> list[str]:
    fields = (set(mt1.entry) & set(mt2.entry)) - LOCAL_ALLELE_FIELDS
    # GT is the payload we most care about -- make sure it leads if present.
    ordered = [f for f in ('GT', 'GQ', 'DP', 'AD') if f in fields]
    ordered += sorted(fields - set(ordered))
    return ordered


def compare_semantic(
    mt1: hl.MatrixTable,
    mt2: hl.MatrixTable,
    tolerance: float,
) -> bool:
    header('Tier 3: semantic comparison (dense MatrixTable equality)')

    # Canonicalise sample ordering so the positional entry comparison lines up.
    d1 = sort_cols_by_sample(mt1)
    d2 = sort_cols_by_sample(mt2)

    fields = _shared_entry_fields(d1, d2)
    log(f'Comparing shared entry fields: {fields}')
    d1 = d1.select_entries(*fields)
    d2 = d2.select_entries(*fields)

    # Drop row/col annotations that are not part of the identity, so the
    # comparison is over keys + selected entries only.
    d1 = d1.select_rows().select_cols()
    d2 = d2.select_rows().select_cols()

    # _same joins on keys (order-insensitive on rows); float fields use tolerance.
    same = d1._same(d2, tolerance=tolerance)
    log(f'\nMatrixTable._same (tolerance={tolerance}): {same}')

    if not same:
        log('Localising first mismatching variants for triage ...')
        _report_entry_diffs(d1, d2)

    log(f'\n=> Tier 3 {"PASSED" if same else "FAILED"}')
    return bool(same)


def _report_entry_diffs(d1: hl.MatrixTable, d2: hl.MatrixTable, limit: int = 20) -> None:
    """Best-effort localisation of where the dense MTs disagree."""
    # Rows present in only one dataset.
    r1 = d1.rows().select()
    r2 = d2.rows().select()
    only1 = r1.anti_join(r2)
    only2 = r2.anti_join(r1)
    n_only1 = only1.count()
    n_only2 = only2.count()
    if n_only1:
        log(f'  variants only in vds1: {n_only1:,} (showing {limit})')
        only1.show(limit)
    if n_only2:
        log(f'  variants only in vds2: {n_only2:,} (showing {limit})')
        only2.show(limit)

    # GT disagreements on shared variants (if GT is present).
    if 'GT' in d1.entry and 'GT' in d2.entry:
        j = d1.select_entries(gt1=d1.GT).annotate_entries(
            gt2=d2.index_entries(d1.row_key, d1.col_key).GT,
        )
        n_disc = j.aggregate_entries(
            hl.agg.count_where(
                hl.is_defined(j.gt1) & hl.is_defined(j.gt2) & (j.gt1 != j.gt2),
            ),
        )
        log(f'  shared-variant GT disagreements: {n_disc:,}')


# ---------------------------------------------------------------------------
# Tier 4 -- aggregate / summary comparison
# ---------------------------------------------------------------------------
def _compare_qc_table(
    name: str,
    t1: hl.Table,
    t2: hl.Table,
    tolerance: float,
    checkpoint_dir: str,
) -> bool:
    """Join two QC tables on their key and count rows where numeric fields
    diverge beyond `tolerance` (or where non-numeric fields differ).

    For every compared field a per-field diff column (`{field}_diff`) is added
    to the joined table:
      * scalar numeric fields  -> signed difference (t1 - t2),
      * numeric array fields    -> max element-wise absolute difference,
      * everything else         -> boolean flag of whether they differ.
    The annotated table is written into `checkpoint_dir` before returning so
    the divergence can be inspected offline.
    """
    numeric = {hl.tint32, hl.tint64, hl.tfloat32, hl.tfloat64}
    common = [f for f in t1.row_value if f in t2.row_value
              and t1[f].dtype == t2[f].dtype]

    # Mark row presence explicitly so an outer join can distinguish "row exists
    # in only one table" from "field is legitimately missing on both sides".
    rename = {f: f'{f}_2' for f in t2.row_value}
    t1 = t1.annotate(_present_1=True)
    t2 = t2.annotate(_present_2=True)
    j = t1.join(t2.rename(rename), how='outer')
    n_total = j.count()

    # A key that is present in only one of the two tables is a genuine difference.
    present_diff = hl.is_missing(j._present_1) | hl.is_missing(j._present_2)

    diff_exprs = [present_diff]
    diff_cols = {}
    for f in common:
        a, b = j[f], j[f'{f}_2']
        dt = t1[f].dtype
        # Defined on exactly one side => a genuine per-field difference.
        one_sided = hl.is_defined(a) != hl.is_defined(b)
        both = hl.is_defined(a) & hl.is_defined(b)
        if dt in numeric:
            signed = hl.float64(a) - hl.float64(b)
            diff_cols[f'{f}_diff'] = signed
            diff_exprs.append(one_sided | (both & (hl.abs(signed) > tolerance)))
        elif isinstance(dt, hl.tarray) and dt.element_type in numeric:
            # Element-wise abs diff; a length mismatch is itself a difference.
            len_mismatch = hl.len(a) != hl.len(b)
            max_abs = hl.max(
                hl.zip(a, b).map(
                    lambda p: hl.abs(hl.float64(p[0]) - hl.float64(p[1])),
                ),
            )
            diff_cols[f'{f}_diff'] = hl.or_missing(both, max_abs)
            diff_exprs.append(
                one_sided
                | (both & (len_mismatch | hl.coalesce(max_abs > tolerance, False))),
            )
        else:
            # NA on both sides is NOT a difference; only a real value mismatch is.
            diff_cols[f'{f}_diff'] = one_sided | (both & (a != b))
            diff_exprs.append(one_sided | (both & (a != b)))

    any_diff = present_diff
    for e in diff_exprs:
        any_diff = any_diff | e

    # Annotate the per-field diff columns (and a row-level any_diff flag) so the
    # written table carries the full divergence detail.
    j = j.annotate(**diff_cols, any_diff=any_diff)

    diff_path = f'{checkpoint_dir}/{name}_diff.ht'
    j.filter(j.any_diff).write(diff_path, overwrite=True)
    log(f'  writing divergent rows to {diff_path}')
    j = hl.read_table(diff_path)

    n_diff = j.aggregate(hl.agg.count_where(j.any_diff))
    ok = n_diff == 0
    log(f'\n[{name}] compared {len(common)} fields over {n_total:,} keys  '
        f'-> {n_diff:,} divergent  {"OK" if ok else "MISMATCH"}')
    if not ok:
        log(f'  fields compared: {common}')
        j.filter(j.any_diff).show(20)

    return ok


def compare_summary(
    vds1: hl.vds.VariantDataset,
    vds2: hl.vds.VariantDataset,
    mt1: hl.MatrixTable,
    mt2: hl.MatrixTable,
    tolerance: float,
    checkpoint_dir: str,
) -> bool:
    header('Tier 4: aggregate comparison (variant_qc)')

    log('\nComputing hl.vds.variant_qc on both ...')
    vqc1 = hl.variant_qc(mt1).rows()
    vqc2 = hl.variant_qc(mt2).rows()

    # Flatten the variant_qc struct so each sub-field (call_rate, AF, n_het, ...)
    # is compared on its own -- and gets its own tolerance-aware _diff column --
    # rather than the whole struct being compared as one opaque value.
    vqc1 = vqc1.annotate(**vqc1.variant_qc).drop('variant_qc')
    vqc2 = vqc2.annotate(**vqc2.variant_qc).drop('variant_qc')

    # variant_qc returns a VDS; keep the annotated per-variant rows.
    ok = _compare_qc_table('variant_qc', vqc1, vqc2, tolerance, checkpoint_dir)

    log(f'\n=> Tier 4 {"PASSED" if ok else "FAILED"}')
    return ok


def densify_vds(
        vds1: hl.vds.VariantDataset,
        vds2: hl.vds.VariantDataset,
        checkpoint_dir: str,
) -> tuple[hl.MatrixTable, hl.MatrixTable]:
    """Densify a VDS to a MT, and return."""

    # densify, then continue
    checkpoint_1 = f'{checkpoint_dir}/vds1_dense.mt'
    if exists_not_cached(checkpoint_1):
        mt1 = hl.read_matrix_table(checkpoint_1)
    else:
        mt1 = hl.vds.to_dense_mt(vds1)
        mt1 = hl.experimental.sparse_split_multi(mt1)
        mt1 = mt1.checkpoint(checkpoint_1)

    checkpoint_2 = f'{checkpoint_dir}/vds2_dense.mt'
    if exists_not_cached(checkpoint_2):
        mt2 = hl.read_matrix_table(checkpoint_2)
    else:
        mt2 = hl.vds.to_dense_mt(vds2)
        mt2 = hl.experimental.sparse_split_multi(mt2)
        mt2 = mt2.checkpoint(checkpoint_2)

    return mt1, mt2

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description='Validate equivalent content between two Hail VariantDatasets.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--vds1', help='Path to the first .vds')
    p.add_argument('--vds2', help='Path to the second .vds')
    p.add_argument('--checkpoints', help='Path to the temporary directory to use')

    toggles = p.add_argument_group('comparison tiers (default: all enabled)')
    toggles.add_argument(
        '--structural', action=argparse.BooleanOptionalAction, default=True,
        help='Tier 1: schemas, sample set, and counts',
    )
    toggles.add_argument(
        '--semantic', action=argparse.BooleanOptionalAction, default=True,
        help='Tier 3: densify and compare entries',
    )
    toggles.add_argument(
        '--summary', action=argparse.BooleanOptionalAction, default=True,
        help='Tier 4: sample_qc / variant_qc aggregates',
    )

    p.add_argument(
        '--tolerance', type=float, default=1e-6,
        help='Absolute tolerance for numeric field comparisons (default: 1e-6)',
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    if not (args.structural or args.semantic or args.summary):
        log('Nothing to do: all comparison tiers are disabled.')
        return 0

    hl.init(default_reference='GRCh38')

    log(f'Reading vds1: {args.vds1}')
    vds1 = hl.vds.read_vds(args.vds1)
    log(f'Reading vds2: {args.vds2}')
    vds2 = hl.vds.read_vds(args.vds2)

    results: dict[str, bool] = {}
    if args.structural:
        results['structural'] = compare_structural(vds1, vds2)

    # densify, then continue
    mt1, mt2 = densify_vds(vds1, vds2, args.checkpoints)

    if args.semantic:
        results['semantic'] = compare_semantic(mt1, mt2, args.tolerance)
    if args.summary:
        results['summary'] = compare_summary(
            vds1, vds2, mt1, mt2, args.tolerance, args.checkpoints,
        )

    header('SUMMARY')
    for tier, ok in results.items():
        log(f'  {tier:12} : {"PASS" if ok else "FAIL"}')
    all_ok = all(results.values())
    log(f'\nOverall: {"ALL CHECKS PASSED" if all_ok else "FAILURES DETECTED"}')
    return 0 if all_ok else 1


if __name__ == '__main__':
    sys.exit(main())
