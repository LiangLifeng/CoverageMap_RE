"""Load immutable alignment or saved-bundle inputs for concordance."""
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import tempfile
import time

from bin.coverage_core import extract_pool, merge_references, runtime_pool

POLICY_DEFAULTS = dict(mate_policy='both', sampling_unit='read', min_mapq=0,
                       exclude_duplicates=False, exclude_qcfail=False)


def file_hash(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def resolve_bundle(path):
    path = Path(path).resolve()
    if path.is_dir():
        return path
    if path.is_file() and path.name == 'storage.json':
        return path.parent
    raise ValueError('reference-bundle must be a complete bundle directory or its storage.json')


def resolve_policies(args, saved=None):
    complete = saved is not None and 'coverage_alignment_policy' in saved
    if saved is not None:
        required = set(POLICY_DEFAULTS)-({'sampling_unit'} if complete else set())
        allowed = set(POLICY_DEFAULTS)|({'coverage_alignment_policy'} if complete else set())
        if not required <= set(saved) <= allowed:
            raise ValueError('Reference bundle has incomplete extraction policies')
        if (saved['mate_policy'] not in ('both', 'read1', 'read2')
                or saved.get('sampling_unit','read') not in ('read', 'fragment')
                or not isinstance(saved['min_mapq'], int) or saved['min_mapq'] < 0
                or not isinstance(saved['exclude_duplicates'], bool)
                or not isinstance(saved['exclude_qcfail'], bool)
                or (saved.get('sampling_unit') == 'fragment' and saved['mate_policy'] != 'both')
                or complete and saved['coverage_alignment_policy'] != 'primary_only'):
            raise ValueError('Invalid saved extraction policies')
    policies = dict(POLICY_DEFAULTS,**(saved or {}))
    for key in policies:
        value = getattr(args, key, None)
        if value is not None:
            if saved is not None and value != policies[key] and not (complete and key=='sampling_unit'):
                raise ValueError(f'Policy {key} conflicts with reference bundle')
            policies[key] = value
    return policies


@dataclass
class ConcordanceInputs:
    reference_pool: object
    test_pool: object
    layout: object
    policies: dict
    identities: dict
    bridge_metadata: list
    bundle_info: object
    comparison_scope: str
    derived_layout: bool
    bundle_directory: object
    layout_comparison: object = None

    def verify(self):
        if any(file_hash(path) != expected for path, expected in self.identities.items()):
            raise ValueError('Input changed during run')
        if self.bundle_directory is not None:
            if ((self.bundle_directory / 'status/FAILED').exists()
                    or (self.bundle_directory / 'status/RUNNING').exists()):
                raise ValueError('Reference bundle status changed during run')


def load_inputs(args):
    identities = {str(p): file_hash(p) for p in
                  (args.alignment, args.test_alignment, args.reference or args.reference_index)
                  if p is not None}
    bundle_info, bundle_directory = None, None
    bridge_metadata = []
    reference_pool, layout = None, None
    saved = None
    source_identity = identities.get(str(args.alignment))
    if args.reference_bundle is not None:
        from bin.reference_store import load_reference_bundle, safe_file
        bundle_directory = resolve_bundle(args.reference_bundle)
        tick = time.perf_counter()
        for name in ('checksums.json', 'status/SUCCESS', 'status/exit_code.txt'):
            path = safe_file(bundle_directory, name)
            identities[str(path)] = file_hash(path)
        reference_pool, layout, manifest = load_reference_bundle(bundle_directory, args.workers)
        for name, expected in json.loads((bundle_directory / 'checksums.json').read_text()).items():
            identities[str(safe_file(bundle_directory, name))] = expected
        saved = manifest['policies']
        sources = manifest['sources']
        source_identity = sources[0].get('sha256') if sources else None
        bundle_info = dict(path=str(bundle_directory), schema=manifest['schema'],
                           sample_id=manifest['sample_id'],
                           manifest_sha256=identities[str(bundle_directory / 'manifest.json')],
                           source_alignment_sha256=source_identity,
                           load_and_verify_seconds=time.perf_counter() - tick)
    policies = resolve_policies(args, saved)
    if reference_pool is not None:
        reference_pool = runtime_pool(reference_pool,saved,policies)
    derived = args.reference is None and args.reference_index is None
    if derived and bundle_info is None:
        raise ValueError('Alignment mode requires --reference or --reference-index')
    def extract(alignment, reference_index):
        ordering = {'sort_contigs': True} if layout is not None and layout.sort_policy else {}
        if args.extraction_python is None:
            return extract_pool(alignment, args.reference, reference_index=reference_index, **ordering, **policies)
        from bin.pool_bridge import extract_with_python
        args.outdir.parent.mkdir(parents=True, exist_ok=True)
        pool, ref, metadata = extract_with_python(args.extraction_python, alignment,
                    args.reference, reference_index, dict(policies, **ordering), args.outdir.parent)
        bridge_metadata.append(metadata)
        return pool, ref
    def extract_both(reference_index):
        nonlocal reference_pool, layout
        if reference_pool is None:
            reference_pool, layout = extract(args.alignment, reference_index)
        if args.alignment is not None and args.alignment.resolve() == args.test_alignment.resolve():
            return reference_pool, layout
        return extract(args.test_alignment, reference_index)
    if derived:
        args.outdir.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='bundle-layout-', dir=args.outdir.parent) as temporary:
            table = Path(temporary) / 'layout.tsv'
            table.write_text(''.join(f'{name}\t{length}\n' for name, length in zip(layout.names, layout.lengths)))
            test_pool, test_layout = extract_both(table)
    else:
        test_pool, test_layout = extract_both(args.reference_index)
    merge_references(layout, test_layout)
    layout_comparison = dict(reference_missing_header_contigs=list(layout.missing_header_contigs),
        test_missing_header_contigs=list(test_layout.missing_header_contigs),
        identity_unverified_contigs=[name for name, first, second in zip(layout.names, layout.md5, test_layout.md5)
                                    if not first or not second],
        policy='full-reference-geometry-known-m5-consistency-v1')
    scope = ('same-source-self-consistency' if source_identity is not None
             and source_identity == identities[str(args.test_alignment)] else 'separate-input-pools')
    result = ConcordanceInputs(reference_pool, test_pool, layout, policies, identities,
                               bridge_metadata, bundle_info, scope, derived, bundle_directory, layout_comparison)
    result.verify()
    return result
