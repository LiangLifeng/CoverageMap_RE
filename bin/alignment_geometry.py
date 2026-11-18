"""Lossless, disk-spooled alignment relationships and runtime unit projections."""
from array import array
from collections import OrderedDict
from contextlib import closing
import errno
import hashlib
from itertools import groupby
import json
import os
from pathlib import Path
import shutil
import sqlite3
import sys
import tempfile
import time

import numpy as np

SCHEMA = 'coverage-geometry-2'
FRAGMENT = np.dtype([('identity', 'S64'), ('rg', '<i8'), ('flags', '<u4'), ('read1', '<i8'), ('read2', '<i8')])
READ = np.dtype([(name, '<i8') for name in ('fragment', 'role', 'flags', 'rg', 'mate_contig', 'mate_start', 'mate_unmapped', 'mate_reverse', 'audit_flags')])
SEGMENT = np.dtype([(name, '<i8') for name in (
    'read', 'type', 'contig', 'start', 'end', 'query_start', 'query_end', 'query_length',
    'reverse', 'mapq', 'flags', 'cigar', 'mate_contig', 'mate_start', 'mate_unmapped',
    'mate_reverse', 'sa', 'alignment_score', 'edit_distance', 'has_alignment_score', 'has_edit_distance',
    'eligible', 'block_begin', 'block_end')])
BLOCK = np.dtype([(name, '<i8') for name in ('segment', 'contig', 'start', 'end')])
STRING = np.dtype([('offset', '<i8'), ('length', '<i8')])
# Nonexclusive fragment audit bits.
SAME_PAIR, CROSS_PAIR, SINGLE_MAPPED, BOTH_UNMAPPED = 1, 2, 4, 8
UNPAIRED, INCOMPLETE, SPLIT_SAME, SPLIT_CROSS, MISSING_RG = 16, 32, 64, 128, 256
MATE_DISAGREEMENT, QUERY_OVERLAP, QUERY_LENGTH_CONFLICT = 512, 1024, 2048


def _spool_settings():
    """Execution-only scratch settings; never part of scientific identities."""
    parent = Path(os.environ.get('COVERAGEMAP_SPOOL_DIR') or
                  os.environ.get('SLURM_TMPDIR') or tempfile.gettempdir())
    if not parent.is_absolute() or not parent.is_dir():
        raise ValueError(f'Capture scratch must be an existing absolute directory: {parent}')
    try:
        reserve_mib = int(os.environ.get('COVERAGEMAP_SPOOL_MIN_FREE_MIB', '1024'))
    except ValueError as error:
        raise ValueError('COVERAGEMAP_SPOOL_MIN_FREE_MIB must be a nonnegative integer') from error
    if reserve_mib < 0:
        raise ValueError('COVERAGEMAP_SPOOL_MIN_FREE_MIB must be a nonnegative integer')
    return parent, reserve_mib * 1024 * 1024


def _check_spool_space(parent, reserve):
    free = shutil.disk_usage(parent).free
    if free < reserve:
        raise OSError(errno.ENOSPC,
            f'Capture scratch has {free} free bytes, below reserve {reserve}; '
            'choose COVERAGEMAP_SPOOL_DIR or reduce concurrent workers', str(parent))


def query_interval(cigar, reverse=False):
    """Original-read coordinates include hard clips; D/N/P consume no query."""
    cigar = cigar or []
    if not cigar:
        return -1, -1, -1
    if any(op not in range(9) or length <= 0 for op, length in cigar):
        raise ValueError('Invalid CIGAR operation')
    total = sum(n for op, n in cigar if op in (0, 1, 4, 5, 7, 8))
    leading = 0
    for op, n in cigar:
        if op not in (4, 5):
            break
        leading += n
    trailing = 0
    for op, n in reversed(cigar):
        if op not in (4, 5):
            break
        trailing += n
    start, end = leading, total - trailing
    if start > end:
        raise ValueError('Invalid clipped query interval')
    return (total - end, total - start, total) if reverse else (start, end, total)


def _string_interner(db, budget_bytes=1024 * 1024, capacity=4096):
    """Bound repeated RG/CIGAR/tag lookups without retaining read identities."""
    cache = OrderedDict()
    next_id, used_bytes = 0, 0
    def intern(value):
        nonlocal next_id, used_bytes
        value = str(value)
        found = cache.get(value)
        if found is not None:
            cache.move_to_end(value)
            return found[0]
        row = db.execute('SELECT id FROM strings WHERE value=?', (value,)).fetchone()
        if row is None:
            number = next_id
            db.execute('INSERT INTO strings VALUES (?,?)', (number, value))
            next_id += 1
        else:
            number = row[0]
        charge = sys.getsizeof(value) + 192
        if charge <= budget_bytes:
            while cache and (len(cache) >= capacity or used_bytes + charge > budget_bytes):
                _, (_, old_charge) = cache.popitem(last=False)
                used_bytes -= old_charge
            cache[value] = (number, charge)
            used_bytes += charge
        return number
    return intern


def _export(db, directory):
    for table, dtype, query in (
        ('fragments', FRAGMENT, 'SELECT identity,rg,flags,read1,read2 FROM fragments ORDER BY id'),
        ('reads', READ, 'SELECT '+','.join(READ.names)+' FROM reads ORDER BY id'),
        ('segments', SEGMENT, 'SELECT ' + ','.join(SEGMENT.names) + ' FROM segments ORDER BY id'),
        ('blocks', BLOCK, 'SELECT segment,contig,start,end FROM blocks ORDER BY id'),
    ):
        count = db.execute(f'SELECT count(*) FROM {table}').fetchone()[0]
        data = np.lib.format.open_memmap(directory / (table + '.npy'), mode='w+', dtype=dtype, shape=(count,))
        for i, row in enumerate(db.execute(query)):
            data[i] = tuple(row)
        data.flush()
        del data
    count = db.execute('SELECT count(*) FROM strings').fetchone()[0]
    data = np.lib.format.open_memmap(directory/'strings.npy', mode='w+', dtype=STRING, shape=(count,))
    with (directory/'strings.utf8').open('wb') as handle:
        handle.write(b'\0')  # Nonempty payload also for an empty alignment file.
        for i, value in db.execute('SELECT id,value FROM strings ORDER BY id'):
            encoded = value.encode('utf-8')
            data[i] = (handle.tell(), len(encoded))
            handle.write(encoded)
    data.flush()


def _audit_relations(db):
    """One ordered segment stream; at most four read roles per fragment."""
    # Pin the read-index loop order; any remaining sort is per read.
    rows = db.execute('''
        SELECT r.fragment,f.flags,r.id,r.role,
               r.mate_contig,r.mate_start,r.mate_unmapped,r.mate_reverse,
               s.type,s.contig,s.flags,s.query_start,s.query_end,s.query_length,s.start
        FROM reads r CROSS JOIN fragments f ON f.id=r.fragment
        CROSS JOIN segments s ON s.read=r.id
        ORDER BY r.fragment,r.role,s.query_start,s.query_end
    ''')
    read_updates, fragment_updates = [], []
    for fid, segments in groupby(rows, key=lambda row: row[0]):
        reads = {}
        for row in segments:
            flags = row[1]
            rid, role = row[2:4]
            if rid not in reads:
                reads[rid] = dict(role=role, declared=row[4:8], primary_count=0,
                    primary=None, mapped=set(), split=set(), supplementary=False,
                    audit=0, last_end=-1, query_length=None)
            read = reads[rid]
            kind, ci, raw, qs, qe, ql, start = row[8:]
            if kind == 0:
                read['primary_count'] += 1
                read['primary'] = (ci, start, raw)
                if not raw & 4:
                    read['mapped'].add(ci)
            if kind != 2 and not raw & 4:
                read['split'].add(ci)
                if qs < read['last_end']:
                    read['audit'] |= QUERY_OVERLAP
                if read['query_length'] is not None and ql != read['query_length']:
                    read['audit'] |= QUERY_LENGTH_CONFLICT
                read['query_length'] = ql
                read['last_end'] = max(read['last_end'], qe)
            read['supplementary'] |= kind == 1
        primary, observed_primary, observed_mates = [], set(), {}
        for read in reads.values():
            role = read['role']
            if read['primary_count'] != 1 or role == 3:
                flags |= INCOMPLETE
                read['audit'] |= INCOMPLETE
            if read['primary_count'] == 1:
                observed_primary.add(role)
                observed_mates[role] = read['primary']
            primary.extend((role, ci) for ci in read['mapped'])
            if read['supplementary']:
                flags |= SPLIT_CROSS if len(read['split']) > 1 else SPLIT_SAME
            flags |= read['audit']
        for rid, read in reads.items():
            role = read['role']
            mate_role = 2 if role == 1 else 1 if role == 2 else None
            actual = observed_mates.get(mate_role)
            if actual is not None:
                ci, start, raw = actual
                declared = read['declared']
                other_unmapped = bool(raw & 4)
                disagree = bool(declared[2]) != other_unmapped
                if not other_unmapped and declared[0] >= 0:
                    disagree |= declared[0] != ci or declared[1] != start or bool(declared[3]) != bool(raw & 16)
                if disagree:
                    flags |= MATE_DISAGREEMENT
                    read['audit'] |= MATE_DISAGREEMENT
            read_updates.append((read['audit'], rid))
        role_set = {read['role'] for read in reads.values()}
        if role_set == {0}:
            flags |= UNPAIRED
        elif role_set != {1, 2}:
            flags |= INCOMPLETE
        if not primary and observed_primary == {1, 2}:
            flags |= BOTH_UNMAPPED
        elif {role for role, _ in primary} == {1, 2}:
            flags |= SAME_PAIR if len({ci for _, ci in primary}) == 1 else CROSS_PAIR
        elif role_set != {0}:
            flags |= SINGLE_MAPPED
        fragment_updates.append((flags, fid))
        if len(read_updates) >= 10000:
            db.executemany('UPDATE reads SET audit_flags=? WHERE id=?', read_updates)
            db.executemany('UPDATE fragments SET flags=? WHERE id=?', fragment_updates)
            read_updates.clear()
            fragment_updates.clear()
    db.executemany('UPDATE reads SET audit_flags=? WHERE id=?', read_updates)
    db.executemany('UPDATE fragments SET flags=? WHERE id=?', fragment_updates)


def capture(alignment, reference, directory, *, namespace='', min_mapq=0,
            exclude_duplicates=False, exclude_qcfail=False, reference_path=None):
    """Stream once into bounded-cache SQLite; persist compact array tables only."""
    import pysam
    from bin.coverage_core import Reference, layout_metadata
    spool_parent, reserve = _spool_settings()
    _check_spool_space(spool_parent, reserve)
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    contigs = {name: i for i, name in enumerate(reference.names)}
    sizes = dict(zip(reference.names, reference.lengths))
    with open(alignment, 'rb') as handle:
        cram = handle.read(4) == b'CRAM'
    # Preserve mate coordinates and AUX. Unused TLEN can trigger external
    # reference loading for multi-reference CRAM slices even without SEQ/QUAL.
    options = dict(format_options=[b'required_fields=6399', b'decode_md=0']) if cram else {}
    if reference_path:
        options['reference_filename'] = str(reference_path)
    started = time.monotonic()
    records = 0
    with tempfile.TemporaryDirectory(prefix='coveragemap-alignment-', dir=spool_parent) as temporary:
        spool_file = Path(temporary)/'spool.sqlite'
        print(f'CAPTURE START alignment={alignment} spool={temporary} reserve_bytes={reserve}', flush=True)
        with closing(sqlite3.connect(str(spool_file))) as db:
            # This private spool is discarded on success/error and is never a
            # recovery artifact. Final arrays/checkpoints retain all validation.
            db.execute('PRAGMA journal_mode=OFF')
            db.execute('PRAGMA synchronous=OFF')
            db.execute('PRAGMA cache_size=-16384')
            db.execute('PRAGMA temp_store=FILE')
            db.executescript('''
                CREATE TABLE strings(id INTEGER PRIMARY KEY,value TEXT UNIQUE NOT NULL);
                CREATE TABLE fragments(id INTEGER PRIMARY KEY,identity TEXT UNIQUE,key TEXT,rg INT,flags INT,read1 INT DEFAULT -1,read2 INT DEFAULT -1);
                CREATE TABLE reads(id INTEGER PRIMARY KEY,fragment INT,role INT,flags INT,rg INT,mate_contig INT,mate_start INT,mate_unmapped INT,mate_reverse INT,audit_flags INT DEFAULT 0,UNIQUE(fragment,role));
                CREATE TABLE segments(id INTEGER PRIMARY KEY,''' + ','.join('"'+n+'" INTEGER' for n in SEGMENT.names) + ''');
                CREATE TABLE blocks(id INTEGER PRIMARY KEY,segment INT,contig INT,start INT,end INT);
                CREATE INDEX segments_read ON segments(read,query_start,query_end);
            ''')
            block_id = 0
            next_fragment, next_read = 0, 0
            intern = _string_interner(db)
            with pysam.AlignmentFile(str(alignment), 'rc' if cram else 'rb', **options) as bam:
                sq = bam.header.to_dict().get('SQ', [])
                names = [r['SN'] for r in sq]
                if len(names) != len(set(names)) or not set(names) <= set(contigs):
                    raise ValueError('Alignment header contigs differ from FAI')
                if reference_path is not None and set(names) != set(contigs):
                    raise ValueError('FASTA mode requires the complete alignment header contig set')
                known = dict(zip(reference.names, reference.md5))
                for row in sq:
                    if row['LN'] != sizes[row['SN']]:
                        raise ValueError('Alignment header length mismatch')
                    md5 = row.get('M5', '')
                    if md5 and known[row['SN']] and md5.lower() != known[row['SN']].lower():
                        raise ValueError('Alignment header M5 mismatch')
                header_md5 = {r['SN']: r.get('M5', '') for r in sq}
                reference = Reference(reference.names, reference.lengths,
                    tuple(known[n] or header_md5.get(n, '') for n in reference.names),
                    tuple(n for n in reference.names if n not in names), reference.sort_policy)
                for sid, read in enumerate(bam.fetch(until_eof=True)):
                    records = sid + 1
                    rg = read.get_tag('RG') if read.has_tag('RG') else ''
                    qname = read.query_name
                    key = json.dumps([namespace, qname, rg, sid if not qname else None], separators=(',', ':'))
                    identity = hashlib.sha256(key.encode()).hexdigest()
                    found = db.execute('SELECT id,key FROM fragments WHERE identity=?', (identity,)).fetchone()
                    if found is None:
                        fid = next_fragment
                        next_fragment += 1
                        db.execute('INSERT INTO fragments(id,identity,key,rg,flags) VALUES (?,?,?,?,?)',
                            (fid, identity, key, intern(rg), (MISSING_RG if not rg else 0) | (INCOMPLETE if not qname else 0)))
                    else:
                        fid, original = found
                        if original != key:
                            raise ValueError('Fragment identity hash collision')
                    role = 1 if read.is_read1 and not read.is_read2 else 2 if read.is_read2 and not read.is_read1 else 0 if not read.is_paired else 3
                    mate = contigs[read.next_reference_name] if read.next_reference_id >= 0 else -1
                    row = db.execute('SELECT id FROM reads WHERE fragment=? AND role=?', (fid, role)).fetchone()
                    if row is None:
                        rid = next_read
                        next_read += 1
                        db.execute('INSERT INTO reads VALUES (?,?,?,?,?,?,?,?,?,?)', (rid, fid, role, read.flag,
                            intern(rg),mate,read.next_reference_start,int(read.mate_is_unmapped),int(read.mate_is_reverse),0))
                        if role in (1,2):
                            db.execute(f'UPDATE fragments SET read{role}=? WHERE id=?',(rid,fid))
                    else:
                        rid = row[0]
                    kind = 2 if read.is_secondary else 1 if read.is_supplementary else 0
                    if read.is_secondary and read.is_supplementary:
                        db.execute('UPDATE fragments SET flags=flags | ? WHERE id=?',(INCOMPLETE,fid))
                    ci = contigs[read.reference_name] if read.reference_id >= 0 else -1
                    if kind == 0:
                        db.execute('UPDATE reads SET flags=?,mate_contig=?,mate_start=?,mate_unmapped=?,mate_reverse=? WHERE id=?',
                            (read.flag,mate,read.next_reference_start,int(read.mate_is_unmapped),int(read.mate_is_reverse),rid))
                    eligible = int(not read.is_unmapped and read.mapping_quality >= min_mapq
                        and not (exclude_duplicates and read.is_duplicate) and not (exclude_qcfail and read.is_qcfail))
                    begin = block_id
                    if not read.is_unmapped:
                        for start, end in read.get_blocks():
                            if ci < 0 or not 0 <= start < end <= reference.lengths[ci]:
                                raise ValueError('Alignment block outside contig')
                            db.execute('INSERT INTO blocks VALUES (?,?,?,?,?)', (block_id, sid, ci, start, end))
                            block_id += 1
                    qs, qe, ql = query_interval(read.cigartuples, read.is_reverse)
                    values = (rid, kind, ci, read.reference_start, read.reference_end if read.reference_end is not None else -1,
                        qs, qe, ql, int(read.is_reverse), read.mapping_quality, read.flag, intern(read.cigarstring or ''),
                        mate, read.next_reference_start, int(read.mate_is_unmapped), int(read.mate_is_reverse),
                        intern(read.get_tag('SA') if read.has_tag('SA') else ''),
                        read.get_tag('AS') if read.has_tag('AS') else -1, read.get_tag('NM') if read.has_tag('NM') else -1,
                        int(read.has_tag('AS')),int(read.has_tag('NM')),eligible, begin, block_id)
                    db.execute('INSERT INTO segments VALUES (' + ','.join('?' for _ in range(1+len(values))) + ')', (sid, *values))
                    if sid and sid % 10000 == 0:
                        db.commit()
                        _check_spool_space(temporary, reserve)
            _audit_relations(db)
            db.commit()
            observed_alignment_counts = {name: db.execute('SELECT count(*) FROM segments WHERE type=?',(number,)).fetchone()[0]
                for number,name in enumerate(('primary','supplementary','secondary'))}
            _export(db, directory)
            print(f'CAPTURE SPOOL_DONE records={records} sqlite_bytes={spool_file.stat().st_size} '
                  f'seconds={time.monotonic()-started:.3f}', flush=True)
    metadata = dict(schema=SCHEMA, layout=layout_metadata(reference),
        missing_header_contigs=list(reference.missing_header_contigs),
        sampling_units=['read', 'fragment'], supplementary_available=True,
        coverage_alignment_policy='primary_only', alignment_types=['primary', 'supplementary', 'secondary'],
        secondary_blocks_available=True, coordinate_convention='0-based half-open; contig local',
        query_coordinate_convention='original-read orientation including hard clips',
        supplementary_overlap_policy='reject overlapping query intervals or inconsistent query length',
        id_policy='first-observation indices + collision-checked namespace/QNAME/RG SHA256', namespace=namespace)
    metadata['audit_flag_bits'] = dict(same_contig_pair=SAME_PAIR,cross_contig_pair=CROSS_PAIR,
        single_mapped_mate=SINGLE_MAPPED,both_unmapped=BOTH_UNMAPPED,unpaired=UNPAIRED,
        incomplete_or_ambiguous=INCOMPLETE,split_read_same_contig=SPLIT_SAME,
        split_read_cross_contig=SPLIT_CROSS,missing_rg=MISSING_RG,
        mate_disagreement=MATE_DISAGREEMENT,query_interval_overlap=QUERY_OVERLAP,
        query_length_conflict=QUERY_LENGTH_CONFLICT)
    metadata['observed_alignment_counts'] = observed_alignment_counts
    (directory/'geometry.json').write_text(json.dumps(metadata, sort_keys=True, indent=2)+'\n')
    return Geometry(directory), reference


class Geometry:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.metadata = json.loads((self.directory/'geometry.json').read_text())
        if self.metadata['schema'] != SCHEMA:
            raise ValueError('Unsupported alignment geometry')
        for name, dtype in (('fragments', FRAGMENT), ('reads', READ), ('segments', SEGMENT), ('blocks', BLOCK), ('strings', STRING)):
            value = np.load(self.directory/(name+'.npy'), mmap_mode='r', allow_pickle=False)
            if value.dtype != dtype or value.ndim != 1:
                raise ValueError('Invalid geometry table '+name)
            setattr(self, name, value)
        if (np.any(self.reads['fragment'] < 0) or np.any(self.reads['fragment'] >= len(self.fragments))
                or np.any(self.segments['read'] < 0) or np.any(self.segments['read'] >= len(self.reads))
                or np.any(self.blocks['segment'] < 0) or np.any(self.blocks['segment'] >= len(self.segments))):
            raise ValueError('Geometry relationship outside table')
        seg = self.segments
        if (np.any(seg['block_begin'] < 0) or np.any(seg['block_end'] < seg['block_begin'])
                or np.any(seg['block_end'] > len(self.blocks))
                or (len(seg) and (seg['block_begin'][0] != 0 or seg['block_end'][-1] != len(self.blocks)
                    or np.any(seg['block_end'][:-1] != seg['block_begin'][1:])))
                or not len(seg) and len(self.blocks)):
            raise ValueError('Invalid segment block pointers')
        for field in ('cigar','sa'):
            if np.any(seg[field] < 0) or np.any(seg[field] >= len(self.strings)):
                raise ValueError('Invalid string dictionary reference')
        text_bytes = (self.directory/'strings.utf8').stat().st_size
        if np.any(self.strings['offset'] < 0) or np.any(self.strings['length'] < 0) or np.any(self.strings['offset']+self.strings['length'] > text_bytes):
            raise ValueError('Invalid UTF8 dictionary bounds')

    def string(self, number):
        offset, length = self.strings[number]
        with (self.directory/'strings.utf8').open('rb') as handle:
            handle.seek(int(offset))
            return handle.read(int(length)).decode('utf-8')

    def project(self, sampling_unit='read', coverage_alignment_policy='primary_only', mate_policy='both'):
        from bin.coverage_core import BlockPool, reference_from_layout
        if sampling_unit not in ('read', 'fragment') or coverage_alignment_policy not in ('primary_only', 'primary_plus_supplementary'):
            raise ValueError('Unsupported runtime geometry capability')
        if mate_policy not in ('both', 'read1', 'read2') or sampling_unit == 'fragment' and mate_policy != 'both':
            raise ValueError('Invalid mate policy for sampling unit')
        ref = reference_from_layout(self.metadata['layout'])
        seg = self.segments
        eligible = (seg['eligible'] == 1) & (seg['block_end'] > seg['block_begin'])
        eligible &= seg['type'] == 0 if coverage_alignment_policy == 'primary_only' else seg['type'] != 2
        if mate_policy != 'both':
            roles = self.reads['role'][seg['read']]
            eligible &= (roles == 0) | (roles == (1 if mate_policy == 'read1' else 2))
        numbers = np.flatnonzero(eligible)
        read_ids = seg['read'][numbers]
        primary_ids = seg['read'][(seg['type']==0)&(seg['eligible']==1)&(seg['block_end']>seg['block_begin'])]
        if np.any(np.bincount(primary_ids,minlength=len(self.reads))>1):
            raise ValueError('Ambiguous read identity has multiple primary alignments; captured evidence retained')
        if coverage_alignment_policy == 'primary_plus_supplementary':
            ordered = numbers[np.argsort(read_ids, kind='stable')]
            # Sort each read's query intervals; one read only held as numeric indices.
            boundaries = np.r_[0, np.flatnonzero(np.diff(seg['read'][ordered]))+1, len(ordered)]
            for left, right in zip(boundaries, boundaries[1:]):
                group = ordered[left:right]
                if not len(group):
                    continue
                if np.sum(seg['type'][group] == 0) != 1 or len(set(seg['query_length'][group])) != 1:
                    raise ValueError('Incomplete/ambiguous split read or inconsistent query length')
                interval_order = group[np.argsort(seg['query_start'][group],kind='stable')]
                lefts,rights = seg['query_start'][interval_order],seg['query_end'][interval_order]
                if np.any(lefts<0) or np.any(rights<=lefts) or np.any(lefts[1:]<np.maximum.accumulate(rights)[:-1]):
                    raise ValueError('Overlapping or invalid split-read query intervals')
        units = read_ids if sampling_unit == 'read' else self.reads['fragment'][read_ids]
        order = np.argsort(units, kind='stable')
        starts, ends, pointers, identities = array('Q'), array('Q'), array('Q', [0]), array('Q')
        last = None
        offsets = list(ref.offsets.values())
        for position in order:
            unit, sid = int(units[position]), int(numbers[position])
            if last is not None and unit != last:
                pointers.append(len(starts))
            if unit != last:
                identities.append(unit)
            last = unit
            row = seg[sid]
            for block in self.blocks[int(row['block_begin']):int(row['block_end'])]:
                ci, start, end = int(block['contig']), int(block['start']), int(block['end'])
                if int(block['segment']) != sid or not 0 <= ci < len(ref.names) or not 0 <= start < end <= ref.lengths[ci]:
                    raise ValueError('Invalid segment block relation or bounds')
                starts.append(offsets[ci]+start)
                ends.append(offsets[ci]+end)
        if last is not None:
            pointers.append(len(starts))
        pool = BlockPool(np.asarray(starts), np.asarray(ends), np.asarray(pointers), ref.length)
        pool.alignment_geometry = self
        pool.geometry_unit_ids = np.asarray(identities)
        pool.sampling_unit = sampling_unit
        pool.coverage_alignment_policy = coverage_alignment_policy
        return pool

    @property
    def nbytes(self):
        return sum(getattr(self,name).nbytes for name in ('fragments','reads','segments','blocks','strings')) + (self.directory/'strings.utf8').stat().st_size


def selection_audit(pool, ids, target_bases):
    ids = np.asarray(ids,dtype=np.intp)
    actual = int(pool.unit_bases[ids].sum())
    unit = getattr(pool,'sampling_unit','read')
    output = dict(sampling_unit=unit,coverage_alignment_policy=getattr(pool,'coverage_alignment_policy','primary_only'),
        target_bases=float(target_bases),actual_bases=actual,target_depth=float(target_bases)/pool.genome_len,
        actual_depth=actual/pool.genome_len,overshoot_bases=actual-float(target_bases),
        selected_units=len(ids),selected_blocks=int((pool.indptr[ids+1]-pool.indptr[ids]).sum()))
    geometry = getattr(pool,'alignment_geometry',None)
    if geometry is None:
        output.update(selected_reads=len(ids) if unit=='read' else None,
                      selected_fragments=len(ids) if unit=='fragment' else None,selected_segments=None)
        return output
    chosen = pool.geometry_unit_ids[ids]
    read_ids = chosen if unit=='read' else np.flatnonzero(np.isin(geometry.reads['fragment'],chosen))
    segments = np.isin(geometry.segments['read'],read_ids)
    segments &= geometry.segments['eligible']==1
    segments &= geometry.segments['type']==0 if output['coverage_alignment_policy']=='primary_only' else geometry.segments['type']!=2
    segments &= geometry.segments['block_end']>geometry.segments['block_begin']
    output.update(selected_reads=len(read_ids),selected_fragments=len(np.unique(geometry.reads['fragment'][read_ids])),
                  selected_segments=int(segments.sum()))
    return output
