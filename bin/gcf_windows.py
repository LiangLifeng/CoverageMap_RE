"""Versioned windows over a GCF's virtual coordinate system."""
from itertools import chain

import numpy as np

from bin.coverage_core import CONTIG_SORT_POLICY, window_edges

GCF_CONCAT = 'gcf-concat-v1'
PER_CONTIG = 'per-contig-fixed-bp-terminal-short-v1'


def window_policy(metadata):
    policy = metadata.get('window_policy', PER_CONTIG)
    if policy not in (GCF_CONCAT, PER_CONTIG):
        raise ValueError(f'Unsupported window policy: {policy}')
    return policy


def require_layout_policy(layout, policy):
    if policy not in (GCF_CONCAT, PER_CONTIG):
        raise ValueError(f'Unsupported window policy: {policy}')
    expected = CONTIG_SORT_POLICY if policy == GCF_CONCAT else None
    if layout.sort_policy != expected:
        raise ValueError(f'Contig layout sort policy does not match {policy}')


def concat_edges(length, *, size=None, count=None):
    if (size is None) == (count is None):
        raise ValueError('Specify exactly one window-size or bin-count')
    if count is not None:
        if type(count) is not int:
            raise ValueError('bin-count must be an integer')
        return window_edges(length, count)
    if type(size) is not int or size < 1:
        raise ValueError('window-size must be a positive integer')
    if not 1 <= length <= np.iinfo(np.int64).max:
        raise ValueError('Genome length must be positive and fit signed 64-bit coordinates')
    return np.fromiter(chain(range(0, length, size), (length,)), dtype=np.int64,
                       count=(length-1)//size+2)


def per_contig_edges(layout, size):
    if type(size) is not int or size < 1:
        raise ValueError('window-size must be a positive integer')
    values = [0]
    for length in layout.lengths:
        offset = values[-1]
        values.extend(offset+end for end in range(size, length, size))
        values.append(offset+length)
    return np.asarray(values, dtype=np.int64)


def edges_for(layout, policy, *, size=None, count=None):
    require_layout_policy(layout, policy)
    if policy == GCF_CONCAT:
        return concat_edges(layout.length, size=size, count=count)
    if (size is None) == (count is None) or (size or count) < 1:
        raise ValueError('Specify exactly one positive window-size or bin-count')
    width = size if size is not None else (layout.length+count-1)//count
    return per_contig_edges(layout, width)


def segments(layout, edges):
    """CSR-style window-to-segment mapping, ordered by global coordinates."""
    offsets = np.asarray(list(layout.offsets.values()), dtype=np.int64)
    cuts = np.union1d(edges, offsets)
    contigs = np.searchsorted(offsets, cuts[:-1], side='right') - 1
    return dict(segment_ptr=np.searchsorted(cuts, edges),
                segment_contig_index=contigs,
                segment_local_start=cuts[:-1]-offsets[contigs],
                segment_local_end=cuts[1:]-offsets[contigs])


def coordinate_rows(layout, edges, nominal_width):
    mapping = segments(layout, edges)
    pointers = mapping['segment_ptr']
    for i, (start, end) in enumerate(zip(edges[:-1], edges[1:])):
        for j in range(pointers[i], pointers[i+1]):
            yield dict(bin_index=i, global_start=start, global_end=end,
                       width=end-start, nominal_width=nominal_width,
                       contig=layout.names[mapping['segment_contig_index'][j]],
                       local_start=mapping['segment_local_start'][j],
                       local_end=mapping['segment_local_end'][j])


COORDINATE_COLUMNS = ('bin_index', 'global_start', 'global_end', 'width', 'nominal_width',
                      'contig', 'local_start', 'local_end')
