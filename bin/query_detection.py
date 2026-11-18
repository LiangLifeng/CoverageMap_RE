"""Genome-length rules and normalized coverage distances for manifest detection."""

import json

import numpy as np

from bin.bin_count_eval_math import METRICS, normalized

DOMAINS = ('bacteria', 'viral')
RULE_KEYS = {'min_genome_bp', 'max_genome_bp_exclusive', 'bin_count', 'distance', 'max_distance'}


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f'Duplicate detection config key: {key}')
        result[key] = value
    return result


def read_rules(path, domain):
    document = json.loads(path.read_text(encoding='utf-8'), object_pairs_hook=unique_object)
    if not isinstance(document, dict) or not document or set(document) - set(DOMAINS):
        raise ValueError('Detection config must have bacteria and/or viral as its top-level keys')
    for name, rules in document.items():
        if not isinstance(rules, list) or not rules:
            raise ValueError(f'{name} must contain a nonempty list of genome-length rules')
        for number, rule in enumerate(rules, 1):
            where = f'{name} rule {number}'
            if not isinstance(rule, dict) or set(rule) != RULE_KEYS:
                raise ValueError(f'{where} requires exactly: {", ".join(sorted(RULE_KEYS))}')
            lower, upper = rule['min_genome_bp'], rule['max_genome_bp_exclusive']
            if type(lower) is not int or lower < 0 or (upper is not None and (type(upper) is not int or upper <= lower)):
                raise ValueError(f'{where} requires a valid half-open genome length interval in bp')
            if rule['distance'] not in METRICS:
                raise ValueError(f'{where} distance must be one of {", ".join(METRICS)}')
            bins = rule['bin_count']
            minimum_bins = 2 if rule['distance'] == 'correlation' else 1
            if type(bins) is not int or bins < minimum_bins or (upper is not None and bins >= upper):
                raise ValueError(f'{where} requires bin_count >= {minimum_bins} and below the upper bp bound')
            threshold = rule['max_distance']
            if threshold is not None and (type(threshold) not in (int, float) or not 0 <= threshold < float('inf')):
                raise ValueError(f'{where} max_distance must be null or a finite nonnegative number')
        rules.sort(key=lambda rule: rule['min_genome_bp'])
        for previous, current in zip(rules, rules[1:]):
            upper = previous['max_genome_bp_exclusive']
            if upper is None or upper > current['min_genome_bp']:
                raise ValueError(f'{name} detection rules have overlapping genome length intervals')
    if domain not in document:
        raise ValueError(f'Detection config has no rules for {domain}')
    return document[domain]


def select_rule(rules, genome_bp):
    for rule in rules:
        upper = rule['max_genome_bp_exclusive']
        if rule['min_genome_bp'] <= genome_bp and (upper is None or genome_bp < upper):
            return rule
    return None


def coverage_profile(mass, widths):
    return normalized(mass / widths)


def reference_center(profiles, count, bins, work_bytes):
    """The grid evaluator's arithmetic mean and 16 MiB reduction blocks."""
    row_bytes = 8 * bins
    block_size = min(count, max(1, (16 * 1024**2) // row_bytes))
    required = (block_size + 5) * row_bytes
    if work_bytes is not None and required > work_bytes:
        raise ValueError(f'Detection profiles exceed work budget: {required} bytes')
    block = np.empty((block_size, bins), dtype=np.float64)
    total = np.zeros(bins, dtype=np.float64)
    for number, profile in enumerate(profiles, 1):
        block[(number - 1) % block_size] = profile
        if number % block_size == 0:
            total += block.sum(axis=0)
    if count % block_size:
        total += block[:count % block_size].sum(axis=0)
    return total / count
