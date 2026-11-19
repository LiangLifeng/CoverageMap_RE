"""Deterministic source replacement using actual eligible whole-read capacity."""


class InsufficientReferenceCapacity(ValueError):
    def __init__(self, message, events):
        super().__init__(message)
        self.events = events


def resolve_reference_choices(spec, arms, capacity, make_rng):
    """Keep every slot and seed; replace only unusable source capacity.

    capacity(reference, arm) returns (aligned_bases, eligible_unit_count).
    It must validate inputs rather than treating decode errors as low capacity.
    """
    cfg = spec['config']
    refs = {r['id']: r for r in spec['references']}
    choices = list(spec['choices'])
    slots = [arm for arm in arms for _ in range(arm.slots)]
    if len(slots) != len(choices):
        raise ValueError('Arm slots must match the original source choices')
    tolerance = cfg.get('sampling_tolerance_bp', 150)
    verified, banned, events = {}, {arm.name: set() for arm in arms}, []

    def inspect(refid, arm):
        key = refid, arm.name
        if key not in verified:
            verified[key] = capacity(refs[refid], arm)
        return verified[key]

    def eligible(refid, arm):
        mass, units = inspect(refid, arm)
        return units > 0 and mass + tolerance >= arm.target_bases

    for slot, arm in enumerate(slots):
        refid, attempt = choices[slot], 0
        while not eligible(refid, arm):
            mass, units = inspect(refid, arm)
            banned[arm.name].add(refid)
            available = [r for r in refs if r not in banned[arm.name]]
            if not available:
                events.append(dict(draw_id=slot, mate_side=arm.mate_side, attempt=attempt+1,
                    rejected_reference_id=refid, replacement_reference_id=None,
                    reason='NO_ELIGIBLE_REFERENCE_REMAINS', available_bases=mass,
                    target_bases=arm.target_bases, tolerance_bp=tolerance, rng_sha256=''))
                raise InsufficientReferenceCapacity(

                    f'No compatible reference can supply {arm.name} target within tolerance', events)
            generator, token = make_rng(cfg['seed'], spec['query_id'],
                                        spec['query']['Sample'], spec['query']['GCFID'],
                                        'replacement-reference', slot, attempt)
            replacement = available[int(generator.integers(0, len(available)))]
            events.append(dict(draw_id=slot, mate_side=arm.mate_side,
                               attempt=attempt+1, rejected_reference_id=refid,
                               replacement_reference_id=replacement,
                               reason='NO_COMPATIBLE_UNITS' if not units else 'INSUFFICIENT_CAPACITY',
                               available_bases=mass, target_bases=arm.target_bases,
                               tolerance_bp=tolerance, rng_sha256=token))
            refid, attempt = replacement, attempt+1
        choices[slot] = refid

    # Prove distinct-sample support after actual capacity exclusions. Probe
    # only as many previously unseen sources as needed, without changing draws.
    for arm in arms:
        samples = {refs[r]['Sample'] for r, name in verified
                   if name == arm.name and eligible(r, arm)}
        for refid in refs:
            if len(samples) > cfg['min_samples']:
                break
            if (refid, arm.name) in verified:
                continue
            if eligible(refid, arm):
                samples.add(refs[refid]['Sample'])
        if len(samples) <= cfg['min_samples']:
            raise InsufficientReferenceCapacity(
                f'Actual eligible distinct reference samples for {arm.name} '
                f'must exceed min_samples={cfg["min_samples"]}', events)
    return choices, events
