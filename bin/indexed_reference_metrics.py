"""Use compatible, unchanged SQL scalar capacities before decoding references."""
import json
from pathlib import Path

FIELDS = ('aligned_bases','read_count','fragment_count','read1_bases','read2_bases',
          'read1_count','read2_count','cram_size','cram_mtime_ns',
          'fna_size','fna_mtime_ns','fai_size','fai_mtime_ns')


def bind_metrics(connection, references, config):
    """Attach only small scalar metrics; historical indexes retain decode fallback."""
    if not references:
        return
    metadata = {r[0]:json.loads(r[1]) for r in connection.execute(
        'SELECT key,value FROM metadata WHERE key IN '
        "('metrics_schema','min_mapq','m5_policy','coverage_alignment_policy')")}
    if (metadata.get('metrics_schema') != 'coverage-reference-metrics-1'
            or metadata.get('min_mapq') != config.get('min_mapq',0)
            or metadata.get('m5_policy') != config.get('m5_policy','strict')
            or metadata.get('coverage_alignment_policy') != config.get('coverage_alignment_policy','primary_only')):
        return
    selected = {r['id']:r for r in references}
    fields = ','.join('m.'+field for field in FIELDS)
    for row in connection.execute('SELECT m.id,'+fields+' FROM reference_metrics m '
            'JOIN reference_records r ON m.id=r.id WHERE r.GCFID=?',(references[0]['GCFID'],)):
        if row[0] in selected:
            selected[row[0]]['sql_capacity'] = dict(zip(FIELDS,row[1:]))


def capacities(reference, arms, spec):
    """Return None on stale source/reference metadata, forcing actual extraction."""
    metric = reference.get('sql_capacity')
    if not metric:
        return None
    paths = [('cram',reference['CramFilePath']),('fna',spec.get('fna')),('fai',spec.get('fai'))]
    try:
        for name,path in paths:
            if not path:
                return None
            stat = Path(path).stat()
            if (stat.st_size,stat.st_mtime_ns) != (metric[name+'_size'],metric[name+'_mtime_ns']):
                return None
    except OSError:
        return None
    result = {}
    for arm in arms:
        if arm.mate_side in ('read1','read2'):
            result[arm.name] = (metric[arm.mate_side+'_bases'],metric[arm.mate_side+'_count'])
        else:
            count = metric['fragment_count'] if spec['config']['sampling_unit']=='fragment' else metric['read_count']
            result[arm.name] = (metric['aligned_bases'],count)
    return result
