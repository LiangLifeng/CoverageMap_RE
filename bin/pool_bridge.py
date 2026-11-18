"""Explicit CRAM extraction bridge between incompatible Python ABIs."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from bin.coverage_core import BlockPool, extract_pool, layout_metadata, reference_from_layout


def extract_with_python(python, alignment, reference, reference_index, policies, temporary_root):
    with tempfile.TemporaryDirectory(prefix='extraction-',dir=temporary_root) as temporary:
        directory=Path(temporary)
        request=dict(alignment=str(alignment),reference=str(reference) if reference else None,
                     reference_index=str(reference_index) if reference_index else None,policies=policies)
        (directory/'request.json').write_text(json.dumps(request))
        command=[str(python),'-B',str(Path(__file__).resolve()),str(directory)]
        subprocess.run(command,check=True)
        metadata=json.loads((directory/'metadata.json').read_text())
        pool=BlockPool(*(np.load(directory/f'{name}.npy',allow_pickle=False) for name in ['starts','ends','indptr']),metadata['genome_len'])
        layout=reference_from_layout(metadata, metadata.get('missing_header_contigs', []))
        metadata['temporary_file_bytes']=sum(p.stat().st_size for p in directory.iterdir())
        metadata['command']=command
    return pool,layout,metadata


def main():
    parser=argparse.ArgumentParser();parser.add_argument('directory',type=Path);a=parser.parse_args()
    request=json.loads((a.directory/'request.json').read_text())
    pool,layout=extract_pool(request['alignment'],request['reference'],reference_index=request['reference_index'],**request['policies'])
    for name in ['starts','ends','indptr']:
        np.save(a.directory/f'{name}.npy',getattr(pool,name),allow_pickle=False)
    import pysam
    (a.directory/'metadata.json').write_text(json.dumps(dict(
        layout_metadata(layout), missing_header_contigs=layout.missing_header_contigs,
        genome_len=layout.length, python=sys.version, numpy=np.__version__, pysam=pysam.__version__)))


if __name__=='__main__':
    main()
