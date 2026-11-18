#!/usr/bin/env python3
"""Run manifest queries against a per-query configurable-reference similarity distribution."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).absolute().parents[1]))
from bin.query_manifest import main


if __name__ == '__main__':
    main(reference_distribution=True)
