# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 openEMS contributors
"""File-based requests for an opt-in running solver (``--dft-snapshots``).

The solver must use this simulation directory as its working directory. Only
one outstanding request is allowed. The returned directory is usable after
``snapshot.xml`` appears there; a sibling ``<id>.error`` reports failure.
These are incomplete accumulated fields, not a solver restart checkpoint.
"""
import os
import re
import tempfile
from pathlib import Path


def request_dft_snapshot(sim_path, request_id):
    """Publish a request atomically, without overwriting a pending request.

    Returns the future immutable snapshot directory. Requires a local filesystem
    supporting hard links (also supported by NTFS). No network service is opened.
    """
    if (not isinstance(request_id, str) or
            not re.fullmatch(r'[A-Za-z0-9_-]{1,64}', request_id)):
        raise ValueError('request_id must contain 1-64 ASCII letters/digits/-/_')
    base = Path(sim_path)
    target = base / 'dft-snapshots' / request_id
    if target.exists() or target.with_name(request_id + '.error').exists():
        raise FileExistsError('request ID has already been used')
    if (base / 'DFT_SNAPSHOT.active').exists():
        raise FileExistsError('a DFT snapshot is being written')
    fd, temporary = tempfile.mkstemp(prefix='.dft-request-', dir=str(base))
    try:
        with os.fdopen(fd, 'w') as stream:
            stream.write(request_id + '\n')
        os.link(temporary, str(base / 'DFT_SNAPSHOT'))
    finally:
        os.unlink(temporary)
    return target
