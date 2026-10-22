# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 openEMS contributors
import tempfile
import unittest
from pathlib import Path
from openEMS.snapshots import request_dft_snapshot


class TestDFTSnapshotRequest(unittest.TestCase):
    def test_atomic_request_and_pending_rejection(self):
        with tempfile.TemporaryDirectory() as name:
            base = Path(name)
            expected = request_dft_snapshot(base, 'preview-1')
            self.assertEqual(expected, base / 'dft-snapshots/preview-1')
            self.assertEqual((base / 'DFT_SNAPSHOT').read_text(), 'preview-1\n')
            with self.assertRaises(FileExistsError):
                request_dft_snapshot(base, 'preview-2')
            self.assertEqual((base / 'DFT_SNAPSHOT').read_text(), 'preview-1\n')
            self.assertEqual(len(list(base.glob('.dft-request-*'))), 0)

    def test_invalid_and_reused_ids(self):
        with tempfile.TemporaryDirectory() as name:
            base = Path(name)
            for bad in ['', '../escape', 'a\nb', 'x' * 65, 'é', '.hidden']:
                with self.assertRaises(ValueError):
                    request_dft_snapshot(base, bad)
            (base / 'dft-snapshots/old').mkdir(parents=True)
            with self.assertRaises(FileExistsError):
                request_dft_snapshot(base, 'old')
            (base / 'DFT_SNAPSHOT.active').write_text('active\n')
            with self.assertRaises(FileExistsError):
                request_dft_snapshot(base, 'next')
