"""Integrity checks prevent source corruption and fabricated feasibility passes."""
import importlib.util
from pathlib import Path
import unittest


spec = importlib.util.spec_from_file_location(
    'small_audit', Path(__file__).resolve().parents[1] / 'src/dataset_search/audit_small_files.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class SourceIntegrityTests(unittest.TestCase):
    def test_git_blob_includes_header_and_exact_bytes(self):
        # Canonical empty git blob: content SHA1 alone is different.
        self.assertEqual(module.verify_blob('', 'e69de29bb2d1d6434b8b29ae775ad8c2e48c5391'), b'')

    def test_corrupted_source_rejected(self):
        with self.assertRaises(ValueError):
            module.verify_blob('changed', 'e69de29bb2d1d6434b8b29ae775ad8c2e48c5391')

    def test_notebook_code_is_separate_from_saved_outputs(self):
        text = '{"cells":[{"source":["x=1"],"outputs":[{"text":"fake cells"}]}]}'
        self.assertEqual(module.notebook_source(text), 'x=1')


if __name__ == '__main__':
    unittest.main()
