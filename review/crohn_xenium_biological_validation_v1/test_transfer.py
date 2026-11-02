import tempfile,unittest
from pathlib import Path
import verify_transfer as v

class TransferTests(unittest.TestCase):
 def fixture(self,root):
  (root/'sample').mkdir();p=root/'sample/matrix.bin';p.write_bytes(b'transfer test')
  return {'files':[{'relative_path':'sample/matrix.bin','bytes':p.stat().st_size,'sha256':v.sha256(p)}]}
 def test_matching_payload(self):
  with tempfile.TemporaryDirectory() as td:
   root=Path(td);m=self.fixture(root);self.assertTrue(v.verify(root,m)[0]['verified'])
 def test_corruption(self):
  with tempfile.TemporaryDirectory() as td:
   root=Path(td);m=self.fixture(root);(root/'sample/matrix.bin').write_bytes(b'transfer TEST');self.assertFalse(v.verify(root,m)[0]['verified'])
 def test_missing(self):
  with tempfile.TemporaryDirectory() as td:
   root=Path(td);m=self.fixture(root);(root/'sample/matrix.bin').unlink();self.assertEqual(v.verify(root,m)[0]['failure'],'missing_file')
 def test_size(self):
  with tempfile.TemporaryDirectory() as td:
   root=Path(td);m=self.fixture(root);m['files'][0]['bytes']+=1;self.assertEqual(v.verify(root,m)[0]['failure'],'byte_count_mismatch')
 def test_traversal(self):
  with tempfile.TemporaryDirectory() as td:
   with self.assertRaises(ValueError):v.resolve_entry(Path(td),'../secret')
 def test_absolute_path(self):
  with tempfile.TemporaryDirectory() as td:
   with self.assertRaises(ValueError):v.resolve_entry(Path(td),'/tmp/secret')
 def test_windows_path(self):
  with tempfile.TemporaryDirectory() as td:
   with self.assertRaises(ValueError):v.resolve_entry(Path(td),'sample\\matrix.bin')
 def test_duplicate(self):
  with tempfile.TemporaryDirectory() as td:
   root=Path(td);m=self.fixture(root);m['files']*=2
   with self.assertRaises(ValueError):v.verify(root,m)
 def test_no_source_modification(self):
  with tempfile.TemporaryDirectory() as td:
   root=Path(td);m=self.fixture(root);before=(root/'sample/matrix.bin').read_bytes();v.verify(root,m);self.assertEqual((root/'sample/matrix.bin').read_bytes(),before)
if __name__=='__main__':unittest.main()
