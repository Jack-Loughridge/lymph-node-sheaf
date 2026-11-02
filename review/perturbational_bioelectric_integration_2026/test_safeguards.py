import unittest,tempfile,zipfile,stat,hashlib
from pathlib import Path
from io import BytesIO
from unittest.mock import patch
import h5py,numpy as np
import cleanup_verified as cleanup
from audit_small_payload import manifest,xlsx_schema,audit,member_schema
from audit_author_exporter import literal_exports
class SafeguardsTests(unittest.TestCase):
 def info(self,name):
  i=zipfile.ZipInfo(name);i.file_size=10;i.external_attr=(stat.S_IFREG|0o644)<<16;return i
 def test_cleanup_exact_regular_path_and_protected_names(self):
  with tempfile.TemporaryDirectory() as d:
   base=Path(d);f=base/'source.dat';f.write_bytes(b'x')
   with patch.object(cleanup,'BASE',base):
    self.assertEqual(cleanup.checked_path(f),f)
    forbidden=base/'humann_databases';forbidden.mkdir();(forbidden/'x').write_bytes(b'x')
    self.assertRaises(ValueError,cleanup.checked_path,forbidden/'x')
    symlink=base/'link';symlink.symlink_to(f);self.assertRaises(ValueError,cleanup.checked_path,symlink)
 def test_cleanup_outside_root_refused(self):
  with tempfile.TemporaryDirectory() as d:
   f=Path(d)/'file';f.write_bytes(b'x')
   with patch.object(cleanup,'BASE',Path(d)/'different'):self.assertRaises(ValueError,cleanup.checked_path,f)
 def test_unsafe_zip_paths(self):
  for name in ['../x','/x','C:/x','a\\x']:
   b=BytesIO()
   with zipfile.ZipFile(b,'w') as z:z.writestr(name,'x')
   with zipfile.ZipFile(BytesIO(b.getvalue())) as z:self.assertRaises(ValueError,manifest,z)
 def test_zip_symlink_refused(self):
  b=BytesIO();i=self.info('link');i.external_attr=(stat.S_IFLNK|0o777)<<16
  with zipfile.ZipFile(b,'w') as z:z.writestr(i,'target')
  with zipfile.ZipFile(BytesIO(b.getvalue())) as z:self.assertRaises(ValueError,manifest,z)
 def test_xlsx_is_schema_not_GO(self):
  b=BytesIO();xml='<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData><row r="1"><c t="inlineStr"><is><t>Time (s)</t></is></c></row><row r="2"><c><v>0.0</v></c></row><row r="3"><c><v>0.1</v></c></row></sheetData></worksheet>'
  with zipfile.ZipFile(b,'w') as z:z.writestr('xl/worksheets/sheet1.xml',xml)
  result=xlsx_schema(b.getvalue());self.assertEqual(result[0]['row_count'],3);self.assertTrue(result[0]['potential_time_labels']);self.assertTrue(result[0]['not_proof_of_simultaneous_cells'])
 def test_hdf_metadata_not_biological_validation(self):
  b=BytesIO()
  with h5py.File(b,'w') as f:f['traces']=np.ones((20,30));f['time']=np.arange(30)
  result=member_schema('sample.h5',b.getvalue());self.assertEqual(result['datasets']['traces']['shape'],[20,30]);self.assertFalse(result['physiological_semantics_verified'])
 def test_hash_mismatch_and_no_automatic_GO(self):
  with tempfile.TemporaryDirectory() as d:
   f=Path(d)/'sample.csv';f.write_text('time,cell1\n0,1\n1,2\n')
   self.assertRaises(ValueError,audit,f,'0'*64)
   result=audit(f,hashlib.sha256(f.read_bytes()).hexdigest());self.assertTrue(result['sender_hash_checked']);self.assertFalse(result['biological_GO'])
 def test_comments_do_not_create_exported_fields(self):
  text='% ret.trace = ignored;\nret.eventFrames = frames;\nret.(key)=value;\n'
  self.assertEqual(literal_exports(text),['eventFrames'])
 def test_no_physiology_inference_from_unknown_binary(self):
  result=member_schema('file.ibw',b'not decoded');self.assertTrue(result['reader_required_before_GO'])
if __name__=='__main__':unittest.main()
