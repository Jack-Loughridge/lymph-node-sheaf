import unittest,tempfile,json,hashlib
from pathlib import Path
from check_transfer import inspect
class TransferTests(unittest.TestCase):
 def test_original_five_receipt_cannot_pass_remaining_cohort(self):
  with tempfile.TemporaryDirectory() as t:
   p=Path(t);f=p/'sample-cudata.h5';f.write_bytes(b'fixture')
   receipt=p/'receipt.json';receipt.write_text(json.dumps({'algorithm':'SHA256','files':[{'filename':f.name,'bytes':7,'sha256':hashlib.sha256(b'fixture').hexdigest()}]}))
   result=inspect(p,receipt,{'remaining-cudata.h5'},{f.name,'remaining-cudata.h5'})
   self.assertEqual(result['receipt_entries_verified'],1)
   self.assertFalse(result['receipt_matches_exact_remaining24_manifest'])
   self.assertNotEqual(result['status'],'PASS')
 def test_full_transfer_and_hash_failure(self):
  with tempfile.TemporaryDirectory() as t:
   p=Path(t);names={f'{i}-cudata.h5' for i in range(24)};rows=[]
   for n in names:
    (p/n).write_bytes(b'x');rows.append({'filename':n,'bytes':1,'sha256':hashlib.sha256(b'x').hexdigest()})
   receipt=p/'receipt.json';receipt.write_text(json.dumps({'algorithm':'SHA256','files':rows}))
   self.assertEqual(inspect(p,receipt,names,names)['status'],'PASS')
   (p/next(iter(names))).write_bytes(b'y')
   self.assertEqual(inspect(p,receipt,names,names)['receipt_hash_or_size_failures'],1)
if __name__=='__main__':unittest.main()
