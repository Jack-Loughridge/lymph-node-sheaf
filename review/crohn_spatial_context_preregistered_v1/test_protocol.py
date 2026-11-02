import unittest
from protocol_checks import paired_folds,check_features,verify_frozen
from background_audit import per_control_background
from retrieve_supporting_files import validate
class ProtocolTests(unittest.TestCase):
 def rows(self):return [dict(patient_id=str(i),specimen_id=f'T{i}{s}',status=s,S1_verified=True) for i in range(8) for s in ['affected','adjacent']]
 def test_pair_holdout(self):
  rows=self.rows();folds=paired_folds(rows);self.assertEqual(len(folds),8)
  for tr,te in folds:
   self.assertEqual(len(tr),14);self.assertEqual(len(te),2);self.assertFalse(set(tr)&set(te))
   self.assertEqual({v['patient_id'] for v in rows if v['specimen_id'] in te}.__len__(),1)
 def test_no_slide_surrogate(self):
  rows=self.rows()
  for r in rows:
   if r['patient_id']=='7':r['patient_id']='6'
  with self.assertRaises(ValueError):paired_folds(rows)
 def test_requires_s1(self):
  rows=self.rows();rows[0]['S1_verified']=False
  with self.assertRaises(ValueError):paired_folds(rows)
 def test_reject_duplicate(self):
  rows=self.rows();rows[0]['specimen_id']=rows[1]['specimen_id']
  with self.assertRaises(ValueError):paired_folds(rows)
 def test_complete_pairs(self):
  rows=self.rows();rows[0]['status']='adjacent'
  with self.assertRaises(ValueError):paired_folds(rows)
 def test_reuse_rejected(self):
  for prohibited in ['Y','MARKER','LIGAND']:
   with self.assertRaises(ValueError):check_features([prohibited],['Y'],['MARKER'],['LIGAND'])
  self.assertTrue(check_features(['REFERENCE'],['Y'],['MARKER'],['LIGAND']))
 def test_probe_denominator(self):
  self.assertEqual(per_control_background(40,100,20),.02)
  with self.assertRaises(ValueError):per_control_background(40,100,0)
 def test_html_not_csv(self):
  with self.assertRaises(ValueError):validate('data.csv',b'<html>download failed</html>')
 def test_original_and_protocol_frozen(self):self.assertTrue(verify_frozen())
if __name__=='__main__':unittest.main()
