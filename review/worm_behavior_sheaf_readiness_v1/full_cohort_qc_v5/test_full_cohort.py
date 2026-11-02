import unittest,zipfile,stat
import numpy as np
from extract_processed import select_members
from audit_cohort import clock_matches,plateau_lengths
class FullCohortTests(unittest.TestCase):
 def info(self,name):
  z=zipfile.ZipInfo(name);z.file_size=128;z.external_attr=(stat.S_IFREG|0o644)<<16;return z
 def test_select_only_allowlist(self):
  infos=[self.info('folder/one.h5'),self.info('README.md')]
  self.assertEqual(set(select_members(infos,{'one.h5'})),{'one.h5'})
 def test_duplicates_fail(self):
  self.assertRaises(ValueError,select_members,[self.info('one.h5'),self.info('other/one.h5')],{'one.h5'})
 def test_missing_fails(self):
  self.assertRaises(ValueError,select_members,[],{'one.h5'})
 def test_path_traversal_absolute_windows_fail(self):
  for name in ['../bad','/bad','C:/bad','a\\bad']:
   self.assertRaises(ValueError,select_members,[self.info(name)],set())
 def test_symlink_fails(self):
  z=self.info('one.h5');z.external_attr=(stat.S_IFLNK|0o777)<<16
  self.assertRaises(ValueError,select_members,[z],{'one.h5'})
 def test_size_and_encryption_fail(self):
  z=self.info('one.h5');z.file_size=20_000_001
  self.assertRaises(ValueError,select_members,[z],{'one.h5'})
  z=self.info('one.h5');z.flag_bits=1
  self.assertRaises(ValueError,select_members,[z],{'one.h5'})
 def test_clock_joint_values_required(self):
  tn=np.array([0.,.1,.2,.3]);tc=np.array([.05,.25]);sx=np.array([1.,2.,4.,8.]);sy=sx*3
  match=clock_matches(tc,tn,sx,sy,np.array([1.5,6.]),np.array([4.5,18.]))
  self.assertEqual(match,{'frames':2,'unique':2,'ambiguous':0,'unmatched':0})
  bad=clock_matches(tc,tn,sx,sy,np.array([99.,6.]),np.array([4.5,18.]))
  self.assertEqual(bad['unmatched'],1)
 def test_plateaus_not_saturation_claims(self):
  self.assertEqual(plateau_lengths(np.array([[1,1,2,3,3,3],[0,1,2,3,4,5]])),[3,1])
if __name__=='__main__':unittest.main()
