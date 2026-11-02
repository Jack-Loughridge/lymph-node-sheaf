import unittest,tempfile,os
from pathlib import Path
import numpy as np
from qc_primitives import fixed_split,valid_windows,common_core_exists,identity_join
from cleanup_exact_downloads import check_path,TARGETS
from acquire_sample import verify_format
class ReadinessTests(unittest.TestCase):
 def test_fixed_split_disjoint(self):
  a=[f'animal{i}' for i in range(27)];s=fixed_split(a);self.assertEqual(s,fixed_split(a[::-1]));self.assertEqual([sum(v==k for v in s.values()) for k in ['train','validation','test']],[15,6,6])
 def test_duplicates_not_animals(self):
  with self.assertRaises(ValueError):fixed_split(['a']*27)
 def test_real_windows(self):self.assertEqual(valid_windows(np.ones(20),np.arange(20),10).sum(),10)
 def test_masked_not_filled(self):
  x=np.ones(20);m=np.ones(20,bool);m[10]=False;self.assertEqual(valid_windows(x,np.arange(20),10,m).sum(),0)
 def test_segment_boundaries(self):self.assertEqual(valid_windows(np.ones(20),np.arange(20),10,segment_starts=[10]).sum(),0)
 def test_gap_not_crossed(self):
  t=np.arange(20,dtype=float);t[10:]+=3;self.assertEqual(valid_windows(np.ones(20),t,10).sum(),0)
 def test_bad_time(self):
  with self.assertRaises(ValueError):valid_windows(np.ones(3),[0,1,1],1)
 def test_common_classes_not_marginal(self):
  c=np.array([[1,1,0],[0,1,1]],bool);self.assertFalse(common_core_exists(c,2,2));self.assertTrue(common_core_exists(np.ones((15,15),bool)))
 def test_class_bits_beyond_64(self):
  c=np.zeros((100,15),bool);c[80:95,:]=True;self.assertTrue(common_core_exists(c))
 def test_bilateral_distinction(self):
  r=[dict(trace_row_1based='1',label='AVAL',neuron_class='AVA',LR='L',confidence='5'),dict(trace_row_1based='2',label='AVAR',neuron_class='AVA',LR='R',confidence='5')]
  z=identity_join(r,2);self.assertEqual(len(z),2);self.assertEqual(len({x['neuron_class'] for x in z}),1)
 def test_identity_duplicates(self):
  r=dict(trace_row_1based='1',label='AVAL',confidence='5')
  with self.assertRaises(ValueError):identity_join([r,r],2)
 def test_cleanup_only_two_named_paths(self):self.assertEqual(len(TARGETS),2)
 def test_symlink_preserved(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'source';p.write_text('source');link=Path(d)/'link';link.symlink_to(p)
   with self.assertRaises(ValueError):check_path(link)
   self.assertTrue(p.exists())
 def test_hardlink_preserved(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'source';p.write_text('source');os.link(p,Path(d)/'copy')
   with self.assertRaises(ValueError):check_path(p)
 def test_error_not_recording(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'bad.h5';p.write_text('<html>Error</html>')
   with self.assertRaises(ValueError):verify_format(p.name,p)
if __name__=='__main__':unittest.main()
