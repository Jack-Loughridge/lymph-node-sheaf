import unittest,numpy as np
from sample_numerical_audit import canonical_traces,canonical_reversal_table,exported_reversal_events
class SampleLayoutTests(unittest.TestCase):
 def test_julia_hdf5_axes(self):
  m=np.arange(24).reshape(8,3);x=canonical_traces(m,8);self.assertEqual(x.shape,(3,8));np.testing.assert_array_equal(x[1],m[:,1])
 def test_reject_axis_guess(self):
  with self.assertRaises(ValueError):canonical_traces(np.ones((3,8)),8)
 def test_event_table_axes(self):
  a=np.array([[1,10,20],[5,14,23]]);self.assertEqual(canonical_reversal_table(a).shape,(3,2))
 def test_source_reversal_convention(self):
  np.testing.assert_array_equal(exported_reversal_events([0,0,1,1,0,0]),[[2,4]])
 def test_initial_left_censor(self):
  np.testing.assert_array_equal(exported_reversal_events([1,1,0,0,1,1,0]),[[1,2],[4,6]])
 def test_no_nonbinary(self):
  with self.assertRaises(ValueError):exported_reversal_events([0,2,0])
 def test_no_uncensored_terminal_guess(self):
  with self.assertRaises(ValueError):exported_reversal_events([0,1,1])
if __name__=='__main__':unittest.main()
