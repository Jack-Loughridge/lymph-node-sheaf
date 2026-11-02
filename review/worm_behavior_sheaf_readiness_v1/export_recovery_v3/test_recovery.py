import unittest
import numpy as np,pandas as pd
from reconcile_native_identities import primary_identities,common_core_witness
from acquire_remaining_processed import discover_links
from audit_behaviour_events import band_events
class RecoveryTests(unittest.TestCase):
 def test_note_and_invalid_holds(self):
  d=pd.DataFrame({'identity_valid':['true','true','true','false','true'],'confidence':[4,5,6,5,3],'region':['',' merge warning','','','']})
  self.assertEqual(primary_identities(d).index.tolist(),[0])
 def test_joint_not_marginal_coverage(self):
  c=np.array([[1,1,0],[1,0,1],[0,1,1]],bool)
  self.assertIsNone(common_core_witness(c,2,2))
 def test_more_than64_classes(self):
  c=np.ones((70,15),bool);w=common_core_witness(c,15,15)
  self.assertEqual(len(w[1]),70)
 def test_real_host_allowlist(self):
  name='example-cudata.h5'
  self.assertRaises(ValueError,discover_links,f'<a href="https://evil.invalid/downloads/file_stream/1">{name}</a>',[name])
  self.assertEqual(discover_links(f'<a href="/downloads/file_stream/1">{name}</a>',[name])[name],'https://datadryad.org/downloads/file_stream/1')
 def test_unknown_links_rejected(self):
  self.assertRaises(ValueError,discover_links,'', ['missing-cudata.h5'])
 def test_event_history_and_right_censor(self):
  t=np.arange(21.);distance=np.zeros(21);distance[8:10]=10000;distance[19:]=10000
  events=band_events(t,distance,10000,.2,np.ones(21),np.zeros(21))
  self.assertEqual([e['frame_zero_based'] for e in events],[8])
if __name__=='__main__':unittest.main()
