import ast,io,json,tarfile,tempfile,unittest
from pathlib import Path
import numpy as np
import reconcile as r
import inspect_exports as e

class ReconciliationTests(unittest.TestCase):
 def test_run_at_boundary_is_accepted(self):
  self.assertEqual(r.max_run([0,1,1,0]),2)
 def test_scattered_missing_is_not_total_completeness(self):
  t=np.arange(100.);y=np.ones(100);y[::5]=np.nan
  z=r.source_gate(t,y,30,-30,70,.05)
  self.assertTrue(z['nan_component_pass']);self.assertAlmostEqual(z['finite_fraction'],.8)
 def test_long_run_fails(self):
  y=np.ones(100);y[10:16]=np.nan;self.assertFalse(r.source_gate(np.arange(100),y,30,-30,70,.05)['nan_component_pass'])
 def test_no_interpolation(self):
  y=np.array([1.,np.nan,2.]);copy=y.copy();r.finite_summary(y);np.testing.assert_array_equal(y,copy)
 def test_stable_session_not_export_ordinal(self):
  self.assertEqual(e.session_keys('/project/36/pumpprobe_20220331_114306/'),{'20220331_114306'})
 def test_ambiguous_session_keys_retained(self):
  self.assertEqual(len(e.session_keys('20220331_114306 and 20230314_164209')),2)
 def test_unsafe_archive_path(self):
  m=tarfile.TarInfo('../x_ds_name.txt')
  with self.assertRaises(ValueError):e.safe_member(m)
 def test_unsafe_symlink(self):
  m=tarfile.TarInfo('0_ds_name.txt');m.type=tarfile.SYMTYPE
  with self.assertRaises(ValueError):e.safe_member(m)
 def test_member_budget(self):
  m=tarfile.TarInfo('0_gcamp.txt');m.size=e.MAX_MEMBER+1
  with self.assertRaises(ValueError):e.safe_member(m)
 def test_source_nan_flag_discrepancy(self):
  path=r.OUT/'responses_detect.py'
  if not path.exists():self.skipTest('Pinned external source cache absent; fetch_sources.py first')
  tree=ast.parse(path.read_text());assign=[n for n in ast.walk(tree) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='nan_thresh' for t in n.targets)]
  self.assertEqual(len(assign),1);self.assertEqual(ast.literal_eval(assign[0].value),.3)
  self.assertNotIn('--nan-th',{n.value for n in ast.walk(tree) if isinstance(n,ast.Constant) and isinstance(n.value,str)})
 def test_plain_text_export_alignment(self):
  with tempfile.TemporaryDirectory() as td:
   out=Path(td);path=out/'fixture.tar.gz';blobs={
    '7_ds_name.txt':'/private/path/pumpprobe_20220331_114306/',
    '7_labels.txt':'AVEL\nAVER\n','7_t.txt':'0\n0.5\n1\n',
    '7_gcamp.txt':'1 2\n2 3\n3 4\n','7_stim_volume_i.txt':'1\n','7_stim_neurons.txt':'0\n'}
   with tarfile.open(path,'w:gz') as tar:
    for name,content in blobs.items():
     data=content.encode();m=tarfile.TarInfo(name);m.size=len(data);tar.addfile(m,io.BytesIO(data))
   rows=[dict(subject='sub-36',event_index=0,nwb_total_events=1,nwb_total_frames=3,nearest_nwb_stimulus_frame=1,target_name='AVEL')]
   status=e.reconcile_text_archive(path,rows,[dict(subject='sub-36',filename='sub-36_ses-20220331-114306_desc-segmentation_ophys+ogen.nwb')],out)
   self.assertEqual(rows[0]['original_export_session_match'],'MATCHED_STABLE_SESSION');self.assertTrue(rows[0]['original_target_label_agrees'])
   self.assertIn('NOT_PRESENT',rows[0]['original_autoresponse_flag']);self.assertIn('NOT_EVALUABLE',rows[0]['original_published_inclusion'])
 def test_version_mismatch_not_silent(self):
  with tempfile.TemporaryDirectory() as td:
   path=Path(td)/'exported_data.tar.gz';path.write_bytes(b'not the reported version')
   with self.assertRaises(ValueError):e.reconcile_text_archive(path,[],[],Path(td))
if __name__=='__main__':unittest.main()
