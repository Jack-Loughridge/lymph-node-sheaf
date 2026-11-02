import unittest,tempfile,json
from pathlib import Path
from unittest.mock import patch
import h5py,numpy as np,pandas as pd
import spatial_audit as a

class SpatialAuditTests(unittest.TestCase):
 def test_complete_complex_is_AND(self):
  x=np.array([[1,0],[0,1],[1,1]]);np.testing.assert_array_equal(a.co_detect(x,[0,1],1),[False,False,True])
 def test_threshold_requires_all_subunits(self):
  x=np.array([[2,1],[2,2]]);np.testing.assert_array_equal(a.co_detect(x,[0,1],2),[False,True])
 def test_missing_category_not_last_label(self):
  np.testing.assert_array_equal(a.label_mask([True,False],np.array([-1,0,1])),[False,True,False])
 def test_bad_categorical_code(self):
  with self.assertRaises(ValueError):a.label_mask([True,False],np.array([2]))
 def test_panel_incompatible_markers_are_unevaluable(self):
  result=a.marker_availability({'PDGFRA'},a.FIBRO_MARKERS);self.assertFalse(result['evaluable']);self.assertEqual(result['available'],['PDGFRA'])
 def test_marker_availability_does_not_impute(self):
  self.assertFalse(a.marker_availability(set(),['FAP','PDPN'])['evaluable'])
 def test_no_sender_no_exposure(self):
  self.assertFalse(a.spatial_exposed(np.array([[0.,0.]]),np.empty((0,2)),5).any())
 def test_physical_radius_boundary(self):
  np.testing.assert_array_equal(a.spatial_exposed(np.array([[0.,0.],[1.,0.]]),np.array([[5.,0.]]),4),[False,True])
 def test_sparse_rows_and_offsets(self):
  with tempfile.TemporaryDirectory() as td:
   with h5py.File(Path(td)/'test.h5','w') as f:
    g=f.create_group('counts');g['indptr']=[0,2,3,4];g['indices']=[2,0,1,0];g['data']=[5.,3.,2.,7.]
    arr=a.read_csr_rows(g,1,3,3).toarray();np.testing.assert_array_equal(arr,[[0,2,0],[7,0,0]])
 def test_shared_slide_resolved_by_explicit_biobank(self):
  wb=pd.DataFrame({'Xenium_slide':['0000042','0000042'],'Biobank ID':['BIOA','BIOB'],'Table S1 patient ID':['1','2']})
  with tempfile.TemporaryDirectory() as td,patch.object(a,'LOCAL',Path(td)):
   rows=a.map_samples(['0000042__bioa_inf','0000042__bioa_nor','0000042__biob_inf','0000042__biob_nor'],wb)
   self.assertEqual([r['patient_id'] for r in rows],['1','1','2','2']);self.assertTrue(all(not r['S1_verified'] for r in rows))
 def test_ambiguous_sample_not_guessed(self):
  wb=pd.DataFrame({'Xenium_slide':['0000042','0000042'],'Biobank ID':['BIOA','BIOB'],'Table S1 patient ID':['1','2']})
  with tempfile.TemporaryDirectory() as td,patch.object(a,'LOCAL',Path(td)):
   with self.assertRaises(ValueError):a.map_samples(['0000042_affected'],wb)
 def test_unknown_condition_not_called_healthy(self):
  wb=pd.DataFrame({'Xenium_slide':['0000042'],'Biobank ID':['BIOA'],'Table S1 patient ID':['1']})
  with tempfile.TemporaryDirectory() as td,patch.object(a,'LOCAL',Path(td)):
   with self.assertRaises(ValueError):a.map_samples(['0000042_unlabelled'],wb)
if __name__=='__main__':unittest.main()
