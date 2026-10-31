import csv
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd
from scipy.linalg import null_space

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src/signalling_sheaf'))
from algebra import operator,solve
from allocation import allocate,cap_flux,INPUT_GENES
from prepare import strict_join,labels,ANNOTATION_GENES,signatures
from run_pilot import permute_bundles,evaluate,score_responses


class AlgebraTests(unittest.TestCase):
 def test_multiple_senders_one_receiver_exact_section(self):
  # One receiver sees SUM 2+3, not two separate requirements r=2 and r=3.
  op=operator(np.array([[2.],[3.]]),[0,1],[0,0],[0,0],[2.,3.],[5.],0,1,1)
  np.testing.assert_allclose(op.D@op.x0-op.h,0,atol=1e-12)
  self.assertLess(solve(op)['energy_relaxed'],1e-12)
  self.assertTrue(set(op.free).isdisjoint({0,1,len(op.x0)-1}))

 def test_eliminated_energy_equals_true_incidence_sheaf(self):
  rng=np.random.default_rng(49)
  for degree in [1,2,5,9]:
   A=rng.normal(size=(1,degree));b=.7;x=rng.normal(size=degree)
   Z=null_space(A);z0=(A.T*b/(A@A.T)[0,0]).ravel()
   u=np.linalg.lstsq(Z,x-z0,rcond=None)[0]
   lifted=float(np.sum((x-z0-Z@u)**2))
   eliminated=float((A@x-b)[0]**2/np.sum(A*A))
   self.assertAlmostEqual(lifted,eliminated,places=11)

 def test_operator_psd_stationarity_and_observation_anchoring(self):
  op=operator(np.array([[2.],[3.]]),[0,1],[0,0],[0,0],[1.,2.],[-2.],1,2,.5)
  stats=solve(op)
  self.assertGreater(stats['energy_relaxed'],0)
  self.assertLessEqual(stats['energy_relaxed'],stats['energy_fixed']+1e-10)
  self.assertLess(stats['stationarity'],1e-8)
  A=op.D[:,op.free].toarray();eig=np.linalg.eigvalsh(A.T@A)
  self.assertGreaterEqual(eig.min(),-1e-10)

 def test_disconnected_receiver_kept_and_intrinsic_baseline(self):
  op=operator(np.array([[1.]]),[],[],[],[],[2.],0,1,0)
  stats=solve(op)
  self.assertAlmostEqual(stats['energy_fixed'],4)
  self.assertAlmostEqual(stats['energy_relaxed'],4)
  self.assertEqual(op.sizes['receivers'],1)

 def test_bmp_routes_sum_into_one_coordinate(self):
  op=operator(np.array([[3.,4.]]),[0,0],[0,0],[0,1],[3.,4.],[7.],0,1,1)
  np.testing.assert_allclose(op.D@op.x0-op.h,0,atol=1e-12)
  self.assertEqual(len(op.blocks['response']),1)


class DataAndAllocationTests(unittest.TestCase):
 def test_composite_join_cannot_cross_slides(self):
  expr=pd.DataFrame({'prefix':['A','B'],'fov':[1,1],'cell':['same','same'],'cell_ID':[1,1]})
  meta=pd.DataFrame({'prefix':['B','A'],'fov':[1,1],'cell':['same','same'],'cell_ID':[1,1],'x':[200,100]})
  self.assertEqual(strict_join(expr,meta).x.tolist(),[100,200])
  with self.assertRaises(ValueError):strict_join(expr,pd.concat([meta,meta.iloc[:1]]))
  with self.assertRaises(ValueError):strict_join(expr,meta.iloc[:1])
  bad=meta.copy();bad.loc[0,'cell_ID']=3
  with self.assertRaises(ValueError):strict_join(expr,bad)

 def test_proxy_annotations_preserve_ambiguous_cells(self):
  genes=ANNOTATION_GENES;ix={g:i for i,g in enumerate(genes)}
  X=np.zeros((5,len(genes)))
  for k in range(4):
   X[k,ix['EPCAM']]=1;X[k,ix['KRT8']]=1
  for g in ['MKI67','TOP2A']:X[[0,2],ix[g]]=1
  for g in ['GUCA2B','FABP1']:X[[1,2],ix[g]]=1
  for g in ['MUC2','TFF3']:X[3,ix[g]]=1
  result,_=labels(X,genes)
  self.assertEqual(result.tolist(),['TA_like','colonocyte_like','epithelial_unresolved','epithelial_unresolved','other_or_mixed'])

 def test_receptor_subunits_spatial_mask_no_self_and_budget(self):
  genes=['IFNG','IFNGR1','IFNGR2']
  X=np.array([[2,0,0],[0,1,0],[0,1,1],[2,1,1]],float)
  coords=np.array([[0,0],[.1,0],[.2,0],[100,0]],float)
  a=allocate(X,genes,coords,[1,2,3],'IFNG',1)
  self.assertTrue((a['rec']==1).all())
  self.assertTrue((a['send']==0).all())
  self.assertLessEqual(a['flux'].sum(),1)
  self.assertEqual(a['exposure'][0],0)
  self.assertEqual(a['exposure'][2],0)

 def test_bmp_shared_receptor_not_double_counted(self):
  genes=['BMP2','BMP6','BMPR1A','ACVR1','BMPR2']
  X=np.array([[5,5,0,0,0],[0,0,1,1,1]],float)
  a=allocate(X,genes,np.array([[0,0],[1,0]],float),[1],'BMP',2)
  self.assertAlmostEqual(a['flux'].sum(),1)
  self.assertEqual(set(a['route']),{0,1})
  self.assertAlmostEqual(a['exposure'][0],1)

 def test_bundle_shuffle_preserves_strata_and_gene_bundles(self):
  X=np.arange(32).reshape(8,4)
  changed,moved=permute_bundles(X,np.array(['a']*4+['b']*4),np.random.default_rng(7))
  self.assertGreater(moved,0)
  for lo,hi in [(0,4),(4,8)]:
   self.assertEqual(sorted(map(tuple,X[lo:hi])),sorted(map(tuple,changed[lo:hi])))

 def test_signed_manifest_keeps_input_exclusion_and_bmp_pooling(self):
  fields=['factor','sheet','gene','padj','log2fc','in_panel','strong','direction','source_usable']
  rows=[]
  for sheet in ['A_TA_DEGs','D_Colonocyte_DEGs']:
   for factor in ['IFNG','TNF','BMP2','BMP6']:
    for gene in ['GENE1','GENE2','TNF','MKI67']:
     fc=-2 if gene=='GENE2' else 2
     rows.append(dict(zip(fields,[factor,sheet,gene,.01,fc,1,1,'down' if fc<0 else 'up',1])))
  with tempfile.TemporaryDirectory() as d:
   path=Path(d)/'sig.csv'
   with path.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
   manifest,weights=signatures(path)
  for w in weights.values():
   self.assertEqual(set(w),{'GENE1','GENE2'})
   self.assertEqual(w['GENE1'],.5);self.assertEqual(w['GENE2'],-.5)
  self.assertTrue(all(r['exclusion_reason']=='ligand_or_receptor_input' for r in manifest if r['gene']=='TNF'))

 def test_response_does_not_enter_allocation_features(self):
  # Genes outside the input/receptor set cannot change allocation.
  X=np.array([[2,0,0,100],[0,1,1,1]],float)
  genes=['IFNG','IFNGR1','IFNGR2','RESPONSE']
  a=allocate(X,genes,np.array([[0,0],[1,0]]),[1],'IFNG',2)
  X[:,-1]=[0,99999]
  b=allocate(X,genes,np.array([[0,0],[1,0]]),[1],'IFNG',2)
  np.testing.assert_array_equal(a['flux'],b['flux'])

 def test_fixed_residual_zero_gain_is_intrinsic(self):
  r=np.array([-2,1,3.])
  mse,e=evaluate(r,np.array([100,0,500]),1,3,0)
  np.testing.assert_array_equal(e,r)
  self.assertAlmostEqual(mse,np.mean(r*r))

 def test_commot_not_silently_substituted(self):
  if importlib.util.find_spec('commot') is None:
   with self.assertRaisesRegex(RuntimeError,'COMMOT backend requires'):
    allocate(np.ones((2,3)),['IFNG','IFNGR1','IFNGR2'],np.array([[0,0],[1,0]]),[1],'IFNG',2,'commot')

if __name__=='__main__':unittest.main()
