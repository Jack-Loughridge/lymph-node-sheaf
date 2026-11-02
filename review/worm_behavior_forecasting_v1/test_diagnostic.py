import unittest,tempfile
from pathlib import Path
import numpy as np
from models import moments,combine,fit,predict,mse_from_moments,cross_rank1
import run_diagnostic as run
class DiagnosticTests(unittest.TestCase):
 def setUp(self):self.rng=np.random.default_rng(9)
 def test_exact_cutoff(self):
  t=np.arange(30.);neural=np.array([t,100+t]);beh=np.array([200+t]*5)
  x,y,idx=run.designs(neural,beh,t,3)
  np.testing.assert_array_equal(x[0,:6],[2,1,0,102,101,100]);np.testing.assert_array_equal(y[0],[3,103]);self.assertEqual(idx[0],3)
 def test_no_window_crosses_gap(self):
  t=np.arange(30.);t[15:]+=10;neural=np.array([np.arange(30.)]*2);beh=np.array([np.arange(30.)]*5)
  x,y,idx=run.designs(neural,beh,t,3)
  self.assertFalse(any(k in idx for k in [15,16,17]))
 def test_heldout_split_access_refused(self):
  for s in ['test','validation']:
   self.assertRaises(ValueError,run.require_training,'example',{'example':s})
 def test_equal_animal_normalisation(self):
  x=np.arange(10.).reshape(5,2);y=x.copy();items=[moments(x,y),moments(np.tile(x+100,(3,1)),np.tile(y+100,(3,1)))]
  f=combine(items);np.testing.assert_allclose(f.mx,x.mean(0)+50)
 def test_moment_MSE_matches_predictions(self):
  x=self.rng.normal(size=(200,14));y=x[:,:2]+self.rng.normal(size=(200,2))*.1;m=moments(x,y);f=combine([m])
  for model in ['AR','VAR','AR_beh','VAR_beh','rank1','graph','persistence']:
   b,o,info=fit(f,2,1,model,(False,True));actual=np.mean((predict(x,f,b,o)-y)**2,axis=0);calculated,_=mse_from_moments(m,f,b,o)
   np.testing.assert_allclose(calculated,actual,rtol=1e-9,atol=1e-10)
 def test_graph_omission_exact_AR(self):
  x=self.rng.normal(size=(100,14));y=x[:,:2];f=combine([moments(x,y)])
  b,o,_=fit(f,2,.1,'AR');g,go,_=fit(f,2,.1,'graph',(False,False));np.testing.assert_allclose(b,g);np.testing.assert_allclose(o,go)
 def test_real_rank_one_and_monotonic_objective(self):
  x=self.rng.normal(size=(400,21));k=np.array([.5,.2,-.1]);y=np.c_[x[:,3:6]@k,x[:,:3]@k*2]+self.rng.normal(size=(400,2))*.03
  f=combine([moments(x,y)]);b,info=cross_rank1(f,3,.1)
  cross=np.stack([b[3:6,0],b[:3,1]]);self.assertLess(np.linalg.svd(cross,compute_uv=False)[1],1e-10);self.assertLessEqual(info['maximum_objective_increase'],1e-9)
 def test_joint_day_group_insufficient_is_not_evaluable(self):
  with tempfile.TemporaryDirectory() as t:
   prior=run.PRIVATE;run.PRIVATE=Path(t);(run.PRIVATE/'checkpoints').mkdir();(run.PRIVATE/'predictions').mkdir()
   try:
    arrays={a:(np.ones((30,14)),np.ones((30,2)),np.arange(30)) for a in ['2000-01-01-a','2000-01-01-b','2000-01-02-c']}
    result=run.nested(('N1','N2'),2,'day',arrays,(True,True),'synthetic');self.assertEqual(result['fit_count'],0);self.assertTrue(result['status'].startswith('NOT EVALUABLE'))
   finally:run.PRIVATE=prior
 def test_nested_outer_outcome_does_not_tune_itself(self):
  arrays={}
  for i in range(3):
   x=self.rng.normal(size=(35,14));y=np.c_[x[:,0],x[:,2]]+self.rng.normal(size=(35,2))*.1;arrays[f'2000-01-0{i+1}-animal']=(x,y,np.arange(35))
  with tempfile.TemporaryDirectory() as t:
   prior=run.PRIVATE
   try:
    def execute(sub,arr):
     run.PRIVATE=Path(t)/sub;(run.PRIVATE/'checkpoints').mkdir(parents=True);(run.PRIVATE/'predictions').mkdir();return run.nested(('N1','N2'),2,'animal',arr,(True,False),'synthetic')
    a=execute('first',arrays);key=sorted(arrays)[0];changed=dict(arrays);x,y,idx=arrays[key];changed[key]=(x,y+1000,idx);b=execute('second',changed)
    fa=next(f for f in a['folds'] if f['animal']==key);fb=next(f for f in b['folds'] if f['animal']==key)
    self.assertNotIn(key,fa['fit_animals']);self.assertEqual(fa['inner_scores'],fb['inner_scores'])
    for model in fa['metrics']:self.assertEqual(fa['metrics'][model]['alpha'],fb['metrics'][model]['alpha'])
   finally:run.PRIVATE=prior
if __name__=='__main__':unittest.main()
