"""Measurement/metadata checks only; no transformations fitted or signals filled."""
import hashlib
from itertools import combinations
import numpy as np

def fixed_split(animal_ids):
 ids=list(animal_ids)
 if len(ids)!=27 or len(set(ids))!=27:raise ValueError('Freeze exactly the 27 distinct manifest IDs before exclusions')
 key=lambda x:hashlib.sha256(('worm-readiness-v1-20261102|'+x).encode()).hexdigest()
 order=sorted(ids,key=key)
 return {a:('train' if i<15 else 'validation' if i<21 else 'test') for i,a in enumerate(order)}

def valid_windows(trace,timestamps,history=10,observed_mask=None,segment_starts=()):
 trace=np.asarray(trace,float);t=np.asarray(timestamps,float)
 if trace.ndim!=1 or t.ndim!=1 or len(trace)!=len(t):raise ValueError('Axes mismatch')
 dt=np.diff(t)
 if not np.isfinite(t).all() or np.any(dt<=0):raise ValueError('Bad timestamps; no silent repair')
 if len(t)<=history:return np.zeros(0,bool)
 finite=np.isfinite(trace)
 if observed_mask is not None:
  obs=np.asarray(observed_mask)
  if obs.dtype!=bool or obs.shape!=trace.shape:raise ValueError('Explicit boolean observation mask required')
  finite &= obs
 breaks=np.zeros(len(t),bool);breaks[1:]=dt>1.5*np.median(dt)
 for k in segment_starts:
  if int(k)!=k or not 0<=k<len(t):raise ValueError('Invalid segment start')
  breaks[int(k)]=True
 good=np.convolve(finite.astype(int),np.ones(history+1,dtype=int),'valid')==history+1
 # A target at k uses history [k-history,k); breaks at its first point are allowed,
 # but any later break would join distinct sessions or cross a missing interval.
 for j in range(len(good)):
  if breaks[j+1:j+history+1].any():good[j]=False
 return good

def common_core_exists(coverage,min_classes=15,min_animals=15):
 """Exact structural gate on shared presence; no outcome-based selection.
 DFS prunes by shared-class intersections. Returns witness metadata only, not
 a predictive feature selector. This can be costly for pathological matrices.
 """
 c=np.asarray(coverage,bool)
 if c.ndim!=2:raise ValueError('classes x animals required')
 if c.shape[0]<min_classes or c.shape[1]<min_animals:return False
 masks=[sum((1<<int(i)) for i in np.flatnonzero(c[:,j])) for j in range(c.shape[1])]
 masks=sorted(masks,key=lambda x:-x.bit_count())
 def visit(start,chosen,mask):
  if mask.bit_count()<min_classes:return False
  if chosen==min_animals:return True
  if len(masks)-start<min_animals-chosen:return False
  for j in range(start,len(masks)-(min_animals-chosen)+1):
   if visit(j+1,chosen+1,mask&masks[j]):return True
  return False
 return visit(0,0,(1<<c.shape[0])-1)

def identity_join(rows,n_neurons):
 seen_rows=set();seen_individual=set();out=[]
 for r in rows:
  k=int(r['trace_row_1based'])-1
  if k<0 or k>=n_neurons or k in seen_rows:raise ValueError('Invalid/duplicate trace row; never substitute ROI index')
  seen_rows.add(k);conf=float(r['confidence'])
  if not np.isfinite(conf) or not 1<=conf<=5:raise ValueError('Bad identity confidence')
  label=r.get('label','');cls=r.get('neuron_class','')
  if conf>=4 and label:
   if label in seen_individual:raise ValueError('Duplicate individual identity')
   seen_individual.add(label)
  out.append(dict(row=k,label=label,neuron_class=cls,confidence=conf,LR=r.get('LR','')))
 return out
