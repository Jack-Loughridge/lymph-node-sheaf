"""Fail-closed preflight; no training or response-score construction."""
from pathlib import Path
from collections import defaultdict
import json,hashlib
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent

def paired_folds(rows):
 groups=defaultdict(list)
 seen=set()
 for r in rows:
  if not r.get('S1_verified'):raise ValueError('S1 verification required')
  key=r['specimen_id']
  if key in seen:raise ValueError('Duplicate specimen')
  seen.add(key);groups[r['patient_id']].append(r)
 if len(groups)!=8:raise ValueError('Eight independently verified patients required')
 for vals in groups.values():
  if len(vals)!=2 or {v['status'] for v in vals}!={'affected','adjacent'}:raise ValueError('Each complete patient pair required')
 return [(sorted(s['specimen_id'] for p,v in groups.items() if p!=hold for s in v),sorted(s['specimen_id'] for s in groups[hold])) for hold in sorted(groups)]

def check_features(inputs,outcomes,annotation_genes,ligands_receptors):
 bad=set(inputs)&(set(outcomes)|set(annotation_genes)|set(ligands_receptors))
 if bad:raise ValueError('Prohibited input reuse: '+','.join(sorted(bad)))
 return True

def verify_frozen():
 for manifest in ['FROZEN_ORIGINAL_SHA256.json','PREREGISTRATION_SHA256.json']:
  hashes=json.loads((HERE/manifest).read_text())
  for name,wanted in hashes.items():
   p=ROOT/name
   if hashlib.sha256(p.read_bytes()).hexdigest()!=wanted:raise ValueError('Frozen file changed: '+name)
 return True

def main():
 verify_frozen()
 needed=['verified_S1_crosswalk_PRIVATE.json','verified_response_dictionary_PRIVATE.json','verified_competence_ledger_PRIVATE.json','verified_anatomy_and_controls_PRIVATE.json']
 private=ROOT/'data/crohn_xenium/spatial_context_v1'
 missing=[name for name in needed if not (private/name).exists()]
 # File presence alone never marks scientifically verified: explicit manual review still required.
 result={'status':'BLOCKED_NOT_EVALUABLE' if missing else 'REVIEW_REQUIRED_NOT_AUTHORIZED_FOR_FIT','missing_reviewed_artifacts':missing,'prediction_fitted':False,'original_NO_GO_unchanged':True,'preregistration_hashes_verified':True}
 (HERE/'READINESS.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
