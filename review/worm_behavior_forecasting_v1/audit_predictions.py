"""Independent numerical audit of saved private outer predictions and split barriers."""
from pathlib import Path
import json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];PRIVATE=ROOT/'data/worm_behavior_forecasting_v1'
def main():
 results=json.loads((PRIVATE/'all_results_PRIVATE.json').read_text());split=json.loads((ROOT/'data/worm_behavior_sheaf_readiness_v1/animal_splits_PRIVATE.json').read_text());counts=0;targetwindows=0;predictionfiles=0
 prefit=json.loads((HERE/'PREFIT_CODE_RECEIPT.json').read_text());assert all(hashlib.sha256((HERE/n).read_bytes()).hexdigest()==h for n,h in prefit['SHA256'].items()),'Pre-fit code/config changed'
 for panel in results:
  for f in panel['folds']:
   assert split[f['animal']]=='train' and all(split[a]=='train' for a in f['fit_animals'])
   assert f['animal'] not in f['fit_animals']
   if panel['grouping']=='day':assert all(a[:10]!=f['acquisition_day'] for a in f['fit_animals'])
   file=PRIVATE/'predictions'/f"{panel['pair'][0]}_{panel['pair'][1]}_H{panel['history']}_{panel['grouping']}_{f['animal']}_PRIVATE.npz"
   assert hashlib.sha256(file.read_bytes()).hexdigest()==f['prediction_sha256']
   with np.load(file) as z:
    assert len(z['targets'])==f['n_windows'];assert np.all(np.diff(z['target_indices'])>0);assert z['target_indices'].min()>=panel['history']
    for name,metrics in f['metrics'].items():
     direct=np.mean((z[name]-z['targets'])**2,axis=0)
     np.testing.assert_allclose(direct,metrics['MSE_raw'],rtol=1e-12,atol=1e-12);counts+=2
    targetwindows+=len(z['targets']);predictionfiles+=1
 frozen=json.loads((ROOT/'review/worm_behavior_sheaf_readiness_v1/FROZEN_PROTOCOL_SHA256.json').read_text());assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in frozen.items())
 receipt=json.loads((PRIVATE/'run_receipt_PRIVATE.json').read_text());assert all(split[a]=='train' for a in receipt['recording_sources_PRIVATE'])
 report={'prediction_files_hash_checked':predictionfiles,'directed_model_MSE_checks_passed':counts,'target_windows_across_all_descriptive_sensitivities':targetwindows,'not_independent_replicates':True,'outer_animal_and_day_exclusions_checked':True,'recording_sources_all_frozen_training':True,'prefit_code_hashes_unchanged':True,'six_original_protocol_hashes_unchanged':True,'validation_test_recordings_opened':False}
 (HERE/'PREDICTION_INTEGRITY.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
