"""Fail-closed cohort preflight. Receipt and individual inventory stay private."""
from pathlib import Path
import csv,json,hashlib
ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent

def inspect(recordings,receipt,expected_remaining,expected_all):
    d=json.loads(receipt.read_text(encoding='utf-8-sig'))
    if d.get('algorithm','').upper().replace('-','')!='SHA256':
        raise ValueError('Expected SHA256 receipt')
    rows=d['files'];names=[r['filename'] for r in rows]
    if any(Path(n).name!=n for n in names) or len(names)!=len(set(names)):
        raise ValueError('Unsafe or duplicate receipt names')
    verified=0;missing=0;bad=0
    for row in rows:
        p=recordings/row['filename']
        if not p.is_file() or p.is_symlink():missing+=1;continue
        if p.stat().st_size!=row['bytes'] or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256'].lower():bad+=1
        else:verified+=1
    actual={p.name for p in recordings.glob('*-cudata.h5') if p.is_file() and not p.is_symlink()}
    ready=set(names)==expected_remaining and actual==expected_all and verified==24 and not(missing or bad)
    return {'status':'PASS' if ready else 'BLOCKED_MISSING_COHORT_OR_WRONG_RECEIPT',
            'receipt_entries':len(rows),'receipt_entries_verified':verified,
            'receipt_missing_files':missing,'receipt_hash_or_size_failures':bad,
            'receipt_matches_exact_remaining24_manifest':set(names)==expected_remaining,
            'expected_independent_animals':len(expected_all),'actual_recordings_present':len(actual),
            'expected_recordings_missing':len(expected_all-actual),'unexpected_recordings':len(actual-expected_all),
            'full_cohort_QC_run':False,'previous_QC_overwritten':False}

def main():
    base=ROOT/'data/worm_behavior_sheaf_readiness_v1';r=HERE.parent
    remaining={row['file'] for row in csv.DictReader((r/'export_recovery_v3/REMAINING_FILES.csv').open())}
    expected={row['file'] for row in csv.DictReader((r/'PUBLIC_FILE_MANIFEST.csv').open()) if row['file'].endswith('-cudata.h5')}
    result=inspect(base/'recordings',base/'recordings/remaining_transfer_PRIVATE.json',remaining,expected)
    journals=list(base.rglob('remaining_acquisition_PRIVATE.json'))
    result['failed_acquisition_journals_present']=len(journals)
    result['existing_journals_preserved']=True
    (HERE/'TRANSFER_CHECKPOINT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    return 0 if result['status']=='PASS' else 2
if __name__=='__main__':raise SystemExit(main())
