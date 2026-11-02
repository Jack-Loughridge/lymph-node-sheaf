"""Allowlisted ZIP extraction, no overwrite. Private manifest and checkpoints."""
from pathlib import Path,PurePosixPath
import zipfile,stat,hashlib,json,csv,shutil,os
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
BASE=ROOT/'data/worm_behavior_sheaf_readiness_v1'
def digest(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def select_members(infos,names):
 selected={n:[] for n in names}
 for i in infos:
  path=PurePosixPath(i.filename)
  if '\\' in i.filename or path.is_absolute() or '..' in path.parts or ':' in i.filename or '\x00' in i.filename:
   raise ValueError('Unsafe archive member path')
  mode=i.external_attr>>16
  if stat.S_ISLNK(mode) or (stat.S_IFMT(mode) not in (0,stat.S_IFREG,stat.S_IFDIR)):
   raise ValueError('Nonregular archive member')
  if path.name in selected:
   if i.is_dir() or i.flag_bits&1 or not 0<i.file_size<=20_000_000:raise ValueError('Unexpected selected member type or size')
   selected[path.name].append(i)
 if any(len(v)!=1 for v in selected.values()):raise ValueError('Each expected processed recording must occur exactly once')
 return {n:v[0] for n,v in selected.items()}
def main():
 archive=BASE/'incoming/doi_10_5061_dryad_w9ghx3g4v__v20260415.zip';dest=BASE/'recordings';private=BASE/'full_cohort_qc_v5'
 if private.exists():raise ValueError('Private checkpoint exists; preserve before a versioned rerun')
 if archive.is_symlink() or dest.is_symlink():raise ValueError('Symlink destination/source refused')
 rows=list(csv.DictReader((HERE.parent/'export_recovery_v3/REMAINING_FILES.csv').open()));names={r['file'] for r in rows}
 if len(names)!=24:raise ValueError('Expected exact missing24 allowlist')
 with zipfile.ZipFile(archive) as z:
  selected=select_members(z.infolist(),names)
  if any((dest/n).exists() or (dest/n).is_symlink() for n in names):raise ValueError('Refusing to overwrite any selected recording')
  if shutil.disk_usage(dest).free<sum(i.file_size for i in selected.values())+1_000_000_000:raise ValueError('Insufficient reserve')
  private.mkdir();historical={str(p.relative_to(BASE)):digest(p) for p in BASE.rglob('*') if p.is_file() and (p.name in ['remaining_acquisition_PRIVATE.json','remaining_transfer_PRIVATE.json','transfer_receipt_PRIVATE.json'])}
  state={'archive_sha256':digest(archive),'archive_bytes':archive.stat().st_size,'archive_manifest':[{'name':i.filename,'bytes':i.file_size,'CRC32':f'{i.CRC:08x}'} for i in z.infolist()],'historical_files_before':historical,'files':[],'windows_archive_integrity':'USER_REPORTED_MATCH; sender digest not independently available'}
  journal=private/'extraction_receipt_PRIVATE.json'
  journal.write_text(json.dumps(state,indent=2)+'\n')
  for n in sorted(names):
   i=selected[n];out=dest/n;temp=dest/(n+'.extracting')
   # ZipFile.open validates CRC at EOF. Exclusive creation prevents clobbering.
   with z.open(i) as src,temp.open('xb') as f:shutil.copyfileobj(src,f,1024*1024)
   with temp.open('rb') as f:magic=f.read(8)
   if magic!=b'\x89HDF\r\n\x1a\n' or temp.stat().st_size!=i.file_size:raise ValueError('HDF5 magic or length failure')
   # Atomic no-overwrite publication; do not replace or follow existing paths.
   os.link(temp,out);temp.unlink()
   state['files'].append({'filename':n,'bytes':out.stat().st_size,'sha256':digest(out),'zip_CRC32':f'{i.CRC:08x}','zip_CRC_verified':True})
   journal.write_text(json.dumps(state,indent=2)+'\n')
  state['historical_files_unchanged']=all(digest(BASE/p)==h for p,h in historical.items())
  state['algorithm']='SHA256';journal.write_text(json.dumps(state,indent=2)+'\n')
  expected={r['file'] for r in csv.DictReader((HERE.parent/'PUBLIC_FILE_MANIFEST.csv').open()) if r['file'].endswith('-cudata.h5')}
  actual={p.name for p in dest.glob('*-cudata.h5')}
  if actual!=expected:raise ValueError('Final cohort membership mismatch')
  for row in state['files']:
   if digest(dest/row['filename'])!=row['sha256']:raise ValueError('Post-extraction hash mismatch')
  aggregate={'archive_bytes':state['archive_bytes'],'archive_sha256':state['archive_sha256'],'manifest_entries':len(z.infolist()),'expected_selected_files':24,'uniquely_matched_files':24,'extracted_files':len(state['files']),'extracted_bytes':sum(r['bytes'] for r in state['files']),'CRC_and_post_extraction_SHA256_checks':'PASS','expected_animal_recordings':27,'actual_animal_recordings':len(actual),'historical_receipts_and_journals_unchanged':state['historical_files_unchanged'],'windows_archive_hash_match':'USER_REPORTED_NOT_INDEPENDENTLY_RECHECKED'}
  (HERE/'EXTRACTION_AGGREGATE.json').write_text(json.dumps(aggregate,indent=2)+'\n');print(json.dumps(aggregate,indent=2))
if __name__=='__main__':main()
