"""Exact allowlist cleanup. Never recursively removes directories or protected data."""
from pathlib import Path
import json,hashlib,subprocess,os,shutil,datetime,stat
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
BASE=ROOT/'data';PRIVATE=BASE/'perturbational_bioelectric_integration_2026'
def digest(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def checked_path(p):
 p=Path(p).absolute()
 if p.resolve()!=p or not p.is_relative_to(BASE):raise ValueError('Symlink or outside authorised project data')
 if any(k in p.parts for k in ['humann_databases','AGORA2_SBML_clean','external_validation','out_graphs','.git','julia_tools']):raise ValueError('Protected path')
 if not p.is_file():raise ValueError('Expected existing regular file')
 return p

def handles(p):
 result=subprocess.run(['lsof','-nP','--',str(p)],capture_output=True,text=True,timeout=20)
 warnings=result.stderr.strip()
 if result.returncode not in (0,1):raise ValueError('Open-process inspection uncertain')
 allowed=["lsof: WARNING: can't stat() hugetlbfs file system /dev/hugepages","lsof: WARNING: can't stat() mqueue file system /dev/mqueue",'Output information may be incomplete.']
 if warnings and any(line.strip() not in allowed for line in warnings.splitlines()):raise ValueError('Unexpected lsof warning')
 if result.stdout.strip():raise ValueError('File has active open handle or mapping')
 # Refuse when another same-user process cannot be inspected.
 inaccessible=0
 for proc in Path('/proc').iterdir():
  if not proc.name.isdigit():continue
  try:
   if proc.stat().st_uid!=os.getuid():continue
   fds=list((proc/'fd').iterdir())
   for fd in fds:
    try:
     if fd.resolve()==p:raise ValueError('Active proc file descriptor')
    except FileNotFoundError:pass
   for line in (proc/'maps').read_text().splitlines():
    fields=line.split(maxsplit=5)
    if len(fields)==6 and fields[-1].removesuffix(' (deleted)')==str(p):raise ValueError('Active memory mapping')
  except FileNotFoundError:continue
  except PermissionError:inaccessible+=1
 if inaccessible:raise ValueError('Same-user process FD inspection incomplete')
 return {'lsof_exit':result.returncode,'open_handles':0,'same_user_inaccessible_fd_directories':0,'unrelated_virtual_filesystem_warnings':warnings,'scope':'Visible process namespace; target descriptors and maps explicitly checked; no claim of privileged host-wide inspection'}

def main():
 journal=PRIVATE/'cleanup_manifest_PRIVATE.json'
 if journal.exists():raise ValueError('Preserve existing cleanup journal; use a new version')
 old=BASE/'worm_behavior_sheaf_readiness_v1';qc=json.loads((old/'full_cohort_qc_v5/cohort_results_PRIVATE.json').read_text());native=json.loads((old/'full_cohort_qc_v5/extraction_receipt_PRIVATE.json').read_text())
 source='https://datadryad.org/dataset/doi:10.5061/dryad.w9ghx3g4v';planned=[]
 for animal in qc['animals']:planned.append((old/'recordings'/f"{animal['animal']}-cudata.h5",animal['sha256'],'Recoverable downloaded source data',source))
 import tomllib
 receipt=tomllib.loads((old/'sample_ingestion_v2/native_identity/native_receipt_PRIVATE.toml').read_text())
 planned.append((old/'recordings/NeuroPAL_labels.jld2',receipt['source_sha256'],'Recoverable downloaded source data',source))
 planned.append((old/'incoming/doi_10_5061_dryad_w9ghx3g4v__v20260415.zip',native['archive_sha256'],'Recoverable downloaded archive',source))
 for name in ['selected_raw_counts.int32.memmap','extra_marker_counts.int32.memmap']:
  planned.append((BASE/'crohn_xenium/biological_validation_v1'/name,None,'Regenerable derived cache','https://singlecell.broadinstitute.org/single_cell/study/SCP3259'))
 # Freeze all prior code/report/specification files and private forecasting results.
 preserve={}
 for directory in ['src','tests','review','data/worm_behavior_forecasting_v1']:
  for p in (ROOT/directory).rglob('*'):
   if p.is_file() and not p.is_symlink() and '__pycache__' not in p.parts and not p.is_relative_to(HERE):preserve[str(p.relative_to(ROOT))]=digest(p)
 state={'start_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'free_before_bytes':shutil.disk_usage(ROOT).free,'targets':[],'preserved_files_sha256_PRIVATE':preserve}
 # Preflight every file before deleting any.
 for p,expected,category,url in planned:
  p=checked_path(p);s=p.stat();h=digest(p)
  if expected is not None and h!=expected:raise ValueError('Prior source checksum mismatch')
  state['targets'].append({'path':str(p),'bytes':s.st_size,'allocated_bytes':s.st_blocks*512,'sha256':h,'prior_checksum':expected or 'Not previously recorded; measured before removal','category':category,'source_url':url,'open_process_check':handles(p),'inode':s.st_ino,'device':s.st_dev,'status':'VERIFIED_PENDING'})
 journal.write_text(json.dumps(state,indent=2)+'\n')
 for row in state['targets']:
  p=checked_path(row['path']);s=p.stat();assert (s.st_ino,s.st_dev,s.st_size)==(row['inode'],row['device'],row['bytes']);assert digest(p)==row['sha256'];handles(p)
  p.unlink();row['status']='REMOVED';row['post_deletion_absent']=not p.exists();journal.write_text(json.dumps(state,indent=2)+'\n')
 assert all(digest(ROOT/p)==h for p,h in preserve.items()),'Preserved-file hash changed'
 state['preserved_files_unchanged']=len(preserve);state['free_after_bytes']=shutil.disk_usage(ROOT).free;state['sum_removed_logical_bytes']=sum(r['bytes'] for r in state['targets']);state['free_space_increase_bytes']=state['free_after_bytes']-state['free_before_bytes']
 journal.write_text(json.dumps(state,indent=2)+'\n')
 summary={'removed_files':len(state['targets']),'removed_logical_bytes':state['sum_removed_logical_bytes'],'free_before_bytes':state['free_before_bytes'],'free_after_bytes':state['free_after_bytes'],'measured_free_space_increase_bytes':state['free_space_increase_bytes'],'preserved_code_reports_specs_and_forecasting_outputs_hash_checked':len(preserve),'Dryad_recoverability':'Public30-file manifest verified; original per-file/ZIP integrity receipts retained privately','Crohn_count_caches_recoverability':'Original H5AD public SCP3259 source and extraction code/shape/reference metadata retained','private_exact_manifest':'data/perturbational_bioelectric_integration_2026/cleanup_manifest_PRIVATE.json','NWB_files':'RETAINED current recovery access unverified','Julia_tools':'RETAINED','microbiome_project':'UNTOUCHED','old_IBD_archive':'RETAINED READ_ONLY; guarded normal-SSH command supplied','Crohn_8_98GB_H5AD':'Already absent; not counted as reclaimed now'}
 (HERE/'CLEANUP_AGGREGATE.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
