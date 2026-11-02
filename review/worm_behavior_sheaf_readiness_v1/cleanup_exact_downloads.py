"""Exact two-file cleanup authorized by user; no recursive deletion.
Inventory/checksum/proc scan first; record every outcome durably. Paths cannot
be supplied via CLI. Run normally for inventory or --execute for exact deletion.
"""
from pathlib import Path
import os,stat,hashlib,json,subprocess,argparse,shutil,datetime
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
PRIVATE=ROOT/'data/worm_behavior_sheaf_readiness_v1';PRIVATE.mkdir(exist_ok=True)
TARGETS=[dict(path='/home/Jack/Real_Data/lymph_node_sheaf/data/crohn_xenium/scp3259/xenium5k_scVI_CellCharter_batch2_reannot_SCP.h5ad',expected_bytes=8984620838,existing_sha256='74759d2d474964dcb537533127f0683bdb9180030f3b94a8f3f8271a061494bb',source_url='https://singlecell.broadinstitute.org/single_cell/study/SCP3259',recoverability='Same named public processed object; source hash and prior source inventory preserved',completed_analysis='review/crohn_xenium_biological_validation_v1/REPORT.md'),dict(path='/home/Jack/Real_Data/spatial_ibd_sheaf/archives/GSE312415_files.tar.gz',expected_bytes=6464275219,existing_sha256=None,source_url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE312415',recoverability='GEO supplementary payloads; archive member list retained in spatial_ibd_sheaf/metadata/GSE312415_archive_members.txt; regenerated bundle bytes may differ',completed_analysis='review/signalling_sheaf/BIOLOGICAL_FINAL_DECISION.txt')]

def check_path(p):
 if not p.is_absolute():raise ValueError('Absolute path required')
 for part in [p,*p.parents]:
  if part.is_symlink():raise ValueError('Symlink component forbidden')
 s=p.lstat()
 if not stat.S_ISREG(s.st_mode) or s.st_nlink!=1:raise ValueError('Only singly-linked regular downloads permitted')
 return s

def open_handles(p):
 found=[];unreadable=0
 for entry in Path('/proc').iterdir():
  if not entry.name.isdigit():continue
  try:
   for fd in (entry/'fd').iterdir():
    try:
     if os.readlink(fd).removesuffix(' (deleted)')==str(p):found.append(int(entry.name))
    except (FileNotFoundError,PermissionError):pass
  except PermissionError:unreadable+=1
  except FileNotFoundError:pass
 r=subprocess.run(['lsof','-t','--',str(p)],capture_output=True,text=True)
 found+= [int(s) for s in r.stdout.split() if s.isdigit()]
 return dict(pids=sorted(set(found)),lsof_exit=r.returncode,lsof_warnings=r.stderr.strip(),inaccessible_proc_fd_directories=unreadable,scope='Accessible proc descriptors and lsof; no claim of privileged inspection')

def write(state):
 p=HERE/'cleanup_manifest.json';temp=p.with_suffix('.tmp');temp.write_text(json.dumps(state,indent=2)+'\n');temp.replace(p)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--execute',action='store_true');args=ap.parse_args()
 if (HERE/'cleanup_manifest.json').exists():raise SystemExit('Existing cleanup checkpoint: do not overwrite; inspect it before any retry')
 state=dict(start_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),free_before_bytes=shutil.disk_usage(ROOT).free,targets=[],policy='Exact user-named obsolete downloads only; preserve derived caches/results, all metadata and worm recordings')
 for t in TARGETS:
  row=dict(t,status='PENDING');state['targets'].append(row);p=Path(t['path'])
  try:
   s=check_path(p);row.update(bytes=s.st_size,allocated_bytes=s.st_blocks*512,inode=s.st_ino,device=s.st_dev)
   if s.st_size!=t['expected_bytes']:raise ValueError('Unexpected source size')
   if not (ROOT/t['completed_analysis']).is_file():raise ValueError('Completed analysis evidence absent')
   row['open_handles_before_hash']=open_handles(p)
   if row['open_handles_before_hash']['pids']:raise ValueError('Active file users detected')
   h=hashlib.sha256()
   with p.open('rb') as f:
    for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
   row['sha256']=h.hexdigest()
   if t['existing_sha256'] and row['sha256']!=t['existing_sha256']:raise ValueError('Frozen checksum mismatch')
   s2=check_path(p)
   if (s.st_ino,s.st_size,s.st_mtime_ns)!=(s2.st_ino,s2.st_size,s2.st_mtime_ns):raise ValueError('File changed during verification')
   row['open_handles_before_delete']=open_handles(p)
   if row['open_handles_before_delete']['pids']:raise ValueError('Active file users detected')
   row['status']='VERIFIED_ELIGIBLE';write(state)
   if args.execute:
    # No globbing, directory removal, symlink traversal or permission escalation.
    try:p.unlink();row['status']='REMOVED'
    except PermissionError as e:row.update(status='FILESYSTEM_PERMISSION_BLOCKED_PRESERVED',error=str(e))
  except (OSError,ValueError) as e:row.update(status='PRESERVED_CHECK_FAILED',error=str(e))
  state['free_after_bytes']=shutil.disk_usage(ROOT).free;write(state)
 state.update(finished_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),removed_logical_bytes=sum(r.get('bytes',0) for r in state['targets'] if r['status']=='REMOVED'),removed_allocated_bytes=sum(r.get('allocated_bytes',0) for r in state['targets'] if r['status']=='REMOVED'),observed_free_increase_bytes=state['free_after_bytes']-state['free_before_bytes'])
 write(state);print(json.dumps({k:state[k] for k in ['removed_logical_bytes','removed_allocated_bytes','observed_free_increase_bytes','free_after_bytes']}));print('\n'.join(r['path']+': '+r['status'] for r in state['targets']))
if __name__=='__main__':main()
