"""Run from a normal tda-vm SSH session only; exact user-authorized archive.
This is not executed through another host to bypass the executor's read-only mount.
"""
from pathlib import Path
import json,hashlib,shutil,datetime
from cleanup_exact_downloads import check_path,open_handles,HERE,ROOT
P=Path('/home/Jack/Real_Data/spatial_ibd_sheaf/archives/GSE312415_files.tar.gz')
SHA='32fa2857563a117f974ed507083a996c34a260b52633b30fc3a9e1bc334c7cc6'
def main():
 if not P.exists():raise SystemExit('Archive absent; no action and no prior receipt overwritten')
 receipt=HERE/'archive_cleanup_host_receipt.json'
 if receipt.exists():raise SystemExit('Existing host receipt; inspect instead of overwrite')
 m=json.loads((HERE/'cleanup_manifest.json').read_text());row=next(x for x in m['targets'] if x['path']==str(P))
 if row.get('sha256')!=SHA:raise SystemExit('Recorded checksum missing or mismatched')
 s=check_path(P)
 if s.st_size!=6464275219:raise SystemExit('Unexpected archive size')
 users=open_handles(P)
 if users['pids']:raise SystemExit('Active users detected; stop')
 h=hashlib.sha256()
 with P.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 if h.hexdigest()!=SHA:raise SystemExit('Checksum changed; stop')
 after=check_path(P)
 if (s.st_ino,s.st_mtime_ns,s.st_size)!=(after.st_ino,after.st_mtime_ns,after.st_size):raise SystemExit('File changed; stop')
 users=open_handles(P)
 if users['pids']:raise SystemExit('Active users detected; stop')
 before=shutil.disk_usage(ROOT).free
 record=dict(path=str(P),bytes=s.st_size,allocated_bytes=s.st_blocks*512,sha256=SHA,source_url=row['source_url'],recoverability=row['recoverability'],open_handle_check=users,status='VERIFIED_BEFORE_REMOVAL')
 receipt.write_text(json.dumps(record,indent=2)+'\n')
 P.unlink()
 record.update(status='REMOVED',UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),free_after_bytes=shutil.disk_usage(ROOT).free)
 record['free_increase_bytes']=record['free_after_bytes']-before
 receipt.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
if __name__=='__main__':main()
