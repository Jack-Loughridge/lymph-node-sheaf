"""Normal SSH only. One exact recoverable archive. Never changes permissions."""
from pathlib import Path
import hashlib,subprocess,json,datetime,os
from cleanup_verified import handles
P=Path('/home/Jack/Real_Data/spatial_ibd_sheaf/archives/GSE312415_files.tar.gz')
EXPECTED_BYTES=6464275219
EXPECTED_SHA='32fa2857563a117f974ed507083a996c34a260b52633b30fc3a9e1bc334c7cc6'
def main():
 if not P.exists():print('Already absent; no action');return
 if P.resolve()!=P or P.is_symlink() or not P.is_file():raise ValueError('Exact regular path required')
 if not os.access(P.parent,os.W_OK):raise PermissionError('Read-only boundary: stop; do not chmod or sudo')
 s=P.stat()
 if s.st_size!=EXPECTED_BYTES:raise ValueError('Size mismatch; preserve file')
 h=hashlib.sha256()
 with P.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 if h.hexdigest()!=EXPECTED_SHA:raise ValueError('Checksum mismatch; preserve file')
 inspection=handles(P)
 report=Path('/home/Jack/Real_Data/lymph_node_sheaf/data/perturbational_bioelectric_integration_2026/normal_ssh_archive_cleanup_PRIVATE.json')
 if report.exists():raise ValueError('Prior receipt exists; preserve and inspect')
 record={'UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'path':str(P),'bytes':s.st_size,'sha256':h.hexdigest(),'source_url':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE312415','status':'VERIFIED_PENDING','open_process_check':inspection}
 report.write_text(json.dumps(record,indent=2)+'\n')
 assert (P.stat().st_ino,P.stat().st_dev)==(s.st_ino,s.st_dev)
 P.unlink();record.update(status='REMOVED',absent_after=not P.exists());report.write_text(json.dumps(record,indent=2)+'\n');print('Removed exact verified archive:',s.st_size,'bytes')
if __name__=='__main__':main()
