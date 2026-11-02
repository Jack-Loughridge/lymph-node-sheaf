"""Bounded anonymous processed-file acquisition after cleanup. Stops on access failure.
Never downloads raw videos; never overwrites originals. No model fitting.
"""
from pathlib import Path
import requests,json,hashlib,shutil,datetime
ROOT=Path(__file__).resolve().parents[2];HERE=Path(__file__).resolve().parent
DATA=ROOT/'data/worm_behavior_sheaf_readiness_v1/recordings';DATA.mkdir(exist_ok=True)
FILES=[('2023-03-30-01-cudata.h5',4707115,9.71),('2023-06-21-01-cudata.h5',4707119,11.38),('2023-06-23-08-cudata.h5',4707120,10.17),('NeuroPAL_labels.jld2',4707122,41.61),('load_data_example.jl',4707124,.00161)]

def verify_format(name,p):
 b=p.open('rb').read(8192)
 if name.endswith(('.h5','.jld2')) and not any(b[k:k+8]==b'\x89HDF\r\n\x1a\n' for k in [0,512,1024,2048,4096]):raise ValueError('Not HDF5/JLD2; possible error page')
 if name.endswith('.jl'):
  if b.lstrip().startswith(b'<'):raise ValueError('Error page instead of Julia source')
  b.decode('utf-8-sig')
 return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
 cleanup=json.loads((HERE/'cleanup_manifest.json').read_text())
 if not cleanup.get('finished_UTC'):raise SystemExit('Cleanup must finish before acquisition')
 if shutil.disk_usage(ROOT).free<1_000_000_000:raise SystemExit('Insufficient reserve')
 rows=[];blocked=False
 for name,fileid,mb in FILES:
  row=dict(file=name,url=f'https://datadryad.org/downloads/file_stream/{fileid}',listed_MB_rounded=mb,repository_checksum=None);p=DATA/name
  if p.exists():
   row.update(status='LOCAL_PRESENT_PENDING_SCIENTIFIC_QC',bytes=p.stat().st_size,sha256=verify_format(name,p));rows.append(row);continue
  if blocked:row['status']='NOT_RETRIED_SAME_HOST_BLOCKED';rows.append(row);continue
  part=p.with_suffix(p.suffix+'.part')
  if part.exists():raise SystemExit('Preexisting partial file; inspect instead of overwrite')
  try:
   response=requests.get(row['url'],timeout=(5,15),stream=True);response.raise_for_status()
   cap=50_000_000 if name.endswith('.jld2') else 15_000_000 if name.endswith('.h5') else 100_000
   n=0
   with part.open('xb') as f:
    for chunk in response.iter_content(65536):
     n+=len(chunk)
     if n>cap:raise ValueError('Payload exceeds processed-file bound')
     f.write(chunk)
   sha=verify_format(name,part);part.rename(p);row.update(status='ACQUIRED_PENDING_SCIENTIFIC_QC',bytes=n,sha256=sha)
  except (requests.RequestException,OSError,ValueError) as e:
   row.update(status='ACCESS_BLOCKED',error=str(e)[:800]);blocked=True
   # Keep any downloaded partial bytes privately for diagnosis, never call them data.
  rows.append(row)
 result=dict(time_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),files=rows,acquired_recordings=sum(r['file'].endswith('.h5') and r['status'].startswith(('ACQUIRED','LOCAL')) for r in rows),all27_download_authorized_only_after_sample_QC=True,network_retries=0)
 (HERE/'acquisition_inventory.json').write_text(json.dumps(result,indent=2)+'\n');print('\n'.join(r['file']+': '+r['status'] for r in rows))
if __name__=='__main__':main()
