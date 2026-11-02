"""Bounded small-file acquisition. No spatial-atlas downloads or fitted analyses."""
from pathlib import Path
import requests,json,hashlib,gzip,io
ROOT=Path(__file__).resolve().parents[2]
DEST=ROOT/'data/crohn_xenium/supporting_files';DEST.mkdir(exist_ok=True)
R=Path(__file__).resolve().parent
PREFIX='https://static-content.springer-cdn.com/esm/art%3A10.1038%2Fs41590-026-02617-0/MediaObjects/'
FILES={
 'GSE307526_processed_bulkRNAseq_cytokine_stimulations.csv.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE307nnn/GSE307526/suppl/',
 '41590_2026_2617_MOESM3_ESM.pdf':PREFIX,
 '41590_2026_2617_MOESM6_ESM.csv':PREFIX,
 '41590_2026_2617_MOESM9_ESM.csv':PREFIX,
 '41590_2026_2617_MOESM12_ESM.csv':PREFIX}
def validate(name,b):
 if b.lstrip().lower().startswith((b'<html',b'<!doctype',b'<?xml')):raise ValueError('Error/document payload instead of data')
 if name.endswith('.gz'):
  with gzip.GzipFile(fileobj=io.BytesIO(b)) as f: b=f.read(1024)
 if name.endswith('.pdf') and not b.startswith(b'%PDF'):raise ValueError('Not PDF')
 if name.endswith('.csv') and b',' not in b[:8192]:raise ValueError('Not comma-delimited content')
 return True

def main():
 rows=[]
 for name,prefix in FILES.items():
  p=DEST/name;row={'file':name,'url':prefix+name}
  if p.exists():
   validate(name,p.read_bytes());row.update(status='LOCAL_PRESENT_NOT_YET_SCIENTIFICALLY_VERIFIED',bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest());rows.append(row);continue
  try:
   response=requests.get(prefix+name,timeout=(5,12),stream=True);response.raise_for_status();chunks=[];n=0
   for b in response.iter_content(65536):
    n+=len(b)
    if n>6_000_000:raise ValueError('Above small-file size cap')
    chunks.append(b)
   b=b''.join(chunks);validate(name,b);p.write_bytes(b);row.update(status='RETRIEVED_NOT_YET_SCIENTIFICALLY_VERIFIED',bytes=n,sha256=hashlib.sha256(b).hexdigest())
  except (requests.RequestException,ValueError,OSError) as e:row.update(status='BLOCKED',error=str(e)[:700])
  rows.append(row)
 (R/'supporting_file_access.json').write_text(json.dumps(rows,indent=2)+'\n')
 print('\n'.join(x['file']+': '+x['status'] for x in rows))
if __name__=='__main__':main()
