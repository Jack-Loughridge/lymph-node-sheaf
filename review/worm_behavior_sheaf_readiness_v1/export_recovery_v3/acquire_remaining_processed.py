"""Normal-VM acquisition of ONLY the 24 manifest-listed processed HDF5 files.
Discover actual links from public HTML; no guessed file IDs, raw imaging or retries.
Run --download only after the three-sample processed-signal QC review.
"""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urljoin,urlparse
import argparse,csv,json,hashlib,shutil,datetime
import requests
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
DEST=ROOT/'data/worm_behavior_sheaf_readiness_v1/recordings'
JOURNAL=ROOT/'data/worm_behavior_sheaf_readiness_v1/export_recovery_v3/remaining_acquisition_PRIVATE.json'
STUDY='https://datadryad.org/dataset/doi:10.5061/dryad.w9ghx3g4v'
class FileLinks(HTMLParser):
 def __init__(self):super().__init__();self.current=None;self.links=[]
 def handle_starttag(self,tag,attrs):
  if tag=='a':self.current=[dict(attrs).get('href',''),[]]
 def handle_data(self,text):
  if self.current:self.current[1].append(text)
 def handle_endtag(self,tag):
  if tag=='a' and self.current:self.links.append((self.current[0],''.join(self.current[1]).strip()));self.current=None

def discover_links(html,names):
 p=FileLinks();p.feed(html);found={}
 for href,text in p.links:
  for name in names:
   if name in text:
    url=urljoin(STUDY,href)
    if urlparse(url).hostname!='datadryad.org' or '/downloads/file_stream/' not in urlparse(url).path:continue
    if name in found and found[name]!=url:raise ValueError('Ambiguous public download link: '+name)
    found[name]=url
 if set(found)!=set(names):raise ValueError('Public HTML does not expose all named processed downloads; use browser transfer, do not guess IDs')
 return found

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--download',action='store_true');args=ap.parse_args()
 manifest=list(csv.DictReader((HERE/'REMAINING_FILES.csv').open()))
 names=[r['file'] for r in manifest]
 if len(names)!=24 or len(set(names))!=24 or any(Path(x).name!=x or not x.endswith('-cudata.h5') for x in names):raise ValueError('Unexpected allowlist')
 if not (ROOT/'data/worm_behavior_sheaf_readiness_v1/sample_ingestion_v2/native_identity/native_receipt_PRIVATE.toml').is_file():raise ValueError('Native validation not completed')
 if shutil.disk_usage(ROOT).free<1_000_000_000:raise ValueError('Insufficient reserve')
 if JOURNAL.exists():raise ValueError('Existing acquisition journal: preserve and review before retry')
 response=requests.get(STUDY,timeout=(5,15));response.raise_for_status()
 urls=discover_links(response.text,names)
 state={'start_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':'Full processed-cohort QC; genuine-observation gate remains unresolved','files':[],'repository_checksums_available':False}
 JOURNAL.write_text(json.dumps(state,indent=2)+'\n')
 for row in manifest:
  name=row['file'];p=DEST/name;item={'file':name,'url':urls[name],'status':'PLAN_ONLY'};state['files'].append(item)
  if p.exists():
   item.update(status='PRESENT_PRESERVED',bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
  elif args.download:
   part=p.with_suffix(p.suffix+'.part')
   if part.exists():raise ValueError('Partial exists; inspect before retry')
   try:
    r=requests.get(urls[name],stream=True,timeout=(5,25));r.raise_for_status();n=0;h=hashlib.sha256()
    with part.open('xb') as out:
     for chunk in r.iter_content(65536):
      n+=len(chunk)
      if n>20_000_000:raise ValueError('Processed-file size cap exceeded')
      h.update(chunk);out.write(chunk)
    with part.open('rb') as inp:magic=inp.read(8)
    if magic!=b'\x89HDF\r\n\x1a\n':raise ValueError('Downloaded object is not HDF5')
    if not .8*float(row['listed_size_MB_rounded'])*1_000_000<=n<=1.2*float(row['listed_size_MB_rounded'])*1_000_000:raise ValueError('Unexpected manifest size discrepancy')
    part.rename(p);item.update(status='DOWNLOADED_PENDING_SCIENTIFIC_QC',bytes=n,sha256=h.hexdigest())
   except Exception as error:
    item.update(status='FAILED_STOPPED_NO_RETRY',error=str(error)[:500]);JOURNAL.write_text(json.dumps(state,indent=2)+'\n');raise
  JOURNAL.write_text(json.dumps(state,indent=2)+'\n')
 print('Processed files planned/checked:',len(names),'; downloaded:',sum(x['status'].startswith('DOWNLOADED') for x in state['files']))
if __name__=='__main__':main()
