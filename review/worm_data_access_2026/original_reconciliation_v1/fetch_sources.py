"""Restore pinned author code to ignored cache; not executed as a pipeline.
Ordinary HTTPS only. Stops on failure; no retries or network configuration changes.
"""
import csv,hashlib
from pathlib import Path
import urllib.request

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CACHE=ROOT/'data/worm_data_access_2026/original_reconciliation_v1'

def main():
 CACHE.mkdir(parents=True,exist_ok=True)
 for row in csv.DictReader((HERE/'source_hashes.csv').open()):
  p=CACHE/row['local_filename']
  if p.exists():data=p.read_bytes()
  else:
   with urllib.request.urlopen(row['url'],timeout=20) as r:data=r.read(1_000_001)
  if len(data)>1_000_000 or hashlib.sha256(data).hexdigest()!=row['sha256']:raise ValueError('Pinned source mismatch: '+row['local_filename'])
  if not p.exists():p.write_bytes(data)
 print('Pinned author sources verified; no source pipeline, classifier or pickle executed.')
if __name__=='__main__':main()
