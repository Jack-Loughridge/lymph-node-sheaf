"""Recreate ignored metadata cache via ordinary permitted HTTPS access.
Only small text blobs are fetched; never recording or raw video payloads.
Requires requests. Stops on failure; does not alter networking configuration.
"""
import base64
import csv
import hashlib
import json
import requests
from audit import CACHE, OUT

BLOBS = {
 'neuron_ids':('francescorandi/wormneuroatlas','57100c3573a1a167a375c2dff82387614db344c2'),
 'positions':('francescorandi/wormneuroatlas','1e05b69bce2d1f6fd10e32f0e4a6d53d505b4445'),
 'connectome':('francescorandi/wormneuroatlas','f312c4db1018229cb5707b259440cce768e1e727'),
 'ganglia':('francescorandi/wormneuroatlas','96cfc1338dd265e50b55302827bfc3429e583cdc'),
 'encore_code':('catalystneuro/encore','95658f2c2131ce7970662edbfe8e2aca50b93701'),
 'encore_attempt':('catalystneuro/encore','b7c8f927e4856e4dfee08ed805e7f4a2044b92f6'),
 'creamer_readme':('Nondairy-Creamer/Creamer_LDS_2026','c7eb0b45a2cdc6b1f4d4e011a4c1278591dacc26'),
 'creamer_loading':('Nondairy-Creamer/Creamer_LDS_2026','c2c5a87e2e6b97fc54732678b44b82803c19023f'),
 'creamer_create':('Nondairy-Creamer/Creamer_LDS_2026','646f2fabd8711df86413f5291b70dda402314c96'),
}
COMMITS={'manifest_pinned':'674cb57b15e94e9bb17abc8b1aab8f718bec89d2',
         'manifest':'132cc50d61550e2983db1529fd1dd5ae8065546d'}

def main():
 expected={r['resource']:r for r in csv.DictReader((OUT/'retrieved_resource_hashes.csv').open())}
 CACHE.mkdir(parents=True,exist_ok=True)
 for name in [*BLOBS,*COMMITS]:
  target=CACHE/(name+'.json')
  if target.exists():
   data=json.loads(target.read_text())['structuredContent']['content'].encode()
  else:
   if name in BLOBS:
    repo,sha=BLOBS[name];url=f'https://api.github.com/repos/{repo}/git/blobs/{sha}'
    r=requests.get(url,timeout=20);r.raise_for_status();data=base64.b64decode(r.json()['content'])
   else:
    url=f'https://raw.githubusercontent.com/dandisets/001075/{COMMITS[name]}/.dandi/assets.json'
    r=requests.get(url,timeout=20);r.raise_for_status();data=r.content
  if hashlib.sha256(data).hexdigest()!=expected[name]['sha256']:raise ValueError('Source checksum mismatch: '+name)
  if not target.exists():target.write_text(json.dumps({'structuredContent':{'content':data.decode('utf-8')}}))
 print('Verified all frozen text sources; no recording data transferred.')
if __name__=='__main__':main()
