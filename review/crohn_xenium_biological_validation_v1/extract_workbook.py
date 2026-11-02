"""Run with system Python/openpyxl; save only a private neutral-format cache."""
import json,hashlib,re
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
source=ROOT/'data/crohn_xenium/261023_Xenium_metadata_linkage.xlsx'
out=ROOT/'data/crohn_xenium/biological_validation_v1';out.mkdir(exist_ok=True)
df=pd.read_excel(source,dtype=str,keep_default_na=False,engine='openpyxl')
df=df[df['Table S1 patient ID'].str.fullmatch(r'\d+(?:\.0)?')]
payload=dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),rows=df.to_dict('records'))
(out/'workbook_linkage_PRIVATE.json').write_text(json.dumps(payload,indent=2)+'\n')
print('Private linkage cache refreshed:',len(df),'patient rows; no identifiers printed.')
