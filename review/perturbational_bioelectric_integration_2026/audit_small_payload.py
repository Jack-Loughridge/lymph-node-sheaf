"""Budgeted private ZIP/XLSX/MAT schema audit. Never declares biological GO.
No untrusted code is executed and no full imaging collection is extracted.
"""
from pathlib import Path,PurePosixPath
from io import BytesIO
import zipfile,stat,hashlib,json,argparse,re,xml.etree.ElementTree as ET
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
PRIVATE=ROOT/'data/perturbational_bioelectric_integration_2026'
LIMIT=50_000_000;MEMBER=20_000_000

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()

def manifest(z):
 infos=z.infolist();seen=set()
 if len(infos)>5000 or sum(i.file_size for i in infos)>150_000_000:raise ValueError('Uncompressed archive budget exceeded; inspect manifest externally first')
 for i in infos:
  p=PurePosixPath(i.filename);mode=i.external_attr>>16
  if p.is_absolute() or '..' in p.parts or '\\' in i.filename or ':' in i.filename:raise ValueError('Unsafe path')
  if i.filename in seen:raise ValueError('Duplicate archive member')
  if stat.S_ISLNK(mode) or stat.S_IFMT(mode) not in (0,stat.S_IFREG,stat.S_IFDIR):raise ValueError('Nonregular member')
  if i.flag_bits&1:raise ValueError('Encrypted archive')
  seen.add(i.filename)
 return infos

N={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
def xlsx_schema(data):
 """OOXML structural interpretation, not ad-hoc biological identity decoding."""
 with zipfile.ZipFile(BytesIO(data)) as z:
  infos=manifest(z)
  if any(i.file_size>MEMBER for i in infos):raise ValueError('Large nested workbook XML; stop rather than allocate')
  strings=[]
  if 'xl/sharedStrings.xml' in z.namelist():
   root=ET.fromstring(z.read('xl/sharedStrings.xml'))
   strings=[''.join(e.itertext()) for e in root.findall('m:si',N)]
  result=[]
  for i in infos:
   if re.fullmatch(r'xl/worksheets/sheet\d+\.xml',i.filename):
    root=ET.fromstring(z.read(i));rows=root.findall('m:sheetData/m:row',N);headers=[];formula_count=0;numeric=0
    for row in rows:
     for cell in row.findall('m:c',N):
      if cell.find('m:f',N) is not None:formula_count+=1
      value=cell.find('m:v',N);kind=cell.get('t','n')
      if kind=='s' and value is not None:
       idx=int(value.text);text=strings[idx];headers.append(text)
      elif kind=='inlineStr':headers.append(''.join(cell.find('m:is',N).itertext()))
      elif value is not None and kind=='n':numeric+=1
    explicit_time=any(re.search(r'\btime\b|timestamp|second|frame|millisecond',x,re.I) for x in headers)
    result.append({'worksheet_XML':i.filename,'row_count':len(rows),'numeric_cell_count':numeric,'formula_cells':formula_count,'text_labels_PRIVATE':headers,'potential_time_labels':explicit_time,'not_proof_of_simultaneous_cells':True})
  return result

def member_schema(name,data):
 suffix=Path(name).suffix.lower()
 if suffix=='.xlsx':return {'kind':'XLSX_OOXML','worksheets':xlsx_schema(data)}
 if suffix in ['.csv','.tsv','.tab','.txt']:
  text=data.decode('utf-8-sig',errors='strict');lines=text.splitlines()
  return {'kind':'TEXT_TABLE','rows':len(lines),'first_lines_PRIVATE':lines[:8],'time_units_and_group_semantics_verified':False}
 if suffix in ['.mat','.h5','.hdf5']:
  import h5py
  try:
   with h5py.File(BytesIO(data),'r') as f:
    ds={};f.visititems(lambda k,v:ds.update({k:{'shape':list(v.shape),'dtype':str(v.dtype)}}) if isinstance(v,h5py.Dataset) else None)
    return {'kind':'HDF5_OR_MAT73','datasets':ds,'physiological_semantics_verified':False}
  except OSError:
   from scipy.io import whosmat
   return {'kind':'MAT_LEVEL5','variables':[{'name':n,'shape':list(shape),'class':cls} for n,shape,cls in whosmat(BytesIO(data))],'nested_cells_need_explicit_author_schema':True}
 return {'kind':'UNSUPPORTED_NOT_DECODED','extension':suffix,'bytes':len(data),'reader_required_before_GO':True}

def audit(input_path,expected_sha=None):
 p=Path(input_path)
 if p.is_symlink() or not p.is_file():raise ValueError('Expected regular local file')
 if p.stat().st_size>LIMIT:raise ValueError('Initial50MB acquisition budget exceeded; no analysis')
 h=sha(p)
 if expected_sha and h.lower()!=expected_sha.lower():raise ValueError('Sender SHA256 mismatch')
 result={'source_name_PRIVATE':p.name,'sha256':h,'bytes':p.stat().st_size,'sender_hash_checked':bool(expected_sha),'members':[],'biological_GO':False}
 if zipfile.is_zipfile(p) and p.suffix.lower()!='.xlsx':
  with zipfile.ZipFile(p) as z:
   infos=manifest(z);read=0
   for i in infos:
    entry={'name_PRIVATE':i.filename,'bytes':i.file_size,'compressed_bytes':i.compress_size,'CRC32':f'{i.CRC:08x}'}
    if not i.is_dir() and i.file_size<=MEMBER and read+i.file_size<=LIMIT:
     data=z.read(i);read+=len(data);entry['CRC_verified']=True;entry['sha256']=hashlib.sha256(data).hexdigest()
     try:entry['schema']=member_schema(i.filename,data)
     except Exception as e:entry['schema_status']='READER_FAILED_STOP_NO_REINTERPRETATION';entry['reason']=str(e)[:300]
    else:entry['schema_status']='SKIPPED_SIZE_OR_DIRECTORY'
    result['members'].append(entry)
 else:result['members'].append({'name_PRIVATE':p.name,'schema':member_schema(p.name,p.read_bytes())})
 return result

def main():
 ap=argparse.ArgumentParser();ap.add_argument('input',type=Path);ap.add_argument('--sha256');ap.add_argument('--tag',default='representative');a=ap.parse_args()
 if not re.fullmatch(r'[A-Za-z0-9_-]+',a.tag):raise ValueError('Unsafe tag')
 dest=PRIVATE/(a.tag+'_payload_audit_PRIVATE.json')
 if dest.exists():raise ValueError('Preserve prior audit; use new tag')
 result=audit(a.input,a.sha256);dest.write_text(json.dumps(result,indent=2)+'\n')
 summary={'tag':a.tag,'bytes':result['bytes'],'sha256':result['sha256'],'source_members':len(result['members']),'sender_hash_checked':result['sender_hash_checked'],'biological_GO':False,'status':'SCHEMA_INSPECTED_REQUIRES_GROUPING_CLOCK_SPATIAL_AND_PERTURBATION_REVIEW'}
 (HERE/(a.tag+'_AUDIT_AGGREGATE.json')).write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
