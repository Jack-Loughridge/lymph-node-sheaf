"""Bounded plain-text export reconciliation; no pickle loading or tar extraction.
Optional: used by reconcile.py --wt-archive PATH after manual transfer.
"""
import csv,hashlib,io,json,re,tarfile
from pathlib import Path,PurePosixPath
import numpy as np

MAX_MEMBER=128*1024*1024
MAX_MEMBERS=10000
KNOWN_VERSION1_MD5={'exported_data.tar.gz':'7ae518fea414f664e2493238afe209bb','exported_data_unc31.tar.gz':'e5c83242baa131bca95ee58bc7b58bec'}


def session_keys(name):
    return set('_'.join(m) for m in re.findall(r'(?<!\d)(20\d{6})[_-](\d{6})(?!\d)',name))


def safe_member(member):
    p=PurePosixPath(member.name)
    if p.is_absolute() or '..' in p.parts or member.issym() or member.islnk():raise ValueError('Unsafe archive path/link')
    if member.size>MAX_MEMBER:raise ValueError('Member exceeds text audit budget')
    return member.isfile()


def read_text(tar,member):
    if not safe_member(member):raise ValueError('Not a regular file')
    return tar.extractfile(member).read(MAX_MEMBER+1).decode('utf-8')


def hashfile(path,alg):
    h=hashlib.new(alg)
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def reconcile_text_archive(path,event_rows,manifest,out):
    path=Path(path);md5=hashfile(path,'md5')
    expected=KNOWN_VERSION1_MD5.get(path.name)
    if expected and expected!=md5:raise ValueError('Reported version-1 MD5 differs; resolve archive version before reconciling')
    archive=dict(path=str(path),bytes=path.stat().st_size,md5=md5,sha256=hashfile(path,'sha256'),reported_version1_md5_matched=bool(expected and expected==md5))
    matches={};inventory=[]
    with tarfile.open(path,'r:gz') as tar:
        members=tar.getmembers()
        if len(members)>MAX_MEMBERS:raise ValueError('Archive file count exceeds budget')
        lookup={m.name:m for m in members if m.isfile()}
        if len(lookup)!=sum(m.isfile() for m in members):raise ValueError('Duplicate archive member names')
        for m in members:
            # Index filenames/sizes, but never extract or execute any payload.
            inventory.append(dict(name=m.name,size=m.size,regular=m.isfile()))
            if m.isfile() and m.name.endswith('_ds_name.txt'):
                ds=read_text(tar,m).strip();keys=session_keys(ds)
                if len(keys)!=1:continue
                matches.setdefault(next(iter(keys)),[]).append(m.name[:-len('ds_name.txt')])
        for entry in manifest:
            key=entry['filename'].split('_ses-')[1].split('_desc')[0].replace('-','_')
            rows=[r for r in event_rows if r['subject']==entry['subject']];prefixes=matches.get(key,[])
            if len(prefixes)!=1:
                for row in rows:row['original_export_session_match']='NOT_FOUND' if not prefixes else 'AMBIGUOUS_SESSION'
                continue
            prefix=prefixes[0];need=['gcamp.txt','t.txt','labels.txt','stim_volume_i.txt','stim_neurons.txt']
            if any(prefix+n not in lookup for n in need):
                for row in rows:row['original_export_session_match']='SESSION_FOUND_REQUIRED_TEXT_FILES_MISSING'
                continue
            loaded={n:read_text(tar,lookup[prefix+n]) for n in need}
            labels=loaded['labels.txt'].splitlines();g=np.loadtxt(io.StringIO(loaded['gcamp.txt']),ndmin=2)
            tt=np.atleast_1d(np.loadtxt(io.StringIO(loaded['t.txt'])));vi=np.atleast_1d(np.loadtxt(io.StringIO(loaded['stim_volume_i.txt'])));si=np.atleast_1d(np.loadtxt(io.StringIO(loaded['stim_neurons.txt'])))
            if g.shape!=(len(tt),len(labels)) or len(si)!=len(vi):raise ValueError('Original export matrix/metadata dimensions disagree')
            for row in rows:
                row['original_export_session_match']='MATCHED_STABLE_SESSION'
                i=row['event_index']
                if len(si)!=row['nwb_total_events'] or len(tt)!=row['nwb_total_frames']:
                    row['original_export_session_match']='MATCHED_SESSION_EVENT_OR_FRAME_COUNTS_DIFFER';continue
                if i>=len(vi) or not np.isfinite(vi[i]) or not float(vi[i]).is_integer() or abs(int(vi[i])-row['nearest_nwb_stimulus_frame'])>1:
                    row['original_export_session_match']='MATCHED_SESSION_EVENT_ALIGNMENT_UNRESOLVED';continue
                raw=float(si[i])
                row['original_export_target_code']=raw
                row['original_export_exclusion']='EXPLICIT_CODE_MINUS2' if raw==-2 else 'UNRESOLVED_CODE_MINUS1' if raw==-1 else 'COMPLEMENTARY_LABEL_MINUS3' if raw==-3 else 'NOT_MARKED_EXCLUDED' if raw>=0 else 'OTHER_CODE'
                if raw==-2:row['original_author_exclusion']='EXCLUDED_IN_THIS_EXPORT'
                elif raw>=0:row['original_author_exclusion']='NOT_MARKED_EXCLUDED_IN_THIS_EXPORT; FULL_INCLUSION_UNKNOWN'
                # Match labelled source columns independently of numeric ROI ordinal.
                hit=[j for j,label in enumerate(labels) if label.strip()==row['target_name']]
                if len(hit)!=1:
                    row['original_target_trace_status']='LABEL_JOIN_AMBIGUOUS_OR_MISSING';continue
                j=hit[0];row['original_target_trace_status']='MATCHED_LABEL; MEASURED_VS_IMPUTED_UNRESOLVED_WITHOUT_MASK'
                row['original_target_trace_finite_fraction']=float(np.isfinite(g[:,j]).mean())
                if raw>=0 and float(raw).is_integer() and raw<len(labels):row['original_target_label_agrees']=labels[int(raw)].strip()==row['target_name']
                # Plain text export contains no response-list flag or original nan mask.
                row['original_autoresponse_flag']='NOT_PRESENT_IN_FLAT_TEXT_EXPORT'
                row['original_published_inclusion']='NOT_EVALUABLE_WITHOUT_ORIGINAL_RESPONSE_FLAGS_AND_MASK'
    (out/'original_export_inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
    (out/'original_archive_integrity.json').write_text(json.dumps(archive,indent=2)+'\n')
    # Ensure CSV schema covers optional fields even when some rows could not be matched.
    fields={k for r in event_rows for k in r}
    for r in event_rows:
        for key in fields:r.setdefault(key,None)
    return 'ARCHIVE_INSPECTED; SEE_PER_EVENT_MATCH_STATUS; AUTORESPONSE_STILL_REQUIRES_ORIGINAL_FLAGS'
