"""Streaming transfer-integrity verification only; never loads biological arrays.
Detailed manifest/receipt must stay in ignored data/. No extraction or analysis.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path,PurePosixPath


def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def resolve_entry(root,relative):
    rel=PurePosixPath(relative)
    if rel.is_absolute() or '..' in rel.parts or '\\' in relative or not rel.parts:
        raise ValueError('Unsafe/nonportable manifest path')
    candidate=root.joinpath(*rel.parts)
    if candidate.is_symlink() or not candidate.resolve().is_relative_to(root.resolve()):
        raise ValueError('Manifest path escapes source directory')
    return candidate


def verify(root,manifest):
    entries=manifest.get('files',[])
    if not entries:raise ValueError('Empty transfer manifest')
    if len({x['relative_path'] for x in entries})!=len(entries):raise ValueError('Duplicate manifest paths')
    results=[]
    for e in entries:
        p=resolve_entry(root,e['relative_path']);reason=None
        if not p.is_file():reason='missing_file'
        elif p.stat().st_size!=int(e['bytes']):reason='byte_count_mismatch'
        elif sha256(p).lower()!=e['sha256'].lower():reason='sha256_mismatch'
        results.append(dict(relative_path=e['relative_path'],bytes=int(e['bytes']),verified=reason is None,failure=reason))
    return results


def main():
    p=argparse.ArgumentParser();p.add_argument('--source-root',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--receipt',type=Path);p.add_argument('--freeze-manifest',action='store_true');args=p.parse_args()
    # The manifest contains source filenames; preserve it privately.
    if 'data' not in args.manifest.resolve().parts:raise SystemExit('Keep source manifest under ignored data/')
    manifest=json.loads(args.manifest.read_text(encoding='utf-8-sig'));results=verify(args.source_root,manifest)
    success=all(r['verified'] for r in results)
    receipt=args.receipt or args.manifest.parent/'transfer_verification_PRIVATE.json'
    if 'data' not in receipt.resolve().parts:raise SystemExit('Keep detailed transfer receipt under ignored data/')
    if receipt.exists():raise SystemExit('Receipt exists; choose a new --receipt to preserve this checkpoint')
    with receipt.open('x') as f:json.dump(dict(manifest_sha256=sha256(args.manifest),files=results,all_verified=success),f,indent=2)
    if success and args.freeze_manifest:args.manifest.chmod(0o444)
    print(json.dumps(dict(files_expected=len(results),files_verified=sum(x['verified'] for x in results),expected_bytes=sum(x['bytes'] for x in results),all_verified=success,biological_analysis_performed=False),indent=2))
    return 0 if success else 1
if __name__=='__main__':raise SystemExit(main())
