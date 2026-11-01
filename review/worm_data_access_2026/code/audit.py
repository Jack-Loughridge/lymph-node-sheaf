"""Aggregate-only audit. Local archive metadata and recordings stay under data/.
No model fitting. Python 3 + numpy; h5py optional for actual NWB payloads.
"""
import argparse
import csv
import hashlib
import io
import json
from collections import Counter
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / 'data/worm_data_access_2026'
OUT = ROOT / 'review/worm_data_access_2026'


def content(name):
    return json.loads((CACHE / (name + '.json')).read_text())['structuredContent']['content']


def digest(data):
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(),
            'git_blob_sha1': hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()}


def summarize_manifest(assets):
    paths = [a['path'] for a in assets]
    if len(set(paths)) != len(paths):
        raise ValueError('Duplicate archive paths')
    subjects = {p.split('/')[0] for p in paths}
    seg = [a for a in assets if '_desc-segmentation_' in a['path']]
    imaging = [a for a in assets if '_desc-imaging_' in a['path']]
    sessions = {(p.split('/')[0], p.split('_ses-')[1].split('_')[0]) for p in paths}
    return dict(asset_count=len(assets), declared_subjects=len(subjects), sessions=len(sessions),
                processed_subjects=len({a['path'].split('/')[0] for a in seg}),
                imaging_assets=len(imaging), processed_assets=len(seg),
                imaging_only_subject_count=len(subjects - {a['path'].split('/')[0] for a in seg}),
                total_bytes=sum(a['size'] for a in assets),
                processed_bytes=sum(a['size'] for a in seg),
                processed_min_bytes=min((a['size'] for a in seg), default=None),
                processed_max_bytes=max((a['size'] for a in seg), default=None),
                genotype_asset_counts=dict(Counter(p.get('genotype', 'unknown') for a in assets
                    for p in a['metadata'].get('wasAttributedTo', []))),
                all_declared_open_access=all(any(x.get('status') == 'dandi:OpenAccess'
                    for x in a['metadata'].get('access', [])) for a in assets),
                all_sizes_consistent=all(a['size'] == a['metadata']['contentSize'] for a in assets))


def integrity(t, traces, starts, stops, targets, labels):
    """No interpolation, target inference, response selection or physiology fitting.
    Targets are ROI row indices; -1/-2/-3 remain unknown/excluded.
    Long gaps >3x median dt are reported, never silently repaired.
    """
    t, traces = np.asarray(t, float), np.asarray(traces, float)
    starts, stops, targets = np.asarray(starts, float), np.asarray(stops, float), np.asarray(targets)
    if t.ndim != 1 or len(t) < 2 or traces.ndim != 2 or traces.shape[0] != len(t):
        raise ValueError('Expected time x neuron traces with >=2 timestamps')
    if len(labels) != traces.shape[1] or not (len(starts) == len(stops) == len(targets)):
        raise ValueError('Trace-label or event-column length mismatch')
    if not np.all(np.isfinite(t)) or np.any(np.diff(t) <= 0):
        raise ValueError('Nonfinite or nonmonotone timestamps')
    if not np.all(np.isfinite(starts)) or not np.all(np.isfinite(stops)) or np.any(stops < starts):
        raise ValueError('Invalid stimulation times')
    dt = np.diff(t); median = float(np.median(dt))
    event_inside = (starts >= t[0]) & (stops <= t[-1])
    target_valid = np.array([float(x).is_integer() and 0 <= x < traces.shape[1] for x in targets])
    named = [str(x).strip() for x in labels if str(x).strip()]
    return dict(frames=len(t), recorded_rois=traces.shape[1], stimulation_trials=len(starts),
                duration_seconds=float(t[-1]-t[0]), median_dt_seconds=median,
                max_gap_seconds=float(dt.max()), gaps_over_3x_median=int((dt > 3*median).sum()),
                missing_trace_fraction=float(np.isnan(traces).mean()),
                infinite_trace_values=int(np.isinf(traces).sum()),
                events_outside_recording=int((~event_inside).sum()),
                known_target_trials=int(target_valid.sum()), unknown_target_trials=int((~target_valid).sum()),
                labelled_rois=len(named), unique_neuron_labels=len(set(named)),
                duplicate_nonempty_label_count=len(named)-len(set(named)),
                # A full 30 s pre-window, no earlier stimulation in it, and no long gap.
                clean_30s_prewindows=sum(bool(s-t[0] >= 30 and
                    not np.any((starts < s) & (stops > s-30)) and
                    not np.any((dt > 3*median) & (t[:-1] < s) & (t[1:] > s-30))) for s in starts))


def audit_nwb(path):
    import h5py
    with h5py.File(path, 'r') as f:
        g = f['processing/ophys/GreenSignals/BaseGreenSignal']
        st = f['intervals/OptogeneticStimulusTable']
        pp = f['processing/ophys/PumpProbeGreenSegmentations/PumpProbeGreenPlaneSegmentation']
        npal = f['processing/ophys/NeuroPALSegmentations/NeuroPALPlaneSegmentation']
        def txt(x): return x.decode() if isinstance(x, bytes) else str(x)
        mapping = dict(zip(npal['id'][:].tolist(), [txt(x) for x in npal['labels'][:]]))
        labels = [mapping.get(int(txt(x)), '') if txt(x).isdigit() else '' for x in pp['neuropal_ids'][:]]
        metrics = integrity(g['timestamps'][:], g['data'][:], st['start_time'][:],
                            st['stop_time'][:], st['target_pumpprobe_id'][:], labels)
        # Per-recording output stays local. No response thresholds or inferred targets.
        metrics['file_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
        return metrics


def run():
    pinned, current = json.loads(content('manifest_pinned')), json.loads(content('manifest'))
    summary = summarize_manifest(pinned)
    summary['pinned_current_asset_identities_equal'] = sorted((a['path'],a['size'],a['asset_id']) for a in pinned) == sorted((a['path'],a['size'],a['asset_id']) for a in current)
    summary['retrieved_recording_count'] = len(list(CACHE.glob('recordings/*.nwb')))
    (OUT/'manifest_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    # Exact file/animal inventory is intentionally ignored, never copied to review/.
    local = [{'path':a['path'], 'bytes':a['size'], 'asset_id':a['asset_id'],
              'url':a['metadata']['contentUrl'][-1],
              'sha256':a['metadata']['digest'].get('dandi:sha2-256')} for a in pinned if '_desc-segmentation_' in a['path']]
    (CACHE/'processed_inventory.json').write_text(json.dumps(local, indent=2)+'\n')
    # Choose by stable filename hash independent of outcomes; no raw imaging.
    sample = sorted(local, key=lambda a: hashlib.sha256(a['path'].encode()).hexdigest())[:3]
    (CACHE/'planned_sample.json').write_text(json.dumps(sample, indent=2)+'\n')
    refs = []
    for name in ['manifest_pinned','manifest','neuron_ids','positions','connectome','ganglia','encore_code','encore_attempt','creamer_readme','creamer_loading','creamer_create']:
        refs.append(dict(resource=name, **digest(content(name).encode())))
    with (OUT/'retrieved_resource_hashes.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(refs[0]));w.writeheader();w.writerows(refs)
    names = [line.split()[1] for line in content('neuron_ids').splitlines() if line.strip()]
    edge = list(csv.DictReader(io.StringIO(content('connectome')), delimiter='\t'))
    endpoints = {e[k] for e in edge for k in ['pre','post']}
    poslines = content('positions').splitlines(); posnames=poslines[0].lstrip('#').split()
    xyz=np.loadtxt(io.StringIO('\n'.join(poslines[1:])))
    if xyz.shape != (len(posnames),3) or not np.isfinite(xyz).all(): raise ValueError('Invalid reference positions')
    anatomy=dict(canonical_neuron_ids=len(names), unique_canonical_ids=len(set(names)),
                 connectome_rows=len(edge), endpoint_neurons=len(endpoints),
                 edge_rows_by_type=dict(Counter(e['type'] for e in edge)),
                 unmapped_endpoint_names=len(endpoints-set(names)),
                 position_rows=len(posnames), unique_position_names=len(set(posnames)),
                 duplicated_position_names=len(posnames)-len(set(posnames)),
                 position_names_outside_canonical=len(set(posnames)-set(names)),
                 physical_units_verified=False, recorded_neuron_overlap=None,
                 note='Reference coordinates and connectome are not measurements from recorded animals. AWC lateral/function aliases need explicit handling.')
    (OUT/'anatomy_summary.json').write_text(json.dumps(anatomy,indent=2)+'\n')
    inventory=[]
    for path in sorted(CACHE.glob('recordings/*.nwb')):
        r=audit_nwb(path);(CACHE/(path.stem+'.audit.json')).write_text(json.dumps(r,indent=2));inventory.append(r)
    (OUT/'recording_audit_status.json').write_text(json.dumps(dict(recordings_audited=len(inventory),
        trial_count=sum(r['stimulation_trials'] for r in inventory) if inventory else None,
        identifiable_neuron_counts=None, repeat_neuron_across_animal_counts=None,
        timestamp_integrity='NOT EVALUABLE' if not inventory else 'See ignored per-recording audits',
        reason='No recording payload retrieved; archive metadata counts are not trace integrity evidence.' if not inventory else ''),indent=2)+'\n')
    print(json.dumps(summary,indent=2)); print(json.dumps(anatomy,indent=2))

if __name__ == '__main__':
    run()
