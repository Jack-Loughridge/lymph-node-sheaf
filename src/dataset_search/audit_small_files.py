"""Audit retrieved source data without claiming cell-level coverage or prediction.

Run from repository root. Input JSON envelopes are exact GitHub tool downloads
under ignored data/. Every decoded byte string is checked against its Git blob.
Outputs contain aggregates only. Requires existing pandas and numpy.
"""
import ast
import hashlib
import io
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd

EXPECTED = {
    'plan_fig1c': '337073e4791b9df9099d12a9a162aa352a6c1431',
    'plan_fig1d': 'ae8609937cbabc5b125c37dc46ea9be7ddb422f6',
    'plan_fig4j': 'ce961c31b31524d57df32dfdc5e961402ce5eaf3',
    'worm_gap_expression': '2ccf9f9eb6c74cace0c14645be24aa0112bf929b',
    'worm_connectome': 'f312c4db1018229cb5707b259440cce768e1e727',
    'crohn_annotation': '8ed30965477b4fb06da117ec96f613e5933e629b',
    'crohn_response_a': '6a3c136fc2090ce5f9d099677d4533b63da8c5c9',
    'crohn_response_b': 'f8ccdcd19470ab4b29170c228c0c95e7c009d4e7',
    'age_metadata_code': '5da32b8233421d862c5ddce1f51e03675ae04d5f',
    'sp_vignette': '66a3c1b56ae2813688e4a8b80bafcd44ef93c86e',
}


def verify_blob(text, expected):
    payload = text.encode('utf-8')
    actual = hashlib.sha1(b'blob ' + str(len(payload)).encode() + b'\0' + payload).hexdigest()
    if actual != expected:
        raise ValueError(f'Git blob mismatch: {actual} != {expected}')
    return payload


def notebook_source(text):
    return '\n'.join(''.join(c.get('source', [])) for c in json.loads(text)['cells'])


def main(output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    inputs = Path('data/dataset_search_2026')
    texts, manifest = {}, []
    for name, expected in EXPECTED.items():
        envelope = json.loads((inputs / f'{name}.json').read_text())
        text = envelope['structuredContent']['content']
        payload = verify_blob(text, expected)
        texts[name] = text
        manifest.append({'name': name, 'git_blob': expected, 'bytes': len(payload),
                         'sha256': hashlib.sha256(payload).hexdigest(),
                         'source_url': envelope['structuredContent'].get('url')})

    # Figure source values are animal-level summaries, not expression matrices.
    plan = pd.read_csv(io.StringIO(texts['plan_fig1c']), sep='\t')
    plan_summary = plan.groupby('timepoint', sort=False).agg(
        source_rows=('timepoint', 'size'),
        segmented_units_qc50=('cells with UMI greater than 50', 'sum'))
    plan_summary.to_csv(output / 'planarian_source_summary.csv', mode='x')
    assert len(plan) == 16 and len(plan_summary) == 8
    assert plan_summary.source_rows.eq(2).all()
    response = []
    programme = None
    for line in texts['plan_fig4j'].splitlines():
        if line.startswith(('smed03831+', 'sfrp1+')):
            programme = line.split('+')[0]
        if re.match(r'^(1\.5|3)dpa ', line):
            fields = line.split('\t')
            values = [float(v) for v in fields[1:5] if v.strip()]
            response.append({'marker': programme, 'condition': fields[0],
                             'reported_replicate_values': len(values),
                             'mean_reported_normalized_value': float(np.mean(values)),
                             'warning': 'FISH aggregate source values; not ten-gene response map'})
    pd.DataFrame(response).to_csv(output / 'planarian_perturbation_summary.csv', index=False, mode='x')

    # CeNGeN averages are cell-class reference expression, NOT single-animal
    # measurements, detection percentages, or physiological exposure.
    expression = pd.read_csv(io.StringIO(texts['worm_gap_expression']))
    reference_types = expression.columns[4:]
    rows = []
    for _, row in expression.iterrows():
        values = row[reference_types].astype(float)
        assert np.isfinite(values).all() and (values >= 0).all()
        rows.append({'gene': row['gene_name'], 'reference_classes': len(values),
                     'classes_nonzero_reference_average': int((values > 0).sum()),
                     'not_a_cell_detection_or_exposure_count': True})
    pd.DataFrame(rows).to_csv(output / 'worm_gap_reference_coverage.csv', index=False, mode='x')
    edges = pd.read_csv(io.StringIO(texts['worm_connectome']), sep='\t')
    assert (edges.synapses > 0).all()
    edge_summary = edges.groupby('type').agg(
        rows=('synapses', 'size'), contacts=('synapses', 'sum'),
        source_ids=('pre', 'nunique'), receiver_ids=('post', 'nunique'))
    edge_summary.to_csv(output / 'worm_anatomical_edge_summary.csv', mode='x')

    # Inspect the recorded data schema, not fabricated contents of the .h5ad.
    annotation = json.loads(texts['crohn_annotation'])
    annotation_outputs = '\n'.join(
        ''.join(o.get('data', {}).get('text/plain', []))
        for c in annotation['cells'] for o in c.get('outputs', []))
    match = re.search(r'n_obs × n_vars = (\d+) × (\d+)', annotation_outputs)
    schema = {'recorded_object_rows': int(match[1]), 'recorded_object_features': int(match[2]),
              'recorded_counts_layer': "layers: 'counts'" in annotation_outputs,
              'recorded_spatial_coordinates': "'spatial'" in annotation_outputs,
              'recorded_author_annotations': "'annot2'" in annotation_outputs,
              'loaded_matrix': False, 'actual_detection_and_exposure': None}
    age_source = notebook_source(texts['age_metadata_code'])
    parsed = ast.parse(age_source)
    samples = next(ast.literal_eval(n.value) for n in parsed.body if isinstance(n, ast.Assign)
                   and any(isinstance(t, ast.Name) and t.id == 'metas' for t in n.targets))
    age_groups = pd.DataFrame([{'age_months': int(s.split('-')[0]),
                                'sex': s.split('-')[2], 'source_sample_entries': 1}
                               for s in samples]).groupby(['age_months', 'sex']).sum()
    age_groups.to_csv(output / 'age_source_sample_summary.csv', mode='x')
    assert len(samples) == len(set(samples)) == 12
    results = {
        'status': 'completed small-file/source-schema audit; biological coverage NOT EVALUABLE',
        'planarian_timepoints': len(plan_summary), 'planarian_source_rows': len(plan),
        'planarian_minimum_rows_per_timepoint': int(plan_summary.source_rows.min()),
        'planarian_total_segmented_units_in_source': int(plan.iloc[:, 1].sum()),
        'worm_gap_reference_genes': len(rows), 'worm_reference_classes': len(reference_types),
        'crohn_schema': schema, 'age_distinct_source_sample_entries': len(samples),
        'no_new_cell_expression_matrix_loaded': True,
        'exposed_receiver_counts': None, 'held_out_mse_improvement': None,
        'interpretation': 'Source tables are real downloads with verified hashes. They do not substitute for matched spatial matrices, independent response maps or donor-held-out prediction.',
        'inputs': manifest,
    }
    with (output / 'SMALL_FILE_AUDIT.json').open('x') as handle:
        json.dump(results, handle, indent=2)
        handle.write('\n')
    print(json.dumps({k: v for k, v in results.items() if k != 'inputs'}, indent=2))


if __name__ == '__main__':
    main(sys.argv[1])
