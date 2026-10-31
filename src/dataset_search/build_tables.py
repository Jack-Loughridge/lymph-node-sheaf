"""Render the curated evidence catalogue and conservative acceptance ledger.

Judgments are explicitly curated feasibility judgments, not statistical tests.
Reproduce into a NEW destination: python3 src/dataset_search/build_tables.py DIR
"""
import csv
import json
from pathlib import Path
import sys


ROOT = Path('review/dataset_search_2026')
CRITERIA = {
    'A1': '>=95% unique QC expression/coordinate joins; cell-resolved observations',
    'A2': 'physical units and boundaries verified',
    'A3_broad': '>=90% credible broad identities',
    'A3_receiver': '>=70% supported receiver subtypes',
    'A4': 'independent marker/protein/anatomy concordance',
    'B1': '>=2 distinct experimentally supported programmes; no feature reuse',
    'B2_panel': '>=10 informative non-input response variables per programme',
    'B2_detection': '>=5 targets detected in >=20% correctly typed receivers',
    'B3': 'independent perturbation reproducibility',
    'C1': '>=100 exposed/route, >=20% appropriate receivers in >=3 specimens',
    'C2': 'signal above prespecified background; dropout and edges audited',
    'C3': 'plausible spatial ranges; complete mechanisms; positive/negative controls',
    'D1': 'frozen external maps; specimen-held-out training-only calibration',
    'D2_intrinsic': '>=5% held-out MSE gain versus strong target-excluded intrinsic baseline',
    'D2_null': '>=5% held-out MSE gain versus geometry/type-preserving null',
    'D3': '>=2 programmes with animal/donor 95% paired interval excluding zero',
    'D4': 'exposed/all receivers; boundary and missingness performance',
    'E1': 'effect direction agrees in >=2/3 plausible radii',
    'E2': '80% direction stability under 20% downsampling; median >=half effect',
    'E3': 'multiple QC/detection floors and composition/geometry matched nulls',
    'F': '>=8 independent biological units and >=3/group for intended comparison',
}

EVIDENCE = {
    'CD5K': 'Downloaded author notebook records 181597 x 5001 mesenchymal object with counts, spatial and annot2; no matrix loaded. GEO lists 44 culture libraries/4 replicates per treatment; independent tissue/culture mapping unresolved.',
    'DSS': 'Dryad lists 940-gene colon resource and micron coordinate schema; 15 mice reported; per-mouse group map, gene list and matrix not transferred. Lymph-node Immune Dictionary cannot establish colon transfer.',
    'PFISH': 'Original scripts and registered multimodal design inspected; spatial count/perturbation tables unavailable through current portal; no independent culture/animal table retrieved.',
    'WORM': 'Exact downloaded anatomical reference: 2204 chemical and 296 electrical edge rows; 21 gap genes x 128 reference classes. Reference expression is not cell-level detection or exposure. Calcium response payload unavailable.',
    'PLAN': 'Exact downloaded source table: 16 rows at 8 stages, 2 per stage; 3508004 QC segmented units in table. FISH source has only two named markers and 2-4 reported replicate values/condition; not a ten-gene response matrix.',
    'SPSEQ': 'GEO lists 3 spatial mouse samples versus 4 in paper. 5 dissociated animals do not expand spatial n. Example RDS 69403001 bytes listed but not transferred/read; vignette is code, not measured response data.',
    'AGE500': 'Downloaded source code has 12 sample entries: age3/24 each 3 female +3 male, corroborating paper. 500-gene counts and segmentation not loaded; microglia-only object cannot supply all senders.',
    'SPAC': 'Original scripts and GSA raw accession verified; no usable processed matrix manifest. Specific spatial screen 3 mice; 5-mouse bulk experiment cannot count as spatial replication.',
    'ZFHEART': 'Measured atlas uses bin50/bin70 about35um spatial units; 95 hearts total across modalities is not verified spatial n. Zenodo explicitly separates imputed from measured objects.',
    'AXOLIMB': 'GEO lists six timepoint samples and 100 measured genes. 14dpa matrix 1.7MB, cell table808.8KB, panel3.4KB; HEAD requests fail DNS; gene overlap/detection unavailable.',
}


def write_csv(path, rows):
    with path.open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    catalogue = json.loads((ROOT / 'catalogue.json').read_text())
    assert len(catalogue) >= 12
    assert len({r['id'] for r in catalogue}) == len(catalogue)
    assert [r['rank'] for r in catalogue] == list(range(1, len(catalogue) + 1))
    write_csv(out / 'ranked_candidates.csv', catalogue)
    ledger = []
    for candidate, evidence in EVIDENCE.items():
        for criterion, threshold in CRITERIA.items():
            status = 'NOT EVALUABLE'
            detail = evidence
            if criterion.startswith(('D', 'E')):
                detail = 'Not run: prerequisites/matched real spatial inputs unavailable or structurally inadequate. No held-out gain, confidence interval or stability is claimed.'
            if criterion == 'F' and candidate in ('PLAN', 'SPSEQ'):
                status = 'FAIL'
                detail = ('16 animals total but 2 per prespecified timepoint <3; do not pool stages after inspection.'
                          if candidate == 'PLAN' else '3 deposited /4 reported spatial animals <8. Separate scRNAseq animals are not extra spatial test units.')
            if criterion == 'F' and candidate == 'AGE500':
                status = 'PASS'
                detail = 'Replication design only: 12 mice >=8, 6/age group >=3; 3/sex-age cell. Paper and downloaded sample-list code agree. This does not pass A-E or establish usable file-level IDs.'
            if criterion == 'F' and candidate == 'WORM':
                status = 'PASS'
                detail = 'Reported experimental design only: 113 WT and18 unc-31 animals >=8 and>=3/group; not all per-animal payloads were inspected. No ageing/regeneration contrast is present.'
            if criterion == 'A1' and candidate == 'ZFHEART':
                status = 'FAIL'
                detail = 'Measured about35um spatial bins are not individually segmented cell expression observations. Imputed cell expression/locations cannot satisfy the specified measurement requirement.'
            row = next(r for r in catalogue if r['id'] == candidate)
            ledger.append({'candidate': candidate, 'criterion': criterion, 'threshold': threshold,
                           'status': status, 'scope': 'feasibility audit; no biological predictive validation',
                           'evidence': detail, 'source': row['source'], 'data_source': row['data']})
    write_csv(out / 'acceptance_ledger.csv', ledger)
    assert not any(r['status'] == 'PASS' and r['criterion'].startswith(('D', 'E')) for r in ledger)
    print(f'Wrote {len(catalogue)} ranked candidates and {len(ledger)} criterion assessments.')


if __name__ == '__main__':
    main(sys.argv[1])
