"""Native TSV -> actual HDF5 rows. Never decode Julia object graphs in Python.
All individual/animal tables private. Aggregate reports do not contain them.
"""
from pathlib import Path
import json,csv,tomllib,hashlib,sys
import pandas as pd,numpy as np,h5py
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
PRIVATE=ROOT/'data/worm_behavior_sheaf_readiness_v1/export_recovery_v3'
NATIVE=ROOT/'data/worm_behavior_sheaf_readiness_v1/sample_ingestion_v2/native_identity'
CFG=json.loads((HERE.parent/'READINESS_PROTOCOL.json').read_text())

def common_core_witness(coverage,min_classes=15,min_animals=15):
    c=np.asarray(coverage,bool);masks=[sum(1<<int(i) for i in np.flatnonzero(c[:,j])) for j in range(c.shape[1])]
    if c.shape[0]<min_classes or c.shape[1]<min_animals:return None
    def visit(start,chosen,mask):
        if mask.bit_count()<min_classes or len(masks)-start<min_animals-len(chosen):return None
        if len(chosen)==min_animals:return chosen,[i for i in range(c.shape[0]) if mask>>i&1]
        for j in range(start,len(masks)-(min_animals-len(chosen))+1):
            out=visit(j+1,chosen+[j],mask&masks[j])
            if out is not None:return out
        return None
    return visit(0,[],(1<<c.shape[0])-1)

def primary_identities(table,minimum_confidence=4):
    """Reject native-invalid rows and hold every nonempty free-text annotation.
    Actual region values contain merge/contamination notes, not reliable anatomy.
    """
    valid=table.identity_valid.astype(str).str.lower()=='true'
    return table[valid & (table.confidence>=minimum_confidence) & (table.confidence<=5) & (table.region.str.strip()=='')].copy()

def main():
    receipt=tomllib.loads((NATIVE/'native_receipt_PRIVATE.toml').read_text())
    source=ROOT/'data/worm_behavior_sheaf_readiness_v1/recordings/NeuroPAL_labels.jld2'
    if hashlib.sha256(source.read_bytes()).hexdigest()!=receipt['source_sha256']:raise ValueError('Native source hash mismatch')
    table=NATIVE/'labels_PRIVATE.tsv'
    if hashlib.sha256(table.read_bytes()).hexdigest()!=receipt['labels_sha256']:raise ValueError('Native output hash mismatch')
    d=pd.read_csv(table,sep='\t',keep_default_na=False)
    valid=d[d.identity_valid.astype(str).str.lower()=='true'].copy();confidence_only=valid[valid.confidence>=CFG['identity_primary_confidence']].copy();high=primary_identities(d,CFG['identity_primary_confidence'])
    ids=sorted(json.loads((ROOT/'data/worm_behavior_sheaf_readiness_v1/animal_splits_PRIVATE.json').read_text()))
    if set(d.animal)!=set(ids):raise ValueError('Native dictionary keys do not exactly match frozen 27-animal manifest')
    cohort=[];samples=[];sample_ids=[];sample_labels={};sample_classes={}
    for animal in ids:
        a=d[d.animal==animal];h=high[high.animal==animal]
        cohort.append(dict(animal=animal,native_high_confidence_before_note_review=len(confidence_only[confidence_only.animal==animal]),annotation_note_holds=int(((confidence_only.animal==animal)&(confidence_only.region.str.strip()!='')).sum()),label_rows=len(a),valid_label_rows=int((a.identity_valid.astype(str).str.lower()=='true').sum()),high_confidence_individuals=len(h),high_confidence_classes=h.neuron_class.nunique(),high_side_resolved=int(h.LR.isin(['L','R']).sum())))
        p=ROOT/'data/worm_behavior_sheaf_readiness_v1/recordings'/f'{animal}-cudata.h5'
        if not p.exists():continue
        with h5py.File(p,'r') as f:
            n=int(f['gcamp/trace_array_original'].shape[1]);t=f['timing/timestamp_confocal'].size
            if f['gcamp/trace_array_original'].shape[0]!=t:raise ValueError('Unexpected Julia-HDF5 storage axes')
            if not a.trace_row_1based.between(1,n).all():raise ValueError('Native identity outside actual recording row bounds')
            if h.label.duplicated().any() or h.trace_row_1based.duplicated().any():raise ValueError('Unquarantined duplicate native identity')
            for _,r in h.iterrows():
                trace=f['gcamp/trace_array_original'][:,int(r.trace_row_1based)-1]
                if not np.isfinite(trace).all() or np.std(trace)==0:raise ValueError('A primary identified trace fails numerical checks')
            sample_ids.append(animal);sample_labels[animal]=set(h.label);sample_classes[animal]=set(h.neuron_class)
            samples.append(dict(animal=animal,total_ROIs=n,total_native_labels=len(a),native_high_confidence_before_note_review=len(confidence_only[confidence_only.animal==animal]),annotation_note_holds=int(((confidence_only.animal==animal)&(confidence_only.region.str.strip()!='')).sum()),confidence_ge4_valid_individuals=len(h),high_identity_fraction=len(h)/n,confidence_ge4_classes=h.neuron_class.nunique(),side_resolved_high_individuals=int(h.LR.isin(['L','R']).sum()),DV_resolved_high_individuals=int(h.DV.isin(['D','V']).sum()),quarantined_rows=int((a.identity_valid.astype(str).str.lower()!='true').sum())))
    (PRIVATE/'sample_identity_QC_PRIVATE.json').write_text(json.dumps(samples,indent=2)+'\n')
    pd.DataFrame(samples).to_csv(PRIVATE/'sample_identity_QC_PRIVATE.csv',index=False)
    pd.DataFrame(cohort).to_csv(PRIVATE/'all27_dictionary_coverage_PRIVATE.csv',index=False)
    for field,name in [('neuron_class','class'),('label','individual')]:
        cov=pd.crosstab(high[field],high.animal).reindex(columns=ids,fill_value=0)>0
        cov.astype(int).to_csv(PRIVATE/f'{name}_by_animal_CONF4_PRIVATE.csv')
    cc=pd.crosstab(high.neuron_class,high.animal).reindex(columns=ids,fill_value=0)>0
    witness=common_core_witness(cc.values)
    ic=pd.crosstab(high.label,high.animal).reindex(columns=ids,fill_value=0)>0
    individual_witness=common_core_witness(ic.values)
    if witness:
        animals,classes=witness;(PRIVATE/'class_gate_metadata_witness_PRIVATE.json').write_text(json.dumps({'animals':[ids[i] for i in animals],'classes':[cc.index[i] for i in classes],'interpretation':'Identity-only upper bound, not a usable measured-cohort gate'},indent=2)+'\n')
    common_individual=set.intersection(*(sample_labels[a] for a in sample_ids));common_classes=set.intersection(*(sample_classes[a] for a in sample_ids))
    common_rows=high[high.label.isin(common_individual)]
    # Pairs require identical class and DV metadata plus explicit L/R; no label stripping.
    pairkeys=[]
    for key,g in common_rows.groupby(['neuron_class','DV']):
        if all({'L','R'}<=set(g[g.animal==a].LR) for a in sample_ids):pairkeys.append(key)
    aggregate=dict(native_dictionary_animals=receipt['animal_dictionary_keys'],native_rows=receipt['rows_exported'],two_view_entries_cross_checked=receipt['entries_cross_checked'],valid_identity_rows=receipt['valid_identity_rows'],quarantined_rows=receipt['quarantined_rows'],duplicate_individual_label_groups=receipt['duplicate_individual_label_groups'],duplicate_rows_quarantined=receipt['duplicate_individual_label_rows'],confidence_counts=receipt['confidence_counts'],observed_confidence_range=[int(d.confidence.min()),int(d.confidence.max())],accepted_confidence_range=[1,5],invalid_confidence_not_clamped=receipt['quarantine_reasons'].get('invalid_confidence',0),class_refinements=receipt['class_lookup_relations'],sample_recordings_inspected=len(samples),sample_primary_individual_count_range=[min(r['confidence_ge4_valid_individuals'] for r in samples),max(r['confidence_ge4_valid_individuals'] for r in samples)],sample_primary_individuals_total=sum(r['confidence_ge4_valid_individuals'] for r in samples),sample_primary_class_count_range=[min(r['confidence_ge4_classes'] for r in samples),max(r['confidence_ge4_classes'] for r in samples)],sample_common_primary_classes=len(common_classes),sample_common_primary_individuals=len(common_individual),sample_common_complete_LR_DV_groups=len(pairkeys),sample_quarantined_rows=sum(r['quarantined_rows'] for r in samples),all27_native_valid_high_confidence_rows=len(confidence_only),all27_annotation_note_holds=len(confidence_only)-len(high),region_field_is_not_reliable_anatomical_annotation=True,all27_primary_unambiguous_high_confidence_rows=len(high),all27_metadata_classes=high.neuron_class.nunique(),all27_metadata_common_classes=len(set.intersection(*(set(high[high.animal==a].neuron_class) for a in ids))),same15_classes_in_same15_animals_METADATA_ONLY=witness is not None,same15_individuals_in_same15_animals_METADATA_ONLY=individual_witness is not None,full_cohort_recording_quality_evaluated=False,julia_version=receipt['julia_version'],JLD2_version=receipt['JLD2_version'],genuine_observation_mask_available=False,full_sheaf_fitted=False)
    (HERE/'IDENTITY_QC_AGGREGATE.json').write_text(json.dumps(aggregate,indent=2)+'\n');print(json.dumps(aggregate,indent=2))
if __name__=='__main__':main()
