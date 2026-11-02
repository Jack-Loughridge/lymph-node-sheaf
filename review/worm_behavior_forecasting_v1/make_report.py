"""Consolidated report from all prespecified outputs; no pair selection or refitting."""
from pathlib import Path
import json,shutil,hashlib
import pandas as pd,numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]

def main():
    c=pd.read_csv(HERE/'DIRECTED_CONTRASTS.csv');p=pd.read_csv(HERE/'POOLED_CONTRASTS.csv');m=pd.read_csv(HERE/'MODEL_COMPARISONS.csv')
    primary=c[(c.history==10)&(c.grouping=='animal')&(c.scale=='normalised')]
    unconditional=primary[primary.comparison=='VAR_vs_AR'].set_index(['sender','receiver'])
    conditional=primary[primary.comparison=='VAR_beh_vs_AR_beh'].set_index(['sender','receiver'])
    table=[]
    for (sender,receiver),row in conditional.iterrows():
        un=unconditional.loc[(sender,receiver)]
        table.append({'sender':sender,'receiver':receiver,'independent_training_animals':int(row.n_animals),'unconditioned_gain_percent':un.relative_MSE_reduction_percent,'unconditioned_CI_low':un.CI95_low_percent,'unconditioned_CI_high':un.CI95_high_percent,'conditional_gain_percent':row.relative_MSE_reduction_percent,'conditional_CI_low':row.CI95_low_percent,'conditional_CI_high':row.CI95_high_percent,'conditional_animals_improved':int(row.animals_improved),'conditional_leave_one_animal_influence_min':row.evaluation_leave_one_animal_out_min_percent,'conditional_leave_one_animal_influence_max':row.evaluation_leave_one_animal_out_max_percent})
    pd.DataFrame(table).to_csv(HERE/'ALL_DIRECTIONS_PRIMARY.csv',index=False)
    rows='\n'.join(f"|{r['sender']}→{r['receiver']}|{r['independent_training_animals']}|{r['unconditioned_gain_percent']:.2f} [{r['unconditioned_CI_low']:.2f}, {r['unconditioned_CI_high']:.2f}]|{r['conditional_gain_percent']:.2f} [{r['conditional_CI_low']:.2f}, {r['conditional_CI_high']:.2f}]|{r['conditional_animals_improved']}/{r['independent_training_animals']}|" for r in table)
    md='''# Candidate C: first training-only neural forecasting diagnostic

**Decision: do not advance to a dynamical sheaf on the current evidence.** There is a small reproducible cross-neuronal forecast association before behavioural adjustment, but no reproducible practically substantial added information beyond receiver history and prior behaviour. Keep the dataset as a **PARTIAL GO for processed-signal research**, not as a validated local compatibility model. No full sheaf, multiscale metric, validation-six or final-test-six analysis was performed.

## Design and execution

The original fifteen training animals, fixed twelve-neuron circuit and ten eligible unordered pairs were retained; both directions are reported. Their overlap involves six actual individual labels (ASHL/R, AVAL/R, AIBL/R), with3–6 animals per pair. Fourteen training animals contribute at least one eligible pair; the fifteenth remains allocated to training but contributes no supported pair. All56 unsupported pairs remain NOT EVALUABLE rather than being replaced. There are45 pair-animal outer folds /90 direction-animal scores at the primary lag, **not90 independent animals**. All transforms, coefficients and regularisation selection use the appropriate inner/outer training animals only.

`IMPLEMENTATION_PREREGISTRATION.md` and its timestamp/hash receipts were written before any real-data fit. They resolve previously unspecified numerical conventions without changing the original thresholds, identities, split, alpha grid or lag grid. The primary history is ten preceding samples, with5/15 sensitivity. All methods share the same complete numerical windows and past-only cutoff. Receiver AR, bivariate ridge VAR, receiver+behaviour, neural+behaviour, persistence, joint rank-one cross-lag dynamics and an independently pinned anatomical-edge comparator were evaluated. Behavioural inputs include past signed velocity, speed, reversal state, radial distance and pumping; uncalibrated head-angle signs and future behaviour are excluded. These are inputs to **neural forecasting**, not an avoidance-outcome validation.

Inner leave-one-animal-out selects alpha0.1/1/10/100 by equal-animal normalised MSE. Rank-one alpha is selected jointly for two outputs to preserve its genuine2×H cross-filter constraint; unrestricted models tune each receiver separately. This is an explicitly declared comparator difference. The anatomical comparator uses chemical direction and undirected electrical support from a checksum-pinned White composite, with no correlation-derived edges and no conversion of synapse counts to conductance. Its omission of a cross-history is a reference-model assumption, not proof of no biological connection.

## Primary pooled results

Below, improvement is100×(1−augmented MSE/baseline MSE), using training-scale-normalised errors. Losses are averaged over available directions within each animal, then equally across animals. The14-animal pooled95% intervals use10,000 paired animal bootstrap draws; every animal's panels/windows travel together.

| Comparison | Reduction in outer MSE | 95% animal interval | Animals improved |
|---|---:|---:|---:|
|Sender history added to receiver AR|1.19%|0.48 to2.40%|11/14|
|**Sender history added beyond receiver+behaviour**|**0.19%**|**−0.53 to1.06%**|**8/14**|
|Prior behaviour added to receiver AR|8.25%|3.99 to13.60%|12/14|
|Rank-one cross-lag model vs AR|1.12%|0.21 to2.34%|11/14|
|Anatomical-edge VAR vs AR|0.71%|−0.10 to1.78%|6/14|
|AR vs persistence|−1.29%|−9.24 to3.85%|7/14|
|Unconditioned VAR vs persistence|−0.09%|−7.66 to4.87%|9/14|
|Neural+behaviour vs persistence|7.24%|3.92 to10.19%|11/14|

Persistence is therefore a strong comparator, and much of the extra usable signal is already represented by previous behaviour. Behavioural predictiveness could reflect neural state, shared drivers or movement-related processing artefacts; these data cannot separate those explanations. The low-rank and anatomical models do not supply independent evidence of conditional information beyond the behavioural baseline.

The small unconditioned gain is not solely a single-animal result: deleting one evaluated animal at a time leaves pooled gains0.97–1.48%. Conditional gains range−0.034 to0.394% under the same influence check and can cross zero. These are evaluation influence calculations, not nested refits with each animal removed from every training fold.

## Every eligible direction, including negative findings

Values are percentage MSE reductions with descriptive, unadjusted animal-bootstrap95% intervals. “Conditional” compares neural+behaviour against receiver+behaviour. No row is filtered for significance or effect direction. The [aggregate forest plot](FORECAST_DIAGNOSTIC.png) and [PDF version](FORECAST_DIAGNOSTIC.pdf) show the complete directional and lag results.

|Sender→receiver|Animals|Before behaviour: gain [interval]|After behaviour: gain [interval]|Conditional improved animals|
|---|---:|---:|---:|---:|
'''+rows+'''

Only8/20 conditional point estimates are positive;9/20 improve a majority of their animals. Two conditional descriptive intervals exclude zero positively, but neither point estimate reaches5%. The largest conditional point improvement (AVAR→AVAL,6.70%) has an interval spanning−5.85 to25.13%. AIBR→AIBL is3.78% [0.96,6.12] in six animals; ASHR→AIBR is1.20% [0.57,1.76] in only three. These examples are contextualised by all20 rows and are not selected for a follow-up panel. The poor AVAL→AVAR result remains negative when any single evaluated animal is removed. There is no basis for discarding that direction or its recordings.

No conditional direction meets both a≥5% point improvement and a positive lower interval in this diagnostic. This does **not** constitute a confirmatory test of the untouched validation/final-test gate. Direction intervals are unadjusted across20 comparisons, contain only3–6 independent animals, condition on the fitted cross-validation predictions, and do not account for all training-set/refitting variability. Overlapping cross-validation training sets induce dependence. Treat individual interval exclusions as descriptive, not formal causal or multiplicity-controlled discoveries.

## Lag and acquisition-day sensitivities

|History samples|Unconditioned gain [animal95% interval]|Conditional gain [animal95% interval]|
|---|---:|---:|
|5|1.29% [0.60,2.62]|−0.27% [−1.11,0.28]|
|**10 primary**|**1.19% [0.48,2.40]**|**0.19% [−0.53,1.06]**|
|15|1.01% [0.29,2.14]|0.03% [−0.64,0.77]|

The primary lag was not changed. Conditional added value remains weak across all three frozen histories.

Nested acquisition-date exclusion is evaluable for nine pairs. ASHL–AIBL has only two dates and is NOT EVALUABLE for nested date folds, despite three animal recordings. It is not quietly treated as three independent dates or removed from the primary analysis. Day-based results involve13 animals/11 dates. Primary day-held-out unconditioned gain is1.59% [0.92,3.07] and conditional gain0.52% [−0.05,1.48] by animal bootstrap. Date-cluster intervals are0.91–3.07% and−0.03–1.54%, respectively. The matched-support animal-fold analysis gives1.61% unconditioned and0.52% conditional: most of the apparent change from the full animal analysis reflects the unevaluable pair/support difference, not the day-exclusion procedure. `DAY_MATCHED_POOL.csv` reports this explicitly. Date-held-out conditional gains remain near zero at5 and15 samples as well.

## Training performance, uncertainty and numerical audit

`MODEL_COMPARISONS.csv` includes original-fluorescence MSE, training-scale-normalised MSE, apparent training error, chosen alpha counts and outer support for all models/directions/lags/groupings. `DIRECTED_CONTRASTS.csv` includes raw as well as normalised paired gains, bootstrap intervals and influence ranges. Absolute raw fluorescence units are uncalibrated and should not be compared as physiological amplitude across neurons. Normalisation never uses held-out animal statistics.

Apparent training error is predictably optimistic: mean outer-fit training normalised MSE is0.06178 for AR,0.05940 for VAR,0.05441 for receiver+behaviour and0.05369 for neural+behaviour. Equal-animal outer losses are0.07884,0.07791,0.07234 and0.07220, respectively. The apparent numbers average outer fits/panels and are descriptive; they are not independent training replications. Sender information reduces apparent error more readily than animal-held-out error. The alpha grid is not uniformly pinned to one boundary, and all rank-one candidate fits converge under the preregistered rule. No solver failure explains the principal negative finding.

There were14,430 cached candidate fits across primary/sensitivity analyses, with zero rank-one nonconvergences and no rank-one objective increases. Independent prediction audit verifies261 private prediction-file hashes and3,654 directed model MSE calculations. Thirteen synthetic tests verify training-only access guards, inner/outer exclusion, cutoff, gap masks, equal-animal scaling, actual rank-one constraint, graph omission equivalence, moment/direct MSE agreement and clustered uncertainty. All original six protocol hashes and the pre-fit code hashes remain unchanged. The saved run receipt lists only frozen training recording sources. The12 validation/test recordings were not opened by this diagnostic.

## Interpretation and next-stage decision

**NO-GO for progressing directly to a multivariate dynamical sheaf on this evidence.** The data contain some cross-neuronal predictive structure, but its average incremental value after behavioural adjustment is small, unstable in sign and well below the intended practical target in this pairwise diagnostic. The stronger behavioural baseline materially changes the interpretation. Do not select the few favourable directions, extend the alpha/lag grids, unlock validation/test animals or calculate compatibility metrics to rescue the result.

Several limitations remain distinguishable. Three-to-six co-observed animals give weak per-edge precision: each outer fit uses only2–5 animals, and the smallest inner fits use one animal, despite many temporally correlated samples. This is computational eligibility rather than adequate biological power. Moreover, the native labels differ across animals and fully observed larger panels are scarce. Accepted side-resolved identities were joined exactly and no identity inconsistency was found in evaluated panels; this audit cannot certify error-free NeuroPAL identification. Sender information appears substantially redundant with observed behavioural state. Collinearity, unmeasured common input, heterogeneous effective dynamics and observation filtering can also suppress or distort coefficients. Numerical convergence and all-window integrity argue against an implementation failure; a nearly null conditional forecast gain cannot demonstrate biological absence of interaction.

The observations are **supplied processed fluorescence**. Original observation masks and exact provider preprocessing remain absent; interpolation or whole-record corrections may incorporate future information. This is retrospective prediction of processed series, not certified genuine optical observations, causal online prediction, direct synaptic action, gap-junction conductance or measured voltage. Neither a listed anatomical edge nor a predictive improvement resolves mechanism.

Weak pairwise effects do not rule out collective/synergistic neural information. If that distinct question is pursued, it needs a **newly preregistered training-only** missing-state-aware multivariate/state-space experiment using the same fixed circuit and behavioural conditioning. It must use explicit native-identity observation operators, preserve bilateral identities, marginalise missing states without inventing observations, and compare joint-sender information against strong autoregressive, behaviour, low-rank and graph-constrained baselines. It must test joint-versus-individual sender ablations with animal/date grouping and report identifiability/coverage. No such experiment is automatically authorised or fitted here; the pairwise results do not currently justify a sheaf-specific claim. A reliable pathway to biological validation still requires independently measured perturbations or better observation provenance.
'''
    (HERE/'REPORT.md').write_text(md)
    eligibility=pd.read_csv(ROOT/'review/worm_behavior_sheaf_readiness_v1/full_cohort_qc_v5/TRAINING_PAIR_PANEL_COVERAGE.csv');eligibility['status']=np.where(eligibility.nested_animal_folds_possible,'EVALUATED','NOT EVALUABLE: fewer than three training animals');eligibility.to_csv(HERE/'ELIGIBILITY_ALL_66.csv',index=False)
    decision={'decision':'NO-GO for direct progression to a dynamical sheaf on this evidence','dataset_scope':'PARTIAL GO for processed-signal research','primary_unconditioned_gain_percent':float(p[(p.history==10)&(p.grouping=='animal')&(p.comparison=='VAR_vs_AR')].relative_MSE_reduction_percent.iloc[0]),'primary_conditional_gain_percent':float(p[(p.history==10)&(p.grouping=='animal')&(p.comparison=='VAR_beh_vs_AR_beh')].relative_MSE_reduction_percent.iloc[0]),'primary_conditional_CI95_percent':[float(p[(p.history==10)&(p.grouping=='animal')&(p.comparison=='VAR_beh_vs_AR_beh')][k].iloc[0]) for k in ['CI95_low_percent','CI95_high_percent']],'conditional_directions_meeting_5percent_and_positive_CI':0,'eligible_directed_relationships':20,'evaluated_training_animals':14,'validation_test_opened':False,'higher_order_effects_ruled_out':False,'full_sheaf_fitted':False,'recommendation':'Do not unlock held-out cohorts or fit sheaf; any multivariate synergy question requires a separate prospective missing-state-aware experiment'}
    (HERE/'DECISION.json').write_text(json.dumps(decision,indent=2)+'\n')
    ledger=[['training-only recording access','PASS','All14 opened recordings allocated train; held-out12 unopened'],['frozen protocol and pre-fit implementation','PASS','Six original hashes and pre-fit code hashes unchanged'],['complete eligible pair reporting','PASS','10 pairs20 directions; all56 unsupported pairs reported'],['nested acquisition-day analysis','PASS WITH LIMIT','Nine pairs evaluable; one has only two dates'],['numerical model and prediction integrity','PASS','13 tests;261 prediction hashes;3654 direct MSE checks'],['reproducible unconditioned information','SMALL DESCRIPTIVE EFFECT','Pooled1.19percent; interval0.48to2.40'],['conditional information beyond prior behaviour','NOT SUPPORTED','Pooled0.19percent; interval minus0.53to1.06'],['practical fivepercent conditional target','NOT MET IN DIAGNOSTIC','Zero of20 directions meets point and positive interval'],['confirmatory held-out biological gate','NOT EVALUABLE','Validation6 and test6 untouched'],['genuine observation provenance','NOT EVALUABLE','Original masks and exact processing settings absent'],['sheaf or higher-order added value','NOT EVALUABLE','No sheaf or higher-order model fitted'],['direct progression to dynamical sheaf','NO-GO','Conditional evidence insufficient; no favourable pair selection']]
    pd.DataFrame(ledger,columns=['criterion','status','evidence']).to_csv(HERE/'acceptance_ledger.csv',index=False)
    print('Consolidated report and complete20-direction table written')
if __name__=='__main__':main()
