# Prospective validation protocol and stopping rules — no fit performed

The old sheaf specifications, transcriptomic criteria and frozen worm experiments are untouched. This protocol concerns the new A–F gates, with eight independent biological units per primary condition as a planning target. It does not silently recast cells/FOVs as units or turn that target into a universal biological minimum.

## Acquisition gates before analysis

1. Obtain one minimal control and disrupted recording, and rescue if genuinely available. Inspect real array headers and data values, clocks, actual simultaneously observed cells, geometry, stimulus/controller logs and observation/provenance flags. Confirm treatment and animal/culture identity joins. No GO from a figure caption, code schema, source-data label or fifty overlaid traces.
2. Distinguish reporter/clone, culture dish, animal, slice, FOV, repeated neurotransmitter trial and cell. No train/test leakage through a shared animal/dish or sequential measurements. Confirm whether a date uniquely identifies an animal from authoritative mapping rather than assume it.
3. Check matched cells/regions across the functional boundary and independent target versus receiving-population state. Two-cell clamp records plus separately acquired three-cell spontaneous records cannot be combined into an observed three-cell input-output trajectory.
4. Require negative/laser/vehicle controls appropriate to the modality and intrinsic-excitability controls. Prefer genetic knockout/re-expression, but do not label Cx36 overexpression in wild-type worms or Nifedipine compensation as rescue of the disrupted junction itself. Record original assignment, blinding, batch balance and exclusions; do not invent randomisation retrospectively.
5. Abort the primary model if times, cell roster, simultaneously measured boundary states or biological-unit identities are irrecoverably absent. A request for an existing minimal author export is appropriate; a huge download cannot create an unmeasured third cell or missing experimental control.

## Biological-unit holdout and estimands

Before examining response effect sizes, freeze actual unit/condition IDs, source hashes, QC rules, an anatomical region definition, stimuli and response windows, model family/tuning grid, normalisation/calibration and nulls. Use animal/dish-grouped outer validation with all slices/FOVs/trials inside that unit; block by acquisition batch/day as a separate sensitivity. Learn normalisation, noise, sensor transforms and parameters within training units only. Do not use inactive versus active responses to define regions or choose cells. Fix anatomical regions independently of the response/integration score.

For a sufficiently replicated public cohort, reserve an untouched final biological-unit set before fitting. If only four WT mice or three-to-six dishes are available, the proposed eight-unit planning target is not met and paired/low-variance justification has not been demonstrated. Use any resulting grouped cross-validation as an explicitly small exploratory pilot; do not claim confirmatory precision. Inspect baseline variability/within-unit pairing and conduct a prospectively documented power/precision calculation before collecting or unlocking further data. No invented sample-size exception after observing a favourable effect.

Primary physiological target: held-out response of unilluminated/unclamped cells/regions to a recorded local command, conditional on their own prior state and known common inputs. Distinguish desired optical command, actual actuator exposure and achieved targeted-cell response. Calibrate sensor gain or retain a labelled optical-state estimand. Use repeated trials within unit for precision, not as new independent units. Counterbalance/interleave sham and stimulation if available; absent controls cannot be manufactured by circular shifts.

Independently compare control-trained regional internal constraints and boundary constraints in control versus disrupted/rescue units. Recalibration in disrupted training data is a separate predeclared analysis. Do not tune maps by minimising the evaluation region's sheaf energy and then call that predictive success. Require a demonstrated mechanism-relevant perturbation effect and a useful held-out response predictor before evaluating any sheaf-derived organisation metric. If propagation is absent in the sampled assay, stop rather than broaden response windows or select successful FOVs.

## Baselines and nulls

- Receiver-only autoregressive and known-input/context models; match stimulus histories and observed windows.
- Regularised VAR or latent state-space dynamics with explicit missing-state observation operators; identical training budget and noise/calibration information.
- Low-rank/graph-constrained dynamics, dynamic correlations/standard functional connectivity, graph connectivity/modularity/Laplacians, and ordinary stimulation-to-response transfer maps.
- Spatial/contact-distance-matched rewiring or pairing nulls preserving cell types, internal regions, missingness and per-unit activity statistics; no biologically impossible long-range edges introduced just to weaken a null.
- Intervention-label permutations at the **independent assignment unit**, respecting paired/pre-post and batch structure. If genotype is perfectly confounded with batch, label shuffling cannot resolve it; record the comparison as confounded. Within-cell time shifts are secondary dynamics nulls and must preserve autocorrelation and observation masks, avoid the stimulus window and never replace sham/excitability controls.

Use identical targets and numerical/observational masks in all comparisons. Bootstrap or randomise independent animals/cultures with all their fields and trials together. Report between-unit uncertainty, sign consistency, raw/normalised effects, missing coverage and all prespecified conditions, including unsuccessful ones. Distinguish uncertainty conditional on fitted maps from uncertainty refitting them. Report multiplicity and do not treat many cell pairs as independent biological replication.

## Go/no-go transitions

**DATA GO** requires real inspected, usable local-state arrays, geometry/identity, timed perturbations, independent units and valid controls. **LOCAL BIOLOGICAL GO** additionally requires a replicated mechanism-dependent propagation effect in our own audit with off-target/excitability alternatives assessed. **MULTISCALE SHEAF GO** further requires a simultaneously observed multiport boundary and independent tests of gluing beyond conventional models. **GOAL-DIRECTEDNESS EXTENSION** requires a separately measured, meaningfully linked outcome; group-level phenotypes cannot be promoted into same-animal prediction targets.

No numerical evidence for our maps or integration hypothesis has been generated in this search. A local physiological benchmark may be valuable even if no multiscale or goal-directedness claim can follow. If only pairwise recordings exist, retain the allowed smaller-system physiological benchmark and single-boundary compatibility test; stop claims of distinct multicellular multiport/hierarchical integration unless the actual observation/validation structure supports them. Do not impose twenty cells, a third stimulated-recorded cell, or absolute voltage as universal lower bounds. If raw masks are absent but processed signals are interpretable, label a narrower prospective optical-state study and explicitly retain the provenance gate as unresolved; no fabricated observations.

No further unrestricted search, final sheaf fit, H0, long-distance integration, regional section-distance calculation or pushforward construction is authorised by this document.
