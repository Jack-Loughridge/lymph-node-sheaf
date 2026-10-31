Minimal IFNG / TNF / shared BMP spatial factor-graph sheaf

Mathematics and parameter provenance were documented before implementation:
  MODEL_SPECIFICATION.txt
The implemented eliminated operator is derived from an actual affine cellular
sheaf on the variable-factor incidence graph. It is not a pairwise rule equating
every sender with a receiver's whole response. See algebra.py and its tests.

Run from /home/Jack/Real_Data/lymph_node_sheaf (Python 3.10.12 in this pilot).
Installed numerical dependency versions are pinned in requirements.txt.
No new datasets or packages were downloaded for the pilot.

python3 -m unittest discover -s tests -p 'test_signalling_sheaf.py' -v
python3 src/signalling_sheaf/prepare.py
python3 src/signalling_sheaf/run_pilot.py --out outputs/signalling_sheaf_new --allocator kernel --nulls 20
python3 src/signalling_sheaf/validate.py outputs/signalling_sheaf_new
python3 src/signalling_sheaf/plot_results.py outputs/signalling_sheaf_new
python3 src/signalling_sheaf/diagnose_results.py outputs/signalling_sheaf_new

Use a fresh --out directory; existing completed results are never overwritten.
Inputs default to the existing sibling spatial_ibd_sheaf/pilot_6k project and the
corrected audit's signed_response_genes.csv in this workspace. Override --root,
--signatures, --cache as needed. Existing inputs are read-only. The .npz cache
contains only selected-FOV cell observations; keep it under ignored data/ and do
not commit or include it in a shared review archive. Use a new --cache path when
changing data selection or annotation logic. Sources are SHA256 fingerprinted.

FOVs are fixed before response modelling: GSM9345416 FOV37; GSM9345417 FOV25;
GSM9345418 FOV21; GSM9345419 FOV6. These were previously selected by median metadata
cell count per specimen. Cell joins use (slide prefix, FOV, cell string) and enforce
one-to-one cardinality, numeric cell-ID consistency and no unmatched expression
rows. Numeric IDs alone must never be used across FOVs or slides.

ANNOTATIONS
No published subtype labels are present in either metadata CSV; 'assay_type' is
not a cell-type annotation. Earlier project assignments are marker proxies.
Current TA-like/colonocyte-like assignments remain provisional. Several preferred
absorptive markers (CA1/SLC26A3/AQP8/MS4A12) are absent from the panel. The covered
absorptive recipe includes less specific maturation/progenitor markers. Do not
claim validated organoid-to-tissue subtype matching from these assignments.
PanCK/CD45/CD68 medians are corroborating summaries, not calibrated subtype gates.
Ambiguous epithelial and other cells remain possible senders but receive no map.

RESPONSE FEATURES AND ALLOCATION
The gene manifest retains raw signed log2FC, adjusted p, organoid source and every
exclusion reason. All annotation markers and all selected input ligand/receptor
genes are removed from response features. In particular, TNF transcript cannot
serve simultaneously as an input and as a TNF response target. BMP2/BMP6 input
routes remain distinct, but use one averaged signed response programme. This is
not an attempt to identify their separate effects.

COMMOT and anndata were not installed. The completed run explicitly uses a
bounded spatial kernel allocation (not COMMOT and not an OT solver). The optional
--allocator commot adapter uses the documented custom LR interface, min-subunit
complexes, a spatial/target/self mask and conservative post-allocation caps. It
fails explicitly if the packages are missing; it was not executed or validated
against an installed COMMOT version here. Do not describe these pilot outputs as
COMMOT results. Test the adapter in an isolated compatible environment before
using it for biological conclusions. No silent backend substitution occurs.

Only conservative routes IFNG/IFNGR1_IFNGR2, TNF/TNFRSF1A, BMP2/BMPR1A_BMPR2 and
BMP6/ACVR1_BMPR2 are used. Other known routes are omitted. Both BMP ligands share
the BMPR2 capacity budget. Proxies use log-normalised transcript abundance, not
protein amounts. Radius is 2.5 median cell diameters in PIXELS, not microns.
No self, cross-FOV or cross-specimen links. Receptor subunits must co-occur in the
same receiver. Allocation uses no response data and can receive multiple senders.

OUTPUTS
report.txt: concise findings and limits.
parameter_provenance.json: external experiments versus inferred versus unknown.
annotation_inspection.json / annotation_summary.csv: available columns and proxies.
signed_response_gene_manifest.csv / response_weights.csv: signed targets and roles.
allocation_summary.csv: receiver coverage, multi-sender counts, spatial bounds.
fixed_allocation_residuals.csv: primary discrepancies across gain grid and nulls.
sheaf_energies.csv: secondary sparse least-squares diagnostic with anchored source
and response observations; allocation/aggregation/response energy components.
residual_summary.csv: weighted FOV pooling and descriptive shuffle intervals.
algebra_checks.csv / validation.json: numerical and export checks.
run_config.json: seed, selected data, software versions, code and data hashes,
and manifests for all 13 numerical CSVs.

The cell-intrinsic baseline is a receiver-type intercept (gain=0). Signature gene
and programme standardisation use pooled selected receivers within each proxy
type: this is transductive, not held-out prediction. No disease labels enter maps.
Gains 0,.25,.5,1 are prespecified assumptions; .5 is arbitrary nominal calibration.
No gain is fitted or selected for best performance. Normalisation parameters are
frozen across nulls. Each spatial null moves complete four-ligand bundles within
FOV/proxy-type strata, holding receptors, responses and geometry fixed. There are
20 seeded draws. Quantiles and rank fractions are descriptive, not disease tests.

Use unweighted fixed-allocation response MSE as the primary comparison. Smaller
relaxed sheaf energy can result from flexibility or its degree/gain-dependent
metric. Relaxed auxiliary fluxes and exposures can be negative; they are algebraic
variables, not biological allocation estimates. Nonnegative q0 is the allocation.
Unexposed receivers remain in all denominators. FOV boundaries truncate signals.
No classifier, learned restriction maps, dose calibration or ageing claim is made.

The review ZIP includes code/tests/aggregate results only. Original data are
preserved, and nothing has been committed or pushed.

residual_decomposition.json separates response/exposure cross moments from exposure
amplitude in the gain-grid MSE. It does not estimate or select a response gain.
