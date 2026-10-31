Reproducible Capeling/CosMx response audit

Run from /home/Jack/Real_Data/lymph_node_sheaf. All source data stay read-only in
/home/Jack/Real_Data/spatial_ibd_sheaf. Use a NEW output directory for each run.
Existing outputs are not overwritten. No network access or downloads occur.

python3 -m unittest discover -s tests -p 'test_organoid_response_audit.py' -v
python3 src/organoid_audit/audit_organoid_response_genes.py --out outputs/organoid_response_audit_new --no-download
python3 src/organoid_audit/validate_audit.py outputs/organoid_response_audit_new --workbook
python3 src/organoid_audit/summarize_audit.py outputs/organoid_response_audit_new

Use --root, --previous, --supplement and --catalog to select alternate local inputs.
The catalogue is fingerprinted for provenance; conservative receptor candidates
use the manually reviewed canonical receptor list, not every catalogue partner.
The previous 79-factor CSV provides the factor inventory and catalogue flags.
No new catalogue-wide biological curation is claimed.

Dependencies used for this run are pinned in requirements.txt; all were already
installed. The main XLSX reader uses only the standard library. Independent
validation uses openpyxl. specimens.py preserves the existing four-sample FOV
mapping. Summary includes source and code SHA256 hashes and CSV row manifests.

Important changes relative to the sibling project's original script:
* Header parsing was already correct. Original CSVs and ZIP members contain data.
* Preserve numeric scientific-notation text and source gene names such as THRA1/BTR.
* Exclude 20 Negative and 324 SystemControl targets from panel size and cell QC.
* Detect identical worksheet values and quarantine BOTH ambiguous sheets.
* Preserve signed gene-level response evidence and source identity.
* Separate catalogue binding-partner coverage from canonical signalling complexes.
* Count receptor subunits together within each cell, never across distinct cells.
* Count target detection within each specimen and receptor-positive cells.
* Distinguish covered targets from observed detection and from response activation.
* Strict finite padj <0.05, abs(log2FC)>1; explicit unknown/control/missing rows.
* Fixed TNFALPHA alias; refuse mismatched panels and overwrite of completed runs.
* Round-trip each exported CSV, then independently verify tables and source data.

The 14-factor list is a permissive feasibility screen. Detection (>0 counts) is
not background-corrected activity. Tissue receiver cell identities have not been
assigned; organoid signatures must not yet be called epithelial tissue responses.
Two identical source worksheets cannot support independent stem/goblet claims.
Signed DEG signatures are not the published cNMF programmes (other supplements).
Use receiver-specific signed programmes for modelling; pooled gene unions are
coverage summaries. There is no spatial-neighbour or sheaf calculation here.

Original files were preserved because the sibling directory is outside this
session's writable roots. No commit or push was performed. Do not add the input
matrices, source workbook, or any cell-level observations to version control.
