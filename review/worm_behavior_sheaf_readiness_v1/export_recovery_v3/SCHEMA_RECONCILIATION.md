# Native identity recovery

The original exporter assigned counters in top-level loops and failed Julia soft scope. The corrected exporter wraps execution in `export_identities` and `main`; no global-counter workaround. The failed source is preserved separately. Existing Julia 1.9.2/JLD2 0.6.7 were reused without installation.

Native JLD2 decoding cross-checks every row against the authors' class-indexed `dict_match_dict` view, including label, class, confidence, laterality, dorsoventral identity and annotation field. Across 27 keys, 2,595 entries agree. Lookup keys need not equal the coarser `neuron_class`: 1,692 exact class keys, 181 individual-label keys, 722 explicit D/V refinements. Unknown relationships abort. Raw keys and classes are retained separately.

Observed confidence scores are 2:183, 3:597, 4:515, 5:1,299 and 6:1. The accepted declared range remains 1–5, primary threshold 4. The confidence-6 record is quarantined, never clamped. Three duplicated individual-label groups contain six records; both members are excluded, without choosing the higher score. Seven records are quarantined in total; 2,588 pass structural checks.

The misleadingly named `region` field contains annotation/QC comments, including merged/contaminated ROIs, rather than reliable anatomical labels. Every nonempty comment is conservatively held from primary analysis pending review blind to outcomes. Of 1,813 structurally valid confidence≥4 records, 120 have comments; 118 explicitly concern merges/contamination. There are 1,693 primary clear-note records. No comment text or identity table is published.

The native receipt and TSV are private, source/output SHA-256 linked. The successful destination did not exist before export; no partial outputs were overwritten. Regression tests use synthetic native JLD2 fixtures and fresh Julia subprocesses to reproduce script scope, reject existing outputs, check finer keys, duplicate quarantine and view disagreement. Twenty Julia assertions pass.

All dictionary-wide coverage is metadata only, not evidence that all 27 recordings have usable measured traces. Coarse-class coverage does not imply identical side-resolved neurons. Do not use this audit to select favourable animals or replace the frozen circuit.
