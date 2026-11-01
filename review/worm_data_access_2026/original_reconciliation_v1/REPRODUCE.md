# Reproduce without publishing observations

The commit contains code, aggregate counts, source hashes and public access metadata. It deliberately excludes the 29-row event ledger, NWBs, OSF archives, author source caches and individual traces.

Existing VM, from repository root:

```bash
bash review/worm_data_access_2026/original_reconciliation_v1/run_all.sh
```

Override the interpreter with `WORM_QC_PYTHON=/path/to/python` if necessary. The existing VM interpreter is `/home/Jack/micromamba/envs/humann_py312/bin/python` (Python 3.12, NumPy 2.5.0, h5py 3.16.0). No packages are installed or updated by the runner. The supplied prior-QC helper code and rules are unchanged dependency copies; their main functions are not executed by this reconciliation.

Fresh checkout:

1. Place the three previously selected, hash-matching NWBs at ignored `data/worm_nwb/sub-36.nwb`, `sub-72.nwb`, `sub-80.nwb`.
2. Restore only public metadata/reference text using `python3 review/worm_data_access_2026/code/fetch_metadata.py` on a machine with permitted HTTPS and Python requests; it checks the supplied source hashes. These are small text resources, not extra recordings.
3. Run the reconciliation runner above with an interpreter having NumPy and h5py. Pinned author sources are restored into ignored data/ by `fetch_sources.py`; failures stop without changing network settings.
4. For originals, follow `OSF_MANUAL_TRANSFER.md`, then pass `--wt-archive data/worm_original_exports/exported_data.tar.gz` to the runner. The version-1 MD5 is checked before parsing that known filename. Different versions must be investigated explicitly.

The local detailed result is `data/worm_data_access_2026/original_reconciliation_v1/event_failure_decomposition.csv`. It is freshly regenerated from measured NWB arrays; no prior animal/event-level checkpoint is required. Original export statuses remain UNKNOWN until actual files pass session, dimension and event-alignment checks. The scanner never executes pickle or extracts tar paths to disk. An archive with only dense processed matrices does not establish that missing measured values were recovered.

Frozen primary QC remains unchanged; component checks and the optional original-data crosswalk have separate statuses. Running with new original files may change the aggregate reconciliation outputs prospectively, so retain the current checkpoint before updating any public report. Do not publish the resulting individual event ledger.
