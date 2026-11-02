# Full-cohort transfer checkpoint

The filesystem visible to this executor contains three HDF5 animal recordings, not27. The receipt supplied at `recordings/remaining_transfer_PRIVATE.json` has five entries matching the original three recordings, NeuroPAL dictionary and Julia loader. All five recorded SHA256 hashes and byte sizes verify. None of the24 expected additional recordings is visible in the project/Real_Data search. This is a workspace/transfer discrepancy, not a biological failure or proof that files are absent from another VM filesystem/session.

One existing `remaining_acquisition_PRIVATE.json` journal is present and preserved byte-for-byte. The preflight does not use that journal as a gate. No acquisition retries, deletions, overwrites of previous QC, or model fits were performed. Full-cohort QC cannot start; earlier PARTIAL GO remains unchanged and full-cohort neural/behavioural criteria are NOT EVALUABLE.

`check_transfer.py` reproducibly checks exact receipt membership, hashes, sizes and cohort presence against the frozen public manifests. It exits2 on this blocked checkpoint. Synthetic tests demonstrate that an old valid receipt cannot pass as the new24-file receipt and that corruption is rejected. Individual hashes and animal inventories remain private; public output is aggregate only.

Please inspect the same absolute path through your normal Windows SSH connection and reconcile its view with this executor. If all27 are visible there, provide that command output rather than downloading/transferring the data again. If only three are visible there too, transfer the already-downloaded remaining24 and the correct Windows-generated receipt using the previous manual procedure; no repeat download is needed.

The next model remains the previously frozen short-history neural forecasting diagnostic within training15, with receiver autoregression, regularised multivariate dynamics, pre-target covariates, and graph/low-rank baselines. It cannot yet be run on the missing cohort. No thresholds or split allocations were modified.
