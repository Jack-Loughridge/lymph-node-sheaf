# Completed cleanup within executor permissions

Removed31 exact files:27 Dryad HDF5 recordings, the source JLD2, their duplicate ZIP, and two Crohn count memmaps. Logical bytes981,031,258; measured free-space increase980,271,104bytes. Free space afterward17,793,581,056bytes. The difference reflects allocation/accounting and concurrent filesystem activity; logical bytes are not substituted for measured reclaimed space. An exact private manifest stores each path, size, SHA256, source URL, inode/device, process checks and deletion verification. It also fingerprints767 preserved source/report/specification/forecast-result files, all unchanged.

One preflight stopped before deletion because lsof warned about unrelated hugetlbfs/mqueue mounts. That failed attempt is preserved privately. The successful preflight explicitly checked visible same-user file descriptors and memory maps and whitelisted only those exact virtual-filesystem warnings. This does not claim privileged host-wide process visibility.

The8.98GB Crohn H5AD was already absent and is not counted as reclaimed. Small linkage files, all private metadata, dictionaries/shape specifications, references, all Julia tools, NWB recordings and55MB of frozen forecasting outputs remain. Original microbiome data/databases/models/graphs were never deletion candidates. No recursive directory removal, permission change or symlink traversal occurred.

The6,464,275,219byte IBD GEO archive remains READ_ONLY here, along with extracted CosMx files. No escalation was attempted. Its previously recorded checksum and recoverable GEO source are in the guarded script. To remove that exact archive from a normal SSH session, with checksum and process guards:

```powershell
ssh tda-vm "cd /home/Jack/Real_Data/lymph_node_sheaf && python3 review/perturbational_bioelectric_integration_2026/cleanup_readonly_archive_NORMAL_SSH.py"
```

This command intentionally refuses altered files, symlinks, open/uncertain process checks, or a read-only directory. It never calls sudo/chmod. Additional extracted IBD matrices and remaining Crohn metadata are retained because they are read-only or their recoverability/private-linkage status needs independent confirmation. The archive cleanup has NOT run and is not included in the reclaimed total.

Dryad's public30-file manifest was reverified before deletion. The DANDI current route could not be verified in this session, so its small recordings were conservatively retained. All frozen conclusions can be reconstructed from source manifests and the retained code/audits plus restored public inputs; rerunning prior ingestion now requires restoring the deliberately removed downloads.
