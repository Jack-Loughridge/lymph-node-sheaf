# Lymph-node sheaf: spatial immunosenescence feasibility study

Research code only. **Do not commit human donor-level data, downloads, or authentication tokens.**

## Dataset

Farzad et al. (2026), *A spatial multi-omics atlas of immunosenescence reveals germinal-center B cell dysfunction in human lymph nodes*.
- [Publication](https://doi.org/10.1016/j.cpblue.2026.100053)
- SenNet collection: `SNT679.VSPM.947`
- Spatial DBiT RNA-seq deposit: `SNT248.HNBH.289`
- Spatial ATAC deposit: `SNT347.QGGP.668`

## Run the first audit on the VM

```bash
cd ~/Real_Data/lymph_node_sheaf
python3 src/sennet_inventory.py --out metadata
```

The script makes read-only metadata requests, preferring the normal SenNet Entity API and falling back to the SenNet public Search API on HTTP 401/403/404. It **does not download sequencing or image data**. It writes a metadata JSON and CSV inventory. Review `metadata/status.json` and `metadata/files.csv`.

A 403 on the Entity API does not establish whether the deposited files themselves are protected. The file index's access levels and the portal determine which files are eligible for public download.

If the public search route also fails, test without sharing credentials:

```bash
curl -sS -o /dev/null -w 'Entity: %{http_code}\n' https://entity.api.sennetconsortium.org/entities/SNT248.HNBH.289
curl -sS -o /dev/null -w 'Search: %{http_code}\n' https://search.api.sennetconsortium.org/indices
```

SenNet documentation: https://docs.sennetconsortium.org/apis/getting-started

**VM space:** 26 GB available at project start; no bulk downloads until sizes and access rights have been verified.
