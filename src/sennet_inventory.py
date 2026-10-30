#!/usr/bin/env python3
"""Read-only Farzad/SenNet public dataset audit (no bulk file downloads).

Run: python3 src/sennet_inventory.py --out metadata
Uses Python standard library only. Never commit metadata or downloaded human data.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import sys
import time
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

RNA_ID = "SNT248.HNBH.289"
COLLECTION_ID = "SNT679.VSPM.947"
ENTITY = "https://entity.api.sennetconsortium.org"
SEARCH = "https://search.api.sennetconsortium.org"
PROCESSED = (".h5ad", ".h5", ".hdf5", ".mtx", ".mtx.gz", ".csv",
             ".csv.gz", ".tsv", ".tsv.gz", ".parquet", ".rds", ".json")
RAW = (".fastq", ".fastq.gz", ".fq", ".fq.gz", ".bam", ".cram", ".sra")


class ApiError(Exception):
    def __init__(self, url, status, detail=""):
        self.status = status
        super().__init__(f"{url}: HTTP {status} {detail[:120]}")


def get_json(url, payload=None):
    body = json.dumps(payload).encode() if payload is not None else None
    headers = {"Accept": "application/json",
               "User-Agent": "Mozilla/5.0 (compatible; spatial-sheaf-audit/0.3)"}
    if body is not None:
        headers["Content-Type"] = "application/json"
    token = os.environ.get("SENNET_TOKEN", "").strip()
    if token:
        headers["Authorization"] = "Bearer " + token
    for attempt in range(3):
        req = Request(url, data=body, headers=headers)
        try:
            with urlopen(req, timeout=60) as response:
                return json.load(response)
        except HTTPError as exc:
            if exc.code in (401, 403, 404) or attempt == 2:
                raise ApiError(url, exc.code) from exc
        except (URLError, TimeoutError) as exc:
            if attempt == 2:
                raise RuntimeError(f"{url}: {exc}") from exc
        time.sleep(2 ** attempt)


def index_lookup(sid):
    for field in ("sennet_id.keyword", "sennet_id"):
        response = get_json(SEARCH + "/entities/search",
                            {"size": 2, "query": {"term": {field: sid}}})
        matches = response.get("hits", {}).get("hits", [])
        if matches:
            return matches[0].get("_source", {})
    raise LookupError(f"{sid} not found in publicly indexed entities")


def get_entity(sid):
    try:
        return get_json(ENTITY + "/entities/" + quote(sid, safe="")), "entity_api"
    except ApiError as exc:
        if exc.status not in (401, 403, 404):
            raise
        print(f"Entity API HTTP {exc.status} for {sid}; trying public Search API.",
              file=sys.stderr)
        return index_lookup(sid), "public_search_fallback"


def member_ids(collection):
    values = collection.get("entities") or collection.get("datasets") or []
    ids = []
    for v in values:
        sid = (v.get("sennet_id") or v.get("sennetId")) if isinstance(v, dict) else v
        if isinstance(sid, str) and sid.startswith("SNT") and sid not in ids:
            ids.append(sid)
    return ids


def list_indexed_files(uuid):
    found = []
    offset = 0
    while offset < 10000:
        query = {
            "query": {"term": {"dataset_uuid.keyword": uuid}},
            "size": 300,
            "from": offset,
            "_source": ["dataset_uuid", "rel_path", "size",
                        "file_extension", "data_access_level",
                        "sha256_checksum", "md5_checksum", "file_uuid"]
        }
        response = get_json(SEARCH + "/files/search", query)
        hits = response.get("hits", {}).get("hits", [])
        found.extend(h.get("_source", {}) for h in hits)
        if len(hits) < 300:
            return found
        offset += len(hits)
    print(f"WARNING: reached 10,000-file pagination limit for {uuid}", file=sys.stderr)
    return found


def csv_write(path, records, keys):
    all_keys = list(dict.fromkeys(keys + [k for row in records for k in row]))
    with open(path, "w", newline="", encoding="utf8") as f:
        writer = csv.DictWriter(f, fieldnames=all_keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def likely_processed(path):
    name = (path or "").lower()
    if name.endswith(RAW):
        return False
    if not name.endswith(PROCESSED):
        return False
    return (any(t in name for t in ("processed", "matrix", "count", "expression",
                                     "spatial", "metadata", "filtered", "anndata"))
            or name.endswith((".h5ad", ".rds", ".mtx", ".mtx.gz")))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="metadata")
    args = parser.parse_args()
    dest = Path(args.out)
    dest.mkdir(parents=True, exist_ok=True)
    status = {"collection_id": COLLECTION_ID, "rna_id": RNA_ID,
              "note": "Inventory only; no files downloaded"}
    rows, files = [], []
    try:
        collection, mode = get_entity(COLLECTION_ID)
        status["collection_lookup"] = mode
        (dest / "collection.json").write_text(json.dumps(collection, indent=2))
        ids = member_ids(collection)
    except Exception as ex:
        status["collection_error"] = str(ex)
        print(f"Collection lookup failed ({ex}); auditing the RNA upload alone.",
              file=sys.stderr)
        ids = []
    if RNA_ID not in ids:
        ids.append(RNA_ID)
    for sid in ids:
        try:
            entity, mode = get_entity(sid)
            uuid = entity.get("uuid")
            access = entity.get("data_access_level", "unknown")
            kind = entity.get("entity_type", "unknown")
            rec = {"sennet_id": sid, "uuid": uuid, "access": access,
                   "entity_type": kind, "dataset_type": entity.get("dataset_type"),
                   "lookup": mode, "title": entity.get("title", "")}
            rows.append(rec)
            print(f"{sid}: {kind} / {rec['dataset_type']} / {access} / via {mode}")
            if sid == RNA_ID:
                status["rna_lookup"] = mode
                status["rna_access"] = access
                status["rna_uuid"] = uuid
                (dest / "rna_upload.json").write_text(json.dumps(entity, indent=2))
            if kind == "Dataset" and uuid:
                for item in list_indexed_files(uuid):
                    item["sennet_id"] = sid
                    item["uuid"] = uuid
                    files.append(item)
        except Exception as ex:
            print(f"WARNING: {sid}: {ex}", file=sys.stderr)
            rows.append({"sennet_id": sid, "error": str(ex)})
            if sid == RNA_ID:
                status["rna_error"] = str(ex)
    csv_write(dest / "datasets.csv", rows,
              ["sennet_id", "uuid", "entity_type", "dataset_type", "access", "lookup", "error"])
    csv_write(dest / "files.csv", files,
              ["sennet_id", "rel_path", "size", "data_access_level", "sha256_checksum"])
    candidates = [f for f in files if f.get("data_access_level") == "public"
                  and likely_processed(f.get("rel_path", ""))]
    csv_write(dest / "candidate_public_processed.csv", candidates,
              ["sennet_id", "rel_path", "size", "data_access_level"])
    status["indexed_files"] = len(files)
    status["candidate_public_processed"] = len(candidates)
    (dest / "status.json").write_text(json.dumps(status, indent=2))
    print(f"Inventory: {len(rows)} entities, {len(files)} indexed files, "
          f"{len(candidates)} candidate public processed files. Output: {dest}/")
    if "rna_error" in status:
        print("RNA deposit metadata still unavailable. This may be a SenNet access "
              "restriction or API issue, not a VM problem.", file=sys.stderr)
        print("Do not paste access tokens into GitHub or ChatGPT.", file=sys.stderr)
        raise SystemExit(2)
    print("No data files were downloaded. Review files.csv before any transfer.")


if __name__ == "__main__":
    main()
