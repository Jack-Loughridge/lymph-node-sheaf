"""Reacquire the ten immutable small Git blobs when executor access is available.

No large matrices or human cell observations are requested. Existing exact files
are reused; corrupt files are rejected, never overwritten. Writes tool-compatible
envelopes under ignored data/ and a separate retrieval log. DNS failure is logged
as failure of this executor, not evidence that the public dataset does not exist.
"""
import argparse
import base64
import datetime
import json
from pathlib import Path
import urllib.request

from audit_small_files import verify_blob


def main(cache, log):
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    records = json.loads(Path(
        'review/dataset_search_2026/small_file_audit/SMALL_FILE_AUDIT.json'
    ).read_text())['inputs']
    attempts = []
    for record in records:
        destination = cache / (record['name'] + '.json')
        outcome = {'name': record['name'], 'url': record['source_url'],
                   'utc': datetime.datetime.now(datetime.timezone.utc).isoformat()}
        try:
            if destination.exists():
                envelope = json.loads(destination.read_text())
                content = envelope['structuredContent']['content']
                outcome['method'] = 'verified existing cache'
            else:
                request = urllib.request.Request(record['source_url'], headers={
                    'Accept': 'application/vnd.github+json',
                    'User-Agent': 'cellular-sheaf-small-source-audit'})
                with urllib.request.urlopen(request, timeout=10) as response:
                    body = response.read(1_000_001)
                if len(body) > 1_000_000:
                    raise ValueError('small source response exceeds 1 MB limit')
                obj = json.loads(body)
                if obj.get('encoding') != 'base64':
                    raise ValueError('not a Git blob JSON response')
                content = base64.b64decode(obj['content']).decode('utf-8')
                envelope = {'structuredContent': {
                    'content': content, 'url': record['source_url']}}
                outcome['method'] = 'GitHub immutable blob'
            verify_blob(content, record['git_blob'])
            if not destination.exists():
                with destination.open('x') as handle:
                    json.dump(envelope, handle)
            outcome['status'] = 'retrieved and verified'
        except Exception as exc:
            outcome.update(status='not retrieved/verified', error=str(exc))
        attempts.append(outcome)
    with Path(log).open('x') as handle:
        json.dump(attempts, handle, indent=2)
    print(f'{sum(r["status"] == "retrieved and verified" for r in attempts)}/{len(attempts)} verified')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cache', default='data/dataset_search_2026')
    parser.add_argument('--log', required=True)
    args = parser.parse_args()
    main(args.cache, args.log)
