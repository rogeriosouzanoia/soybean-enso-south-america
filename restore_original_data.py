"""Download, verify and restore the frozen original project data (standard library)."""
from pathlib import Path, PurePosixPath
import argparse, hashlib, json, shutil, urllib.request, zipfile

ROOT = Path(__file__).resolve().parent
RELEASE = 'https://github.com/rogeriosouzanoia/soybean-enso-south-america/releases/download/v2026.11.10/'

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--download', action='store_true', help='Download the three public release assets when absent')
    parser.add_argument('--archive-dir', type=Path, default=ROOT / '.data_archives')
    parser.add_argument('--target', type=Path, default=ROOT / 'original_project')
    args = parser.parse_args()
    manifest = json.loads((ROOT / 'ORIGINAL_DATA_MANIFEST.json').read_text())
    args.archive_dir.mkdir(parents=True, exist_ok=True)
    target = args.target.resolve()
    target.mkdir(parents=True, exist_ok=True)
    total = 0
    for archive in manifest['archives']:
        local = args.archive_dir / archive['file']
        if not local.is_file() and args.download:
            print('Downloading', archive['file'], flush=True)
            part = local.with_suffix(local.suffix + '.partial')
            with urllib.request.urlopen(RELEASE + archive['file'], timeout=120) as source, part.open('wb') as sink:
                shutil.copyfileobj(source, sink, length=1024 * 1024)
            if part.stat().st_size != archive['bytes'] or digest(part) != archive['sha256']:
                raise SystemExit('Downloaded archive failed integrity check: ' + archive['file'])
            part.replace(local)
        if not local.is_file():
            raise SystemExit('Missing ' + str(local) + '; use --download or --archive-dir with the release assets.')
        if local.stat().st_size != archive['bytes'] or digest(local) != archive['sha256']:
            raise SystemExit('Archive checksum mismatch: ' + archive['file'])
        expected = {x['file']: x for x in archive['members']}
        with zipfile.ZipFile(local) as z:
            if len(z.namelist()) != len(expected) or set(z.namelist()) != set(expected):
                raise SystemExit('Unexpected archive contents: ' + archive['file'])
            for name, entry in expected.items():
                rel = PurePosixPath(name)
                dest = (target / name).resolve()
                if rel.is_absolute() or '..' in rel.parts or '\\' in name or not dest.is_relative_to(target):
                    raise SystemExit('Unsafe archive member: ' + name)
                raw = z.read(name)
                if len(raw) != entry['bytes'] or hashlib.sha256(raw).hexdigest() != entry['sha256']:
                    raise SystemExit('Member checksum mismatch: ' + name)
                if dest.exists():
                    if digest(dest) != entry['sha256']:
                        raise SystemExit('Refusing to overwrite changed file: ' + str(dest))
                else:
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    dest.write_bytes(raw)
                total += 1
        print('Verified and restored', archive['file'], '(' + str(len(expected)) + ' files)', flush=True)
    print('Restored', total, 'files to', target)

if __name__ == '__main__':
    main()
