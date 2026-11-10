"""Verify and unpack the research data archives; uses only the standard library."""
from pathlib import Path, PurePosixPath
import json, hashlib, zipfile, shutil
ROOT=Path(__file__).resolve().parent
manifest=json.loads((ROOT/'DATA_MANIFEST.json').read_text())
for archive in manifest['archives']:
    path=ROOT/archive['file']
    if not path.is_file():raise SystemExit(f"Missing {path.name}; download the complete repository.")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=archive['sha256']:
        raise SystemExit(f"Checksum mismatch: {path.name}")
    expected={r['file']:r for r in archive['members']}
    with zipfile.ZipFile(path) as z:
        if set(z.namelist())!=set(expected):raise SystemExit(f"Unexpected archive members: {path.name}")
        for name in z.namelist():
            rel=PurePosixPath(name)
            if rel.is_absolute() or '..' in rel.parts or '\\' in name:
                raise SystemExit(f"Unsafe archive path: {name}")
            dest=(ROOT/name).resolve()
            if not dest.is_relative_to(ROOT):raise SystemExit(f"Unsafe destination: {name}")
            raw=z.read(name)
            if hashlib.sha256(raw).hexdigest()!=expected[name]['sha256']:
                raise SystemExit(f"Member checksum mismatch: {name}")
            dest.parent.mkdir(parents=True,exist_ok=True)
            if dest.exists() and dest.read_bytes()!=raw:
                raise SystemExit(f"Refusing to overwrite changed file: {name}. Preserve or move it first.")
            if not dest.exists():dest.write_bytes(raw)
    print(f"Verified and unpacked {path.name}: {len(expected)} files")
