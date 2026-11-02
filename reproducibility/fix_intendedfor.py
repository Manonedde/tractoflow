#!/usr/bin/env python3
"""Add IntendedFor to reverse phase-encoded b0 sidecars, so that Tractoflow
(scil_bids_validate) finds the reverse b0 used by topup.

For every subject/session, each DWI (dwi/*_dwi.nii.gz) gets exactly one
reverse b0: a file with the opposite PhaseEncodingDirection among
fmap/*_epi, fmap/*_sbref, dwi/*_sbref (in this order of preference, same run
entity first). Its JSON gets "IntendedFor": "<ses-X/>dwi/<dwi file>" (path
relative to the subject directory, as Tractoflow expects). Existing entries
are kept. Full reverse DWIs (*_dwi with bval) need no IntendedFor.

  fix_intendedfor.py <bids_dir>                    # report only
  fix_intendedfor.py <bids_dir> --out_dir <dir>    # BIDS mirror: symlinks to
                                                   # the original files, edited
                                                   # JSONs written (source untouched)
  fix_intendedfor.py <bids_dir> --in_place         # edit JSONs, keep <json>.orig

With a single DWI, Tractoflow only accepts a reverse b0 stored as fmap/*_epi
(a reverse sbref is rejected). --sbref_as_epi (mirror only) exposes such a
reverse sbref as fmap/<name>_epi.nii.gz (symlink) with its IntendedFor.
"""
import argparse
import glob
import json
import os
import re
import shutil

OPPOSITE = {"i": "i-", "i-": "i", "j": "j-", "j-": "j", "k": "k-", "k-": "k"}
REV_PATTERNS = ["fmap/*_epi.nii*", "fmap/*_sbref.nii*", "dwi/*_sbref.nii*"]


def sidecar(nii):
    return re.sub(r"\.nii(\.gz)?$", ".json", nii)


def load(path):
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def entity(name, key):
    m = re.search(rf"_{key}-([a-zA-Z0-9]+)", os.path.basename(name))
    return m.group(1) if m else None


def plan(bids, sbref_as_epi=False):
    """{json_path: [intended paths]}, {epi_to_create: (src_nii, json)}, report."""
    edits, as_epi, rows = {}, {}, []
    for sub_dir in sorted(glob.glob(os.path.join(bids, "sub-*"))):
        ses_dirs = sorted(glob.glob(os.path.join(sub_dir, "ses-*"))) or [sub_dir]
        for d in ses_dirs:
            label = os.path.relpath(d, bids)
            dwis = [f for f in sorted(glob.glob(os.path.join(d, "dwi", "*_dwi.nii*")))
                    if os.path.exists(re.sub(r"\.nii(\.gz)?$", ".bval", f))]
            revs = [f for p in REV_PATTERNS for f in sorted(glob.glob(os.path.join(d, p)))]
            pe = {f: load(sidecar(f)).get("PhaseEncodingDirection") for f in dwis + revs}
            # DWIs with a full reverse DWI acquisition are already handled
            rev_dwi_pe = {pe[f] for f in dwis}
            for dwi in dwis:
                rel = os.path.relpath(dwi, sub_dir)
                if pe[dwi] is None:
                    rows.append((label, os.path.basename(dwi), "-", "NO PhaseEncodingDirection in DWI json"))
                    continue
                if OPPOSITE[pe[dwi]] in rev_dwi_pe and len(dwis) > 1:
                    rows.append((label, os.path.basename(dwi), "reverse DWI", "ok (no IntendedFor needed)"))
                    continue
                cands = [r for r in revs if pe[r] == OPPOSITE[pe[dwi]]]
                if not cands:
                    rows.append((label, os.path.basename(dwi), "-", "NO reverse b0 found"))
                    continue
                same_run = [r for r in cands if entity(r, "run") == entity(dwi, "run")]
                rev = (same_run or cands)[0]
                if rev.endswith(("_sbref.nii", "_sbref.nii.gz")) and sbref_as_epi:
                    epi = os.path.join(d, "fmap", re.sub(r"_sbref\.nii(\.gz)?$", "_epi.nii.gz",
                                                         os.path.basename(rev)))
                    data = as_epi.setdefault(epi, (rev, {**load(sidecar(rev)), "IntendedFor": []}))[1]
                    data["IntendedFor"].append(rel)
                    rows.append((label, os.path.basename(dwi), os.path.relpath(epi, d),
                                 "sbref exposed as fmap epi + IntendedFor"))
                    continue
                js = sidecar(rev)
                current = load(js).get("IntendedFor", [])
                current = [current] if isinstance(current, str) else current
                if rel in current:
                    status = "ok (already set)"
                else:
                    edits.setdefault(js, list(current)).append(rel)
                    status = "IntendedFor added"
                if len(cands) > 1:
                    status += f" ({len(cands)} candidates, first kept)"
                if "_sbref" in rev and len(dwis) == 1:
                    status += " | WARNING: Tractoflow ignores a reverse sbref, use --sbref_as_epi"
                rows.append((label, os.path.basename(dwi), os.path.relpath(rev, d), status))
    return edits, as_epi, rows


def mirror(bids, out, edits, as_epi):
    for root, dirs, files in os.walk(bids):
        dst_root = os.path.join(out, os.path.relpath(root, bids))
        os.makedirs(dst_root, exist_ok=True)
        for f in files:
            src, dst = os.path.join(root, f), os.path.join(dst_root, f)
            if os.path.lexists(dst):
                os.remove(dst)
            if src in edits:
                data = load(src)
                data["IntendedFor"] = edits[src]
                with open(dst, "w") as fh:
                    json.dump(data, fh, indent=4)
            else:
                os.symlink(os.path.abspath(src), dst)
    for epi, (src, data) in as_epi.items():
        dst = os.path.join(out, os.path.relpath(epi, bids))
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        if os.path.lexists(dst):
            os.remove(dst)
        os.symlink(os.path.abspath(src), dst)
        with open(sidecar(dst), "w") as fh:
            json.dump(data, fh, indent=4)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("bids_dir")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--out_dir", help="Write a BIDS mirror with the edited JSONs.")
    g.add_argument("--in_place", action="store_true",
                   help="Edit the JSONs in place (backup: <json>.orig).")
    p.add_argument("--sbref_as_epi", action="store_true",
                   help="With --out_dir: expose reverse sbref as fmap/*_epi.")
    args = p.parse_args()
    if args.sbref_as_epi and not args.out_dir:
        p.error("--sbref_as_epi requires --out_dir")

    bids = os.path.abspath(args.bids_dir)
    edits, as_epi, rows = plan(bids, args.sbref_as_epi)
    w = max([len(r[0]) for r in rows] + [10])
    for r in rows:
        print(f"{r[0]:{w}s}  {r[1]:45s}  {r[2]:40s}  {r[3]}")
    print(f"\n{len(rows)} DWI, {len(edits)} JSON to edit, {len(as_epi)} sbref as epi, "
          f"{sum('NO ' in r[3] for r in rows)} problem(s)")

    if args.out_dir:
        mirror(bids, os.path.abspath(args.out_dir), edits, as_epi)
        print(f"Mirror written: {args.out_dir}")
    elif args.in_place:
        for js, intended in edits.items():
            shutil.copy2(js, js + ".orig")
            data = load(js)
            data["IntendedFor"] = intended
            with open(js, "w") as fh:
                json.dump(data, fh, indent=4)
        print("JSONs edited in place (.orig backups).")
    else:
        print("Report only: use --out_dir (recommended) or --in_place to write.")


if __name__ == "__main__":
    main()
