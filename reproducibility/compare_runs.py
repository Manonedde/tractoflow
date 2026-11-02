#!/usr/bin/env python3
"""Compare two runs of the same pipeline, file by file and step by step
(computational run-retest reproducibility).

Each run directory must contain the Nextflow trace (trace.txt) and work/.
Every file written by a task (outputs, not staged inputs) is attributed to its
process (= pipeline step) and matched between runs by (step, file name), file
names carrying the subject/session prefix. Works for Tractoflow and
sf-tractomics alike.

Per file:
  - NIfTI : bitwise identical data?, max |diff|, mean |diff| relative to the
            mean |value|, Pearson r (non-zero voxels); Dice for binary masks.
  - TRK/TCK : number of streamlines, identical points?
  - other : sha256 identical?
Writes one CSV row per file, with the step order of the pipeline.
"""
import argparse
import csv
import glob
import hashlib
import os

import numpy as np
import nibabel as nib

SKIP = (".command", ".exitcode", "versions.yml", ".nextflow")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def task_outputs(run_dir):
    """{(step, file): path} and the step order, from trace.txt."""
    files, order = {}, []
    with open(os.path.join(run_dir, "trace.txt")) as f:
        for t in csv.DictReader(f, delimiter="\t"):
            if t["status"] not in ("COMPLETED", "CACHED"):
                continue
            step = t["name"].split(" (")[0].split(":")[-1]
            if step not in order:
                order.append(step)
            for wd in glob.glob(os.path.join(run_dir, "work", t["hash"] + "*")):
                for root, _, names in os.walk(wd):
                    for n in names:
                        p = os.path.join(root, n)
                        if n.startswith(SKIP) or os.path.islink(p):
                            continue  # staged inputs are symlinks
                        files[(step, os.path.relpath(p, wd))] = p
    return files, order


def compare_nifti(a, b):
    ia, ib = nib.load(a), nib.load(b)
    row = {"kind": "nifti"}
    if ia.shape != ib.shape:
        row.update(identical=False, note=f"shape {ia.shape} vs {ib.shape}")
        return row
    da, db = np.asanyarray(ia.dataobj), np.asanyarray(ib.dataobj)
    row["identical"] = bool(np.array_equal(da, db, equal_nan=True)
                            and np.allclose(ia.affine, ib.affine))
    if np.issubdtype(da.dtype, np.integer) and da.max() <= 1 and db.max() <= 1:
        tot = (da > 0).sum() + (db > 0).sum()
        row["dice"] = 2 * np.logical_and(da > 0, db > 0).sum() / tot if tot else 1.0
    fa = np.nan_to_num(da.astype(np.float64))
    fb = np.nan_to_num(db.astype(np.float64))
    diff = np.abs(fa - fb)
    row["max_abs_diff"] = float(diff.max())
    scale = np.abs(fa).mean()
    row["rel_mean_abs_diff"] = float(diff.mean() / scale) if scale else 0.0
    nz = (fa != 0) | (fb != 0)
    if nz.sum() > 2 and fa[nz].std() > 0 and fb[nz].std() > 0:
        row["pearson_r"] = float(np.corrcoef(fa[nz], fb[nz])[0, 1])
    return row


def compare_tractogram(a, b):
    ta = nib.streamlines.load(a, lazy_load=False).streamlines
    tb = nib.streamlines.load(b, lazy_load=False).streamlines
    row = {"kind": "tractogram",
           "n_streamlines_a": len(ta), "n_streamlines_b": len(tb)}
    row["identical"] = (ta._data.shape == tb._data.shape
                        and np.array_equal(ta._data, tb._data))
    return row


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("run_a", help="Run directory (trace.txt + work/).")
    p.add_argument("run_b")
    p.add_argument("out_csv")
    args = p.parse_args()

    fa, order = task_outputs(args.run_a)
    fb, _ = task_outputs(args.run_b)
    rows = []
    for (step, name), pa in sorted(fa.items(), key=lambda x: (order.index(x[0][0]), x[0][1])):
        base = {"step": step, "step_order": order.index(step), "file": name}
        pb = fb.get((step, name))
        if pb is None:
            rows.append({**base, "note": "missing in run_b"})
            continue
        try:
            if name.endswith((".nii.gz", ".nii")):
                res = compare_nifti(pa, pb)
            elif name.endswith((".trk", ".tck")):
                res = compare_tractogram(pa, pb)
            else:
                res = {"kind": "other", "identical": sha256(pa) == sha256(pb)}
        except Exception as e:  # keep going, report the failure
            res = {"note": f"error: {e}"}
        rows.append({**base, **res})

    fields = ["step", "step_order", "file", "kind", "identical", "dice",
              "max_abs_diff", "rel_mean_abs_diff", "pearson_r",
              "n_streamlines_a", "n_streamlines_b", "note"]
    with open(args.out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    n_id = sum(1 for r in rows if r.get("identical") is True)
    print(f"{n_id}/{len(rows)} files identical")
    for step in order:
        rs = [r for r in rows if r["step"] == step]
        if rs:
            k = sum(1 for r in rs if r.get("identical") is True)
            print(f"  {step:45s} {k:4d}/{len(rs):<4d} identical")


if __name__ == "__main__":
    main()
