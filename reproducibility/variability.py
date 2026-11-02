#!/usr/bin/env python3
"""Intra- and inter-subject variability of pipeline outputs, step by step.

Works on Tractoflow (results/<sub>_<ses>/<Process>/...) and sf-tractomics
(outdir/<sub>/<ses>/{anat,dwi}/...) outputs, for any number of sessions per
subject (e.g. Myelo-inferno: 20 subjects x 5-6 sessions).

Per session (native space):
  scalars  : mean of DTI / fODF maps in the session's own WM mask,
             mask volumes, FRF values, number of streamlines.
Per subject, every session is rigidly registered (ANTs, FA) onto the first
session, then compared to it:
  masks       : Dice
  scalar maps : Pearson r and relative mean |diff| in the reference WM mask
  tractograms : Dice of binarized density maps
Group (one row per pipeline x step x measure):
  CV_intra   : mean over subjects of SD/mean across sessions (%)
  CV_inter   : SD/mean of the subject means (%)
  ICC(1)     : one-way random ICC, unbalanced design (sessions per subject vary)
  ratio      : CV_inter / CV_intra (>1: subjects are distinguishable)
  agreement  : mean (SD) of Dice / r to the reference session

Outputs: sessions.csv, agreement.csv, agreement_summary.csv, group.csv.
Run inside the scilus image (ANTs, scilpy, nibabel, pandas).
"""
import argparse
import os
import re
import subprocess
from concurrent.futures import ProcessPoolExecutor

import nibabel as nib
import numpy as np
import pandas as pd

# canonical name -> (kind, step label, regex on the path relative to the
# session directory)
LAYOUTS = {
    "tractoflow": {
        "fa": ("scalar", "DTI", r"DTI_Metrics/.*__fa\.nii\.gz$"),
        "md": ("scalar", "DTI", r"DTI_Metrics/.*__md\.nii\.gz$"),
        "ad": ("scalar", "DTI", r"DTI_Metrics/.*__ad\.nii\.gz$"),
        "rd": ("scalar", "DTI", r"DTI_Metrics/.*__rd\.nii\.gz$"),
        "ga": ("scalar", "DTI", r"DTI_Metrics/.*__ga\.nii\.gz$"),
        "mode": ("scalar", "DTI", r"DTI_Metrics/.*__mode\.nii\.gz$"),
        "afd_max": ("scalar", "fODF", r"FODF_Metrics/.*__afd_max\.nii\.gz$"),
        "afd_total": ("scalar", "fODF", r"FODF_Metrics/.*__afd_total\.nii\.gz$"),
        "afd_sum": ("scalar", "fODF", r"FODF_Metrics/.*__afd_sum\.nii\.gz$"),
        "nufo": ("scalar", "fODF", r"FODF_Metrics/.*__nufo\.nii\.gz$"),
        "b0_mask": ("mask", "DWI brain mask", r".*__b0_mask_resampled\.nii\.gz$"),
        "t1_mask_dwi": ("mask", "T1 registration", r"Register_T1/.*__t1_mask_warped\.nii\.gz$"),
        "mask_wm": ("mask", "Tissue segmentation", r"Segment_Tissues/.*__mask_wm\.nii\.gz$"),
        "mask_gm": ("mask", "Tissue segmentation", r"Segment_Tissues/.*__mask_gm\.nii\.gz$"),
        "mask_csf": ("mask", "Tissue segmentation", r"Segment_Tissues/.*__mask_csf\.nii\.gz$"),
        "map_include": ("prob", "PFT maps", r"PFT_Tracking_Maps/.*__map_include\.nii\.gz$"),
        "map_exclude": ("prob", "PFT maps", r"PFT_Tracking_Maps/.*__map_exclude\.nii\.gz$"),
        "pft_seeding_mask": ("mask", "PFT seeding", r"PFT_Seeding_Mask/.*\.nii\.gz$"),
        "local_tracking_mask": ("mask", "Local masks", r"Local_Tracking_Mask/.*\.nii\.gz$"),
        "local_seeding_mask": ("mask", "Local masks", r"Local_Seeding_Mask/.*\.nii\.gz$"),
        "frf": ("frf", "FRF", r"Compute_FRF/.*__frf\.txt$"),
        "tracking_pft": ("tractogram", "PFT tracking", r"PFT_Tracking/.*\.trk$"),
        "tracking_local": ("tractogram", "Local tracking", r"Local_Tracking/.*\.trk$"),
    },
    "sftractomics": {
        **{m: ("scalar", "DTI",
               rf"^(?!.*space-)(?!.*desc-fwc)dwi/.*_model-tensor_param-{m}_dwimap\.nii\.gz$")
           for m in ["fa", "md", "ad", "rd", "ga", "mode"]},
        **{m: ("scalar", "fODF",
               rf"^(?!.*space-)dwi/.*_model-csd_param-{m.replace('_', '')}_dwimap\.nii\.gz$")
           for m in ["afd_max", "afd_total", "afd_sum", "nufo"]},
        "b0_mask": ("mask", "DWI brain mask", r"^(?!.*space-)dwi/.*_desc-brain_mask\.nii\.gz$"),
        "t1_mask_dwi": ("mask", "T1 registration", r"^(?!.*space-)anat/.*_desc-brain_mask\.nii\.gz$"),
        "mask_wm": ("mask", "Tissue segmentation", r"anat/.*_label-WM_desc-fastseg_mask\.nii\.gz$"),
        "mask_gm": ("mask", "Tissue segmentation", r"anat/.*_label-GM_desc-fastseg_mask\.nii\.gz$"),
        "mask_csf": ("mask", "Tissue segmentation", r"anat/.*_label-CSF_desc-fastseg_mask\.nii\.gz$"),
        "map_include": ("prob", "PFT maps", r"dwi/.*_label-include_desc-pft_probseg\.nii\.gz$"),
        "map_exclude": ("prob", "PFT maps", r"dwi/.*_label-exclude_desc-pft_probseg\.nii\.gz$"),
        "pft_seeding_mask": ("mask", "PFT seeding", r"dwi/.*_label-seeding_desc-pft_mask\.nii\.gz$"),
        "local_tracking_mask": ("mask", "Local masks", r"dwi/.*_label-tracking_desc-local_mask\.nii\.gz$"),
        "local_seeding_mask": ("mask", "Local masks", r"dwi/.*_label-seeding_desc-local_mask\.nii\.gz$"),
        "frf": ("frf", "FRF", r"dwi/.*_frf\.txt$"),
        "tracking_pft": ("tractogram", "PFT tracking", r"dwi/.*_desc-pft_tractogram\.trk$"),
        "tracking_local": ("tractogram", "Local tracking", r"dwi/.*_desc-local_tractogram\.trk$"),
    },
}


def run(cmd):
    subprocess.run(cmd, shell=True, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def list_sessions(root, layout, sid_regex):
    """Return {subject: {session: session_dir}}."""
    out = {}
    if layout == "tractoflow":
        rx = re.compile(sid_regex)
        for d in sorted(os.listdir(root)):
            m = rx.match(d)
            if m and os.path.isdir(os.path.join(root, d)):
                out.setdefault(m.group(1), {})[m.group(2)] = os.path.join(root, d)
    else:
        for sub in sorted(os.listdir(root)):
            if not sub.startswith("sub-"):
                continue
            for ses in sorted(os.listdir(os.path.join(root, sub))):
                if ses.startswith("ses-"):
                    out.setdefault(sub, {})[ses] = os.path.join(root, sub, ses)
    return {s: v for s, v in out.items() if len(v) > 1}


def find_files(ses_dir, layout):
    rels = []
    for r, _, files in os.walk(ses_dir, followlinks=True):
        rels += [os.path.relpath(os.path.join(r, f), ses_dir) for f in files]
    found = {}
    for name, (_, _, rx) in LAYOUTS[layout].items():
        hits = sorted(f for f in rels if re.search(rx, f))
        if hits:
            found[name] = os.path.join(ses_dir, hits[0])
    return found


def load(path):
    return nib.load(path).get_fdata(dtype=np.float32)


def dice(a, b):
    a, b = a > 0, b > 0
    tot = a.sum() + b.sum()
    return 2 * np.logical_and(a, b).sum() / tot if tot else np.nan


def session_values(files, layout, ref_shape):
    """Native-space summary values of one session."""
    vals = {}
    wm = None
    if "mask_wm" in files and load(files["mask_wm"]).shape == ref_shape:
        wm = load(files["mask_wm"]) > 0
    for name, path in files.items():
        kind = LAYOUTS[layout][name][0]
        if kind == "scalar" and wm is not None:
            d = load(path)
            vals[name] = float(np.nanmean(d[wm]))
        elif kind == "mask":
            d = nib.load(path)
            vals[f"{name}_volume_ml"] = float(
                (np.asanyarray(d.dataobj) > 0).sum() * np.prod(d.header.get_zooms()[:3]) / 1000)
        elif kind == "frf":
            for i, v in enumerate(np.loadtxt(path).ravel()):
                vals[f"frf_{i}"] = float(v)
        elif kind == "tractogram":
            vals[f"{name}_n_streamlines"] = float(
                nib.streamlines.load(path, lazy_load=True).header["nb_streamlines"])
    return vals


def process_subject(job):
    pipeline, layout, sub, sessions, tmp = job
    os.makedirs(tmp, exist_ok=True)
    os.environ["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = "1"
    ses_files = {ses: find_files(d, layout) for ses, d in sessions.items()}
    ref_ses = sorted(sessions)[0]
    ref = ses_files[ref_ses]
    if "fa" not in ref:
        return [], []
    ref_img = nib.load(ref["fa"])
    ref_wm = load(ref["mask_wm"]) > 0 if "mask_wm" in ref and \
        load(ref["mask_wm"]).shape == ref_img.shape else load(ref["fa"]) > 0.2

    def density(path, tag):
        out = os.path.join(tmp, f"{tag}_density.nii.gz")
        if not os.path.exists(out):
            run(f"scil_tractogram_compute_density_map {path} {out} --binary -f")
        return out

    sess_rows, agree_rows = [], []
    for ses, files in ses_files.items():
        if "fa" not in files:
            continue
        shape = nib.load(files["fa"]).shape
        sess_rows.append({"pipeline": pipeline, "subject": sub, "session": ses,
                          **session_values(files, layout, shape)})
        if ses == ref_ses:
            continue
        pre = os.path.join(tmp, f"{ses}_to_{ref_ses}_")
        if not os.path.exists(pre + "0GenericAffine.mat"):
            run(f"antsRegistrationSyNQuick.sh -d 3 -f {ref['fa']} -m {files['fa']} "
                f"-t r -n 1 -o {pre}")

        def warp(img, interp, tag):
            out = os.path.join(tmp, f"{ses}_{tag}_in_ref.nii.gz")
            run(f"antsApplyTransforms -d 3 -i {img} -r {ref['fa']} -o {out} "
                f"-t {pre}0GenericAffine.mat -n {interp}")
            return load(out)

        for name, path in files.items():
            if name not in ref:
                continue
            kind, step, _ = LAYOUTS[layout][name]
            row = {"pipeline": pipeline, "subject": sub, "session": ses,
                   "reference": ref_ses, "step": step, "measure": name}
            try:
                if kind == "mask":
                    rd = load(ref[name])
                    if rd.shape != ref_img.shape:
                        continue  # not in DWI space
                    row["dice"] = dice(rd, warp(path, "NearestNeighbor", name))
                elif kind in ("scalar", "prob"):
                    a, b = load(ref[name]), warp(path, "Linear", name)
                    sel = ref_wm & np.isfinite(a) & np.isfinite(b)
                    row["pearson_r"] = float(np.corrcoef(a[sel], b[sel])[0, 1])
                    row["rel_mean_abs_diff"] = float(
                        np.abs(a[sel] - b[sel]).mean() / np.abs(a[sel]).mean())
                elif kind == "tractogram":
                    a = load(density(ref[name], f"{ref_ses}_{name}"))
                    b = warp(density(path, f"{ses}_{name}"), "NearestNeighbor", name)
                    row["dice"] = dice(a, b)
                else:
                    continue
            except Exception as e:  # keep the other measures
                row["note"] = f"error: {e}"
            agree_rows.append(row)
    return sess_rows, agree_rows


def icc1_unbalanced(df, col):
    g = df.groupby("subject")[col]
    n_i, m_i = g.count(), g.mean()
    keep = n_i > 1
    n_i, m_i = n_i[keep], m_i[keep]
    k, N = len(n_i), n_i.sum()
    if k < 3:
        return np.nan
    x = df[df["subject"].isin(n_i.index)]
    gm = x[col].mean()
    msb = (n_i * (m_i - gm) ** 2).sum() / (k - 1)
    msw = ((x[col] - x["subject"].map(m_i)) ** 2).sum() / (N - k)
    n0 = (N - (n_i ** 2).sum() / N) / (k - 1)
    if msb + (n0 - 1) * msw == 0:
        return np.nan  # constant measure
    return (msb - msw) / (msb + (n0 - 1) * msw)


def group_stats(sess, agree):
    rows = []
    for pipeline, df in sess.groupby("pipeline"):
        for col in [c for c in df.columns
                    if c not in ("pipeline", "subject", "session")]:
            x = df[["subject", col]].dropna()
            g = x.groupby("subject")[col]
            m, s = g.mean(), g.std(ddof=1)
            cv_intra = 100 * (s / m.abs()).mean()
            cv_inter = 100 * m.std(ddof=1) / abs(m.mean())
            rows.append({"pipeline": pipeline, "measure": col,
                         "n_subjects": g.ngroups, "n_sessions": len(x),
                         "mean": x[col].mean(), "CV_intra_%": cv_intra,
                         "CV_inter_%": cv_inter,
                         "ratio_inter_intra": cv_inter / cv_intra if cv_intra else np.nan,
                         "ICC1": icc1_unbalanced(x, col)})
    out = pd.DataFrame(rows)
    if len(agree):
        per_sub = agree.groupby(["pipeline", "step", "measure", "subject"]).mean(
            numeric_only=True)
        agg = per_sub.groupby(["pipeline", "step", "measure"]).agg(["mean", "std"])
        agg.columns = [f"{a}_{b}" for a, b in agg.columns]
        agg = agg.reset_index()
        agg.to_csv("agreement_summary.csv", index=False)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--run", nargs=3, action="append", required=True,
                   metavar=("NAME", "LAYOUT", "DIR"),
                   help="Pipeline output to analyse: name, layout "
                        "(tractoflow|sftractomics), results directory. Repeatable.")
    p.add_argument("--out_dir", required=True)
    p.add_argument("--sid_regex", default=r"^(sub-[^_]+)_(ses-[^_]+)",
                   help="Tractoflow layout: regex with groups (subject, session).")
    p.add_argument("--nproc", type=int, default=1)
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    jobs = []
    for name, layout, root in args.run:
        for sub, sessions in list_sessions(root, layout, args.sid_regex).items():
            jobs.append((name, layout, sub, sessions,
                         os.path.join(args.out_dir, "work", name, sub)))
    print(f"{len(jobs)} subject x pipeline jobs")

    sess, agree = [], []
    with ProcessPoolExecutor(args.nproc) as ex:
        for s, a in ex.map(process_subject, jobs):
            sess += s
            agree += a
    sess, agree = pd.DataFrame(sess), pd.DataFrame(agree)
    os.chdir(args.out_dir)
    sess.to_csv("sessions.csv", index=False)
    agree.to_csv("agreement.csv", index=False)
    group = group_stats(sess, agree)
    group.to_csv("group.csv", index=False)
    print(group.to_string(index=False))


if __name__ == "__main__":
    main()
