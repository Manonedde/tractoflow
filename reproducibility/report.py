#!/usr/bin/env python3
"""QC images and figures for the Tractoflow / sf-tractomics study.

  report.py qc    --tractoflow DIR --sftractomics DIR --out_dir DIR [--nproc N]
      Per session, Tractoflow vs sf-tractomics side by side (FA, MD, AFD,
      NuFO, masks, tractogram densities) + difference maps, and agreement
      between pipelines (Dice, r, Bland-Altman of WM means) in pipelines.csv.
      sf-tractomics images are resampled on the Tractoflow grid (same DWI, same
      world space: no registration).

  report.py plots --analysis_dir DIR
      Figures + index.html from the outputs of compare_runs.py (run-retest),
      variability.py (intra/inter-subject) and `report.py qc`.

Subjects are shown as S01, S02, ... (mapping in subject_codes.csv).
Run inside the scilus image (matplotlib, nibabel, pandas, scilpy).
"""
import argparse
import glob
import html
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import nibabel as nib  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from nibabel.processing import resample_from_to  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from variability import find_files, list_sessions  # noqa: E402

# Reference palette (dataviz skill): light surface, ink, categorical slots 1-2
SURFACE, INK, INK2, MUTED = "#f3f372", "#0b0b0b", "#52514e", "#4D411E"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
COLORS = {"tractoflow": "#204d85", "sftractomics": "#d47e0e"}
LABELS = {"tractoflow": "Tractoflow", "sftractomics": "sf-tractomics"}
DIVERGING = LinearSegmentedColormap.from_list(
    "div", ["#2a78d6", "#f0efec", "#e34948"])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "text.color": INK, "axes.labelcolor": INK2,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.edgecolor": AXIS,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "grid.linestyle": "-", "axes.axisbelow": True, "axes.spines.top": False,
    "axes.spines.right": False, "lines.linewidth": 2, "font.size": 10,
    "axes.titlesize": 11, "axes.titleweight": "bold", "legend.frameon": False,
    "font.family": "sans-serif",
})

MAPS = [("fa", "FA", "gray", (0, 1)), ("md", "MD", "gray", None),
        ("afd_total", "AFD total", "gray", None), ("nufo", "NuFO", "gray", (0, 5))]
MASKS = ["b0_mask", "mask_wm", "mask_gm", "mask_csf", "pft_seeding_mask",
         "local_tracking_mask", "local_seeding_mask", "t1_mask_dwi"]
TRACKINGS = ["tracking_pft", "tracking_local"]


# ----------------------------------------------------------------- qc ------
def on_grid(path, ref, order):
    img = nib.load(path)
    if img.ndim > 3:
        img = img.slicer[..., 0]
    return np.asarray(resample_from_to(img, ref, order=order).dataobj, dtype=np.float32)


def density(trk, out):
    if not os.path.exists(out):
        subprocess.run(f"scil_tractogram_compute_density_map {trk} {out} -f",
                       shell=True, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    return out


def n_streamlines(path):
    return int(nib.streamlines.load(path, lazy_load=True).header["nb_streamlines"])


def dice(a, b):
    a, b = a > 0, b > 0
    tot = a.sum() + b.sum()
    return 2 * np.logical_and(a, b).sum() / tot if tot else np.nan


def qc_session(job):
    sub, ses, code, tf_dir, sf_dir, out_dir = job
    tf, sf = find_files(tf_dir, "tractoflow"), find_files(sf_dir, "sftractomics")
    if "fa" not in tf or "fa" not in sf:
        return []
    ref_img = nib.load(tf["fa"])
    ref = (ref_img.shape[:3], ref_img.affine)
    wm = nib.load(tf["mask_wm"]).get_fdata() > 0 if "mask_wm" in tf else \
        ref_img.get_fdata() > 0.2
    tmp = os.path.join(out_dir, "work", f"{sub}_{ses}")
    os.makedirs(tmp, exist_ok=True)

    rows, panels = [], []
    base = {"subject": code, "session": ses}
    for key, label, cmap, lim in MAPS:
        if key not in tf or key not in sf:
            continue
        a, b = on_grid(tf[key], ref, 1), on_grid(sf[key], ref, 1)
        sel = wm & np.isfinite(a) & np.isfinite(b)
        rows.append({**base, "kind": "map", "measure": key,
                     "pearson_r_wm": float(np.corrcoef(a[sel], b[sel])[0, 1]),
                     "wm_mean_tractoflow": float(a[sel].mean()),
                     "wm_mean_sftractomics": float(b[sel].mean())})
        panels.append((label, a, b, cmap, lim))
    for key in MASKS:
        if key in tf and key in sf:
            rows.append({**base, "kind": "mask", "measure": key,
                         "dice": dice(on_grid(tf[key], ref, 0), on_grid(sf[key], ref, 0))})
    dens = []
    for key in TRACKINGS:
        if key not in tf or key not in sf:
            continue
        a = on_grid(density(tf[key], f"{tmp}/{key}_tf.nii.gz"), ref, 0)
        b = on_grid(density(sf[key], f"{tmp}/{key}_sf.nii.gz"), ref, 0)
        rows.append({**base, "kind": "tractogram", "measure": key, "dice": dice(a, b),
                     "n_tractoflow": n_streamlines(tf[key]),
                     "n_sftractomics": n_streamlines(sf[key])})
        dens.append((key.replace("tracking_", "").upper() + " density (log)",
                     np.log1p(a), np.log1p(b)))

    # Figure: one row per map, columns Tractoflow | sf-tractomics | difference
    fa = on_grid(tf["fa"], ref, 1)
    mask = on_grid(tf["b0_mask"], ref, 0) > 0 if "b0_mask" in tf else fa > 0
    z = int(np.round(np.argwhere(mask).mean(0)[2])) if mask.any() else fa.shape[2] // 2
    n_rows = len(panels) + len(dens) + 1
    fig, axes = plt.subplots(n_rows, 3, figsize=(9, 3 * n_rows))
    for ax in axes.ravel():
        ax.axis("off")
    for c, t in enumerate(["Tractoflow", "sf-tractomics", "sf-tractomics − Tractoflow"]):
        axes[0, c].set_title(t, color=INK)
    r = 0
    for label, a, b, cmap, lim in panels + [(l, a, b, "magma", None) for l, a, b in dens]:
        vmin, vmax = lim or (0, np.nanpercentile(np.r_[a[mask], b[mask]], 99))
        for c, d in enumerate([a, b]):
            axes[r, c].imshow(np.rot90(d[:, :, z]), cmap=cmap, vmin=vmin, vmax=vmax)
        diff = np.where(mask, b - a, np.nan)
        m = np.nanpercentile(np.abs(diff[mask]), 99) or 1
        im = axes[r, 2].imshow(np.rot90(diff[:, :, z]), cmap=DIVERGING, vmin=-m, vmax=m)
        fig.colorbar(im, ax=axes[r, 2], fraction=0.046, pad=0.02)
        axes[r, 0].text(-0.05, 0.5, label, transform=axes[r, 0].transAxes,
                        rotation=90, ha="right", va="center", color=INK2)
        r += 1
    # Masks: FA underlay, Tractoflow contour (blue) vs sf-tractomics (orange)
    for c, key in enumerate(["b0_mask", "mask_wm", "local_seeding_mask"]):
        ax = axes[r, c]
        ax.imshow(np.rot90(fa[:, :, z]), cmap="gray", vmin=0, vmax=1)
        for p, files, lw in (("tractoflow", tf, 3), ("sftractomics", sf, 1)):
            if key in files:
                ax.contour(np.rot90(on_grid(files[key], ref, 0)[:, :, z]),
                           levels=[0.5], colors=COLORS[p], linewidths=lw)
        ax.set_title(key, color=INK2, fontsize=9)
    axes[r, 0].text(-0.05, 0.5, "Masks", transform=axes[r, 0].transAxes, rotation=90,
                    ha="right", va="center", color=INK2)
    fig.suptitle(f"{code} {ses} — axial slice {z}\nmask contours: Tractoflow (blue, wide), "
                 f"sf-tractomics (orange, thin)", color=INK)
    fig.tight_layout(rect=(0.03, 0, 1, 0.97))
    os.makedirs(os.path.join(out_dir, "sessions"), exist_ok=True)
    fig.savefig(os.path.join(out_dir, "sessions", f"{code}_{ses}.png"), dpi=110)
    plt.close(fig)
    return rows


def cmd_qc(args):
    tf = list_sessions(args.tractoflow, "tractoflow", args.sid_regex)
    sf = list_sessions(args.sftractomics, "sftractomics", None)
    subs = sorted(set(tf) & set(sf))
    codes = {s: f"S{i + 1:02d}" for i, s in enumerate(subs)}
    os.makedirs(args.out_dir, exist_ok=True)
    pd.DataFrame({"subject": subs, "code": [codes[s] for s in subs]}).to_csv(
        os.path.join(args.out_dir, "subject_codes.csv"), index=False)
    jobs = [(s, ses, codes[s], tf[s][ses], sf[s][ses], args.out_dir)
            for s in subs for ses in sorted(set(tf[s]) & set(sf[s]))]
    rows = []
    with ProcessPoolExecutor(args.nproc) as ex:
        for r in ex.map(qc_session, jobs):
            rows += r
    pd.DataFrame(rows).to_csv(os.path.join(args.out_dir, "pipelines.csv"), index=False)
    print(f"{len(jobs)} sessions compared")


# -------------------------------------------------------------- plots ------
def save(fig, out_dir, name):
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, name), dpi=130)
    plt.close(fig)
    return name


def legend_pipelines(ax, pipelines):
    for p in pipelines:
        ax.scatter([], [], color=COLORS.get(p, INK2), s=40, label=LABELS.get(p, p))
    ax.legend(loc="best", fontsize=9)


def plot_runretest(A, out):
    """Per step: % identical files and size of the differences, A vs B."""
    figs = []
    csvs = sorted(glob.glob(os.path.join(A, "runretest_*_*.csv")))
    if not csvs:
        return figs, None
    df = pd.concat([pd.read_csv(f).assign(
        pipeline=os.path.basename(f).split("_")[1],
        config=os.path.basename(f).split("_")[2].replace(".csv", "")) for f in csvs])
    df["identical"] = df["identical"].astype(str).str.lower() == "true"
    pipelines = [p for p in COLORS if p in df.pipeline.unique()]
    configs = [c for c in ("repro", "mt") if c in df.config.unique()]
    cfg_label = {"repro": "reproducible profile", "mt": "default multi-threading"}

    summ = (df.groupby(["pipeline", "config", "step", "step_order"])
            .agg(n_files=("file", "size"), pct_identical=("identical", "mean"),
                 max_rel_diff=("rel_mean_abs_diff", "max"))
            .reset_index())
    summ["pct_identical"] *= 100

    n_steps = summ.groupby("pipeline").step.nunique().max()
    fig, axes = plt.subplots(len(configs), len(pipelines), squeeze=False,
                             figsize=(5.5 * len(pipelines), 0.28 * n_steps * len(configs) + 1.5))
    for i, c in enumerate(configs):
        for j, p in enumerate(pipelines):
            ax = axes[i, j]
            s = summ[(summ.pipeline == p) & (summ.config == c)].sort_values("step_order")
            y = np.arange(len(s))
            ax.barh(y, s.pct_identical, height=0.6, color=COLORS[p])
            for yy, v in zip(y, s.pct_identical):
                if v < 100:  # direct label only where reproducibility breaks
                    ax.text(v + 1, yy, f"{v:.0f}%", va="center", fontsize=8, color=INK)
            ax.set_yticks(y, s.step, fontsize=8)
            ax.invert_yaxis()
            ax.set_xlim(0, 110)
            ax.set_xlabel("% of files bitwise identical (run A vs run B)")
            ax.set_title(f"{LABELS[p]} — {cfg_label[c]}")
            ax.grid(axis="y", visible=False)
    figs.append(("Bitwise identical files per step", save(fig, out, "rr_identical.png"),
                 "Horizontal bars: percentage of the files of each step that are "
                 "bitwise identical between two runs on the same data."))

    nd = df[(~df.identical) & (df.rel_mean_abs_diff > 0)]
    if len(nd):
        fig, axes = plt.subplots(len(configs), len(pipelines), squeeze=False,
                                 figsize=(5.5 * len(pipelines), 0.28 * n_steps * len(configs) + 1.5))
        for i, c in enumerate(configs):
            for j, p in enumerate(pipelines):
                ax = axes[i, j]
                order = summ[(summ.pipeline == p) & (summ.config == c)].sort_values(
                    "step_order").step.tolist()
                s = nd[(nd.pipeline == p) & (nd.config == c)]
                for k, step in enumerate(order):
                    v = s[s.step == step].rel_mean_abs_diff
                    if len(v):
                        ax.scatter(v, np.full(len(v), k) + np.random.uniform(-.2, .2, len(v)),
                                   s=14, color=COLORS[p], alpha=0.7, linewidths=0)
                ax.set_yticks(range(len(order)), order, fontsize=8)
                ax.invert_yaxis()
                ax.set_xscale("log")
                ax.set_xlabel("relative mean |A − B| (non-identical files)")
                ax.set_title(f"{LABELS[p]} — {cfg_label[c]}")
        figs.append(("Size of run-to-run differences", save(fig, out, "rr_differences.png"),
                     "One dot per non-identical image: mean absolute difference between "
                     "runs relative to the mean value, log scale."))
    summ.to_csv(os.path.join(out, "runretest_summary.csv"), index=False)
    return figs, summ


def plot_variability(A, out, codes):
    figs = []
    V = os.path.join(A, "variability")
    if not os.path.exists(os.path.join(V, "group.csv")):
        return figs, None
    group = pd.read_csv(os.path.join(V, "group.csv"))
    group = group[(group["CV_intra_%"] > 0) | (group["CV_inter_%"] > 0)]  # drop constants
    sess = pd.read_csv(os.path.join(V, "sessions.csv"))
    agree = pd.read_csv(os.path.join(V, "agreement.csv"))
    sess["code"] = sess.subject.map(codes).fillna(sess.subject)
    pipelines = [p for p in COLORS if p in group.pipeline.unique()]
    measures = list(dict.fromkeys(group.measure))
    off = {p: (k - (len(pipelines) - 1) / 2) * 0.25 for k, p in enumerate(pipelines)}

    # CV intra vs inter: two panels, one axis each, same measures
    fig, axes = plt.subplots(1, 2, figsize=(11, 0.32 * len(measures) + 1.5), sharey=True)
    for ax, col, title in zip(axes, ["CV_intra_%", "CV_inter_%"],
                              ["Intra-subject CV (across sessions)", "Inter-subject CV"]):
        for p in pipelines:
            g = group[group.pipeline == p].set_index("measure").reindex(measures)
            ax.scatter(g[col], np.arange(len(measures)) + off[p], s=40, color=COLORS[p],
                       edgecolors=SURFACE, linewidths=2, zorder=3)
        ax.set_xlim(left=0)
        ax.set_xlabel("CV (%)")
        ax.set_title(title)
    axes[0].set_yticks(range(len(measures)), measures, fontsize=8)
    axes[0].invert_yaxis()
    legend_pipelines(axes[1], pipelines)
    figs.append(("Intra- vs inter-subject variability", save(fig, out, "var_cv.png"),
                 "Coefficient of variation of each measure, across the sessions of a "
                 "subject (left) and between subject means (right)."))

    # ICC(1)
    fig, ax = plt.subplots(figsize=(6.5, 0.32 * len(measures) + 1.5))
    for p in pipelines:
        g = group[group.pipeline == p].set_index("measure").reindex(measures)
        ax.scatter(g.ICC1, np.arange(len(measures)) + off[p], s=40, color=COLORS[p],
                   edgecolors=SURFACE, linewidths=2, zorder=3)
    ax.axvline(0.75, color=MUTED, linewidth=1)
    ax.text(0.75, -0.9, " 0.75", color=MUTED, fontsize=8)
    ax.set_yticks(range(len(measures)), measures, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(min(0, group.ICC1.min()) - 0.05, 1.02)
    ax.set_xlabel("ICC(1) — share of variance explained by the subject")
    ax.set_title("Reliability across sessions")
    legend_pipelines(ax, pipelines)
    figs.append(("ICC per measure", save(fig, out, "var_icc.png"),
                 "One-way random ICC (unbalanced, 5-6 sessions per subject)."))

    # Boxplots per subject (sessions as dots): intra = box height, inter = spread
    box_measures = [m for m in ["fa", "md", "ad", "rd", "afd_total", "nufo",
                                "b0_mask_volume_ml", "mask_wm_volume_ml",
                                "tracking_pft_n_streamlines", "tracking_local_n_streamlines",
                                "frf_0"] if m in sess.columns]
    subj = sorted(sess.code.unique())
    for m in box_measures:
        fig, ax = plt.subplots(figsize=(max(7, 0.45 * len(subj) * len(pipelines)), 3.6))
        for p in pipelines:
            s = sess[sess.pipeline == p]
            data = [s[s.code == c][m].dropna().values for c in subj]
            pos = np.arange(len(subj)) + off[p] * 1.4
            ax.boxplot(data, positions=pos, widths=0.3, showfliers=False,
                       medianprops=dict(color=COLORS[p], linewidth=2),
                       boxprops=dict(color=COLORS[p], linewidth=1),
                       whiskerprops=dict(color=COLORS[p], linewidth=1),
                       capprops=dict(color=COLORS[p], linewidth=1))
            for x, d in zip(pos, data):
                ax.scatter(np.full(len(d), x), d, s=10, color=COLORS[p], alpha=0.6,
                           linewidths=0, zorder=3)
        ax.set_xticks(range(len(subj)), subj, fontsize=8)
        ax.set_ylabel(m)
        ax.set_title(f"{m}: sessions per subject")
        ax.grid(axis="x", visible=False)
        legend_pipelines(ax, pipelines)
        figs.append((f"Per-subject distribution — {m}", save(fig, out, f"var_box_{m}.png"),
                     "Box per subject over its sessions (dots); box height shows "
                     "intra-subject variability, spread of boxes inter-subject variability."))

    # Histograms of within-subject agreement (Dice / r) per step
    for col, xl in (("dice", "Dice with the subject's first session"),
                    ("pearson_r", "voxelwise r with the subject's first session (WM)")):
        a = agree.dropna(subset=[col]) if col in agree else pd.DataFrame()
        if not len(a):
            continue
        ms = list(dict.fromkeys(a.measure))
        ncol = 3
        nrow = int(np.ceil(len(ms) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(11, 2.4 * nrow), squeeze=False)
        for ax in axes.ravel()[len(ms):]:
            ax.axis("off")
        for ax, m in zip(axes.ravel(), ms):
            vals = a[a.measure == m][col]
            bins = np.linspace(vals.min(), 1, 25) if vals.min() < 1 else 10
            for p in pipelines:
                v = a[(a.measure == m) & (a.pipeline == p)][col]
                ax.hist(v, bins=bins, histtype="step", color=COLORS[p], linewidth=2)
            ax.set_title(m, fontsize=9)
            ax.set_xlabel(col, fontsize=8)
        legend_pipelines(axes.ravel()[len(ms) - 1], pipelines)
        name = f"var_hist_{col}.png"
        figs.append((f"Histogram — {xl}", save(fig, out, name),
                     "Distribution over all (subject, session) pairs; masks and "
                     "tractograms use Dice, scalar maps use r."))

    # Session trajectories (FA WM mean)
    if "fa" in sess:
        fig, axes = plt.subplots(1, len(pipelines), figsize=(5.5 * len(pipelines), 3.4),
                                 sharey=True, squeeze=False)
        for ax, p in zip(axes[0], pipelines):
            s = sess[sess.pipeline == p].sort_values("session")
            for c, g in s.groupby("code"):
                ax.plot(range(len(g)), g.fa, color=AXIS, linewidth=1)
            idx = s.groupby("code").cumcount()
            ax.plot(sorted(idx.unique()), s.groupby(idx).fa.mean(), color=COLORS[p])
            ax.set_xlabel("session index")
            ax.set_title(f"{LABELS[p]}: FA (WM mean) per session")
        axes[0, 0].set_ylabel("FA")
        figs.append(("Session trajectories — FA", save(fig, out, "var_traj_fa.png"),
                     "Gray: one line per subject; colored: mean over subjects."))
    return figs, group


def plot_pipelines(A, out):
    figs = []
    f = os.path.join(A, "qc", "pipelines.csv")
    if not os.path.exists(f):
        return figs, None
    df = pd.read_csv(f)
    for col, kinds, xl in (("dice", ["mask", "tractogram"], "Dice Tractoflow vs sf-tractomics"),
                           ("pearson_r_wm", ["map"], "voxelwise r in WM, Tractoflow vs sf-tractomics")):
        s = df[df.kind.isin(kinds)].dropna(subset=[col])
        if not len(s):
            continue
        ms = list(dict.fromkeys(s.measure))
        fig, ax = plt.subplots(figsize=(7, 0.4 * len(ms) + 1.5))
        ax.boxplot([s[s.measure == m][col] for m in ms], vert=False, widths=0.5,
                   showfliers=False, medianprops=dict(color=INK, linewidth=2),
                   boxprops=dict(color=INK2), whiskerprops=dict(color=INK2),
                   capprops=dict(color=INK2))
        for k, m in enumerate(ms):
            v = s[s.measure == m][col]
            ax.scatter(v, np.full(len(v), k + 1) + np.random.uniform(-.15, .15, len(v)),
                       s=10, color=INK2, alpha=0.5, linewidths=0)
        ax.set_yticks(range(1, len(ms) + 1), ms, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel(xl)
        ax.set_title("Agreement between pipelines (one dot per session)")
        figs.append((xl, save(fig, out, f"pipe_{col}.png"), "Same session, both pipelines."))
    maps = df[df.kind == "map"]
    if len(maps):
        ms = list(dict.fromkeys(maps.measure))
        fig, axes = plt.subplots(1, len(ms), figsize=(3.4 * len(ms), 3.2), squeeze=False)
        for ax, m in zip(axes[0], ms):
            s = maps[maps.measure == m]
            mean = (s.wm_mean_tractoflow + s.wm_mean_sftractomics) / 2
            d = s.wm_mean_sftractomics - s.wm_mean_tractoflow
            ax.scatter(mean, d, s=14, color=INK2, linewidths=0)
            for y in (d.mean(), d.mean() - 1.96 * d.std(), d.mean() + 1.96 * d.std()):
                ax.axhline(y, color=MUTED if y != d.mean() else INK, linewidth=1)
            ax.axhline(0, color=AXIS, linewidth=1)
            ax.set_title(m)
            ax.set_xlabel("mean of pipelines (WM)")
        axes[0, 0].set_ylabel("sf-tractomics − Tractoflow")
        figs.append(("Bland-Altman of WM means", save(fig, out, "pipe_bland_altman.png"),
                     "Bias (dark line) and 95% limits of agreement (gray lines)."))
    return figs, df


def table(df, max_rows=200):
    return df.head(max_rows).to_html(index=False, float_format=lambda x: f"{x:.4g}",
                                     border=0, classes="t")


def cmd_plots(args):
    A = args.analysis_dir
    out = os.path.join(A, "report")
    os.makedirs(out, exist_ok=True)
    cf = os.path.join(A, "qc", "subject_codes.csv")
    codes = dict(pd.read_csv(cf).values) if os.path.exists(cf) else {}

    sections = []
    figs, summ = plot_runretest(A, out)
    sections.append(("1. Computational reproducibility (run A vs run B)", figs,
                     table(summ) if summ is not None else ""))
    figs, group = plot_variability(A, out, codes)
    sections.append(("2. Intra- and inter-subject variability (Myelo-inferno)", figs,
                     table(group) if group is not None else ""))
    figs, pipe = plot_pipelines(A, out)
    qc_imgs = sorted(glob.glob(os.path.join(A, "qc", "sessions", "*.png")))
    links = "".join(f'<li><a href="../qc/sessions/{os.path.basename(p)}">'
                    f'{html.escape(os.path.basename(p)[:-4])}</a></li>' for p in qc_imgs)
    sections.append(("3. Tractoflow vs sf-tractomics", figs,
                     (table(pipe.groupby(["kind", "measure"]).mean(numeric_only=True)
                            .reset_index()) if pipe is not None else "")
                     + (f"<h3>QC images per session</h3><ul class='qc'>{links}</ul>"
                        if links else "")))

    body = ""
    for title, figs, extra in sections:
        body += f"<section><h2>{html.escape(title)}</h2>"
        if not figs and not extra:
            body += "<p class='muted'>No data yet.</p>"
        for t, img, alt in figs:
            body += (f"<figure><img src='{img}' alt='{html.escape(alt)}'>"
                     f"<figcaption><b>{html.escape(t)}</b> — {html.escape(alt)}"
                     f"</figcaption></figure>")
        body += f"{extra}</section>"
    page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Reproducibility report</title><style>
:root{{--bg:#f9f9f7;--surface:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--muted:#898781;--line:#e1e0d9}}
body{{margin:0;padding:24px 16px;background:var(--bg);color:var(--ink);
font:15px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif}}
main{{max-width:1150px;margin:auto}} section{{background:var(--surface);border:1px solid var(--line);
border-radius:8px;padding:16px 20px;margin:20px 0}} h2{{margin-top:0}}
figure{{margin:16px 0}} img{{max-width:100%;height:auto}} figcaption{{color:var(--ink2);font-size:13px}}
table.t{{border-collapse:collapse;font-size:12px;font-variant-numeric:tabular-nums;display:block;overflow-x:auto}}
table.t td,table.t th{{padding:3px 8px;border-bottom:1px solid var(--line);text-align:right}}
.muted{{color:var(--muted)}} ul.qc{{columns:4;font-size:13px}}
</style></head><body><main><h1>Tractoflow / sf-tractomics — reproducibility report</h1>
<p class="muted">Blue: Tractoflow (scilpy 2.2.2 port). Orange: sf-tractomics. Same image
(scilus/scilus:2.2.2), reproducible profiles. Subjects coded S01…</p>{body}</main></body></html>"""
    with open(os.path.join(out, "index.html"), "w") as f:
        f.write(page)
    print(os.path.join(out, "index.html"))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawTextHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)
    q = sp.add_parser("qc")
    q.add_argument("--tractoflow", required=True)
    q.add_argument("--sftractomics", required=True)
    q.add_argument("--out_dir", required=True)
    q.add_argument("--sid_regex", default=r"^(sub-[^_]+)_(ses-[^_]+)")
    q.add_argument("--nproc", type=int, default=1)
    g = sp.add_parser("plots")
    g.add_argument("--analysis_dir", required=True)
    args = p.parse_args()
    cmd_qc(args) if args.cmd == "qc" else cmd_plots(args)


if __name__ == "__main__":
    main()
