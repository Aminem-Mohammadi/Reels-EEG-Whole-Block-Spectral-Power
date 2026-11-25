#!/usr/bin/env python3
"""
F01_publication_figures.py  --  intervention EEG, publication figures

Draws the intervention EEG figures with MNE + matplotlib from the MATLAB outputs.
No statistics are re-estimated: every model value (EMM, estimate, SE, p, FDR, Holm)
is read from the MATLAB outputs. Only descriptive quantities (means and 95% CIs of
spectra, paired differences) and t = estimate / SE (used for colouring) are computed.

Contrast maps show t = estimate / SE with robust colour limits (98th percentile).
Markers: filled = omnibus BH-FDR across electrodes AND Holm across the 3 contrasts;
open = p < 0.05 uncorrected only.

RUN ORDER
  MATLAB: X03, X04, X07, X09 (both RUN_MODEs), X05 (both), X25 (both),
          X12 (both), then F00_export_for_python.m
  Python: python3 F01_publication_figures.py  [--root <OUT_DIR>]

REQUIREMENTS  python >= 3.9, mne >= 1.4, matplotlib, numpy, pandas, scipy, openpyxl

OUTPUT  <OUT_DIR>/31_publication_figures/
          main/          Fig1..Fig6 (PDF + PNG)
          topo_panels/   one 2x3 panel per band x metric x RUN_MODE
          check_montage.png   electrode positions (orientation check)
          figure_log.txt
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# ---------------------------------------------------------------------------
# 0. SETTINGS (edit here)
# ---------------------------------------------------------------------------
DEFAULT_ROOT = "results"

RUN_MODES = {"full_sample": "Full sample (n = 30)",
             "sensitivity_excl_noisy": "Sensitivity (excl. 7, 23; n = 28)"}
CONDITIONS = ["E", "N", "R"]
COND_NAMES = {"E": "Reels", "N": "Neutral", "R": "Rest"}
COND_COLORS = {"E": "#D55E00", "N": "#0072B2", "R": "#6E6E6E"}   # Okabe-Ito, colour-blind safe
CONTRASTS = [("E - N", "EminusN"), ("E - R", "EminusR"), ("N - R", "NminusR")]

BANDS = ["Delta", "Theta", "Alpha", "Beta", "Beta_low_exploratory", "Beta_high_exploratory",
         "Gamma_low_exploratory", "Gamma_high1_exploratory", "Gamma_high2_exploratory"]
BAND_LABEL = {"Delta": "Delta 1-4 Hz", "Theta": "Theta 4-8 Hz", "Alpha": "Alpha 8-13 Hz",
              "Beta": "Beta 13-30 Hz", "Beta_low_exploratory": "Beta low 13-20 Hz*",
              "Beta_high_exploratory": "Beta high 20-30 Hz*",
              "Gamma_low_exploratory": "Gamma 30-45 Hz*", "Gamma_high1_exploratory": "Gamma 55-80 Hz*",
              "Gamma_high2_exploratory": "Gamma 80-100 Hz*"}
EXPLORATORY = {b for b in BANDS if b.endswith("_exploratory")}
METRICS = {"LogAbsolutePower": ("Absolute power", "ln(µV²)"),
           "RelativePower": ("Relative power", "fraction of 1-45 Hz")}

ROBUST_PCT = 98            # colour limits = this percentile of the plotted values
T_MIN_LIMIT = 3.0          # t-maps never use a narrower scale than +/- this
ALPHA = 0.05
FIG_W_CM = 17.8            # two-column journal width
FONT = "DejaVu Sans"       # always available with matplotlib; change to "Arial" if installed
CMAP_SEQ = "viridis"
CMAP_DIV = "RdBu_r"

POSTERIOR_ROIS = ["Parietal", "Occipital"]        # ROI v2 names used for spectra / paired plot
KEY_BAND, KEY_METRIC = "Theta", "RelativePower"   # band and metric of the key E - N map (Fig1)

plt.rcParams.update({
    "font.family": FONT, "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
    "ps.fonttype": 42, "savefig.dpi": 300, "figure.dpi": 100})

CM = 1 / 2.54
LOG = []


def log(msg: str) -> None:
    print(msg)
    LOG.append(msg)


# ---------------------------------------------------------------------------
# 1. INPUTS
# ---------------------------------------------------------------------------
class Paths:
    def __init__(self, root: Path):
        self.root = root
        self.inputs = root / "30_python_inputs"
        self.out = root / "31_publication_figures"
        self.main = self.out / "main"
        self.panels = self.out / "topo_panels"
        for d in (self.out, self.main, self.panels):
            d.mkdir(parents=True, exist_ok=True)

    def elec(self, mode, name):
        return self.root / "09_electrode_stats" / mode / name

    def roi2(self, mode):
        return self.root / "25_roi_v2_stats" / mode / "roi_v2_all_contrasts.csv"

    def roi2_def(self):
        return self.root / "24_roi_v2_definitions" / "roi_v2_definitions.csv"

    def global_contrasts(self, mode):
        return self.root / "05_primary_stats" / mode / "phase4_all_contrasts_ungated.xlsx"

    def ratio_pairwise(self, mode):
        return self.root / "12_band_ratio_indices" / mode / "phase12_band_ratio_electrode_pairwise.csv"


def need(path: Path, hint: str) -> Path:
    if not path.exists():
        sys.exit(f"ERROR: missing input {path}\n  -> {hint}")
    return path


def read_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    for c in df.columns:
        if df[c].dtype == object:
            df[c] = df[c].astype(str).str.strip()
    return df


def uc(s):
    return s.astype(str).str.upper().str.strip()


# ---------------------------------------------------------------------------
# 2. ELECTRODE POSITIONS
#    EEGLAB / X07 convention: +X nose, +Y left ear, +Z up.  We project the 3D
#    positions azimuthally (distance from the vertex = colatitude), the same
#    projection EEGLAB topoplot uses: an electrode on the equator (Fpz, T7,
#    Oz ...) lands on the head circle, lower ones (P9, P10) slightly outside.
# ---------------------------------------------------------------------------
HEAD_R = 0.095


def load_positions(paths: Paths):
    loc = read_csv(need(paths.inputs / "channel_locations.csv", "run F00_export_for_python.m in MATLAB"))
    xyz = loc[["X", "Y", "Z"]].to_numpy(float)
    xyz = xyz / np.linalg.norm(xyz, axis=1, keepdims=True)
    colat = np.arccos(np.clip(xyz[:, 2], -1, 1))            # 0 at vertex
    az = np.arctan2(xyz[:, 1], xyz[:, 0])                     # 0 = nose, + toward left ear
    r = colat / (np.pi / 2)
    r = r / max(r.max() / 0.92, 1.0) * HEAD_R   # equator on the circle unless electrodes lie below it,
                                                # then shrink so the lowest (P9/P10) sit just inside
    # plot coords: +x = subject's right, +y = nose
    pos = np.column_stack([-r * np.sin(az), r * np.cos(az)])
    labels = loc["Label"].astype(str).str.strip().tolist()
    return labels, pos


def check_montage(labels, pos, paths: Paths):
    fig, ax = plt.subplots(figsize=(9 * CM, 9 * CM))
    _draw_head(ax)
    ax.scatter(pos[:, 0], pos[:, 1], s=6, c="k")
    for (x, y), lab in zip(pos, labels):
        ax.text(x, y + 0.004, lab, ha="center", va="bottom", fontsize=5)
    ax.set_title("Montage check: nose up, C3 on the LEFT, O1/O2 at the back")
    ax.set_aspect("equal"); ax.axis("off")
    lim = 1.25 * HEAD_R
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    f = paths.out / "check_montage.png"
    fig.savefig(f, bbox_inches="tight"); plt.close(fig)
    i = {l.upper(): k for k, l in enumerate(labels)}
    ok = all(k in i for k in ("C3", "C4", "FZ", "OZ")) and pos[i["C3"], 0] < 0 < pos[i["C4"], 0] \
        and pos[i["FZ"], 1] > 0 > pos[i["OZ"], 1]
    log(f"Montage check written to {f}  ->  orientation {'OK' if ok else 'SUSPICIOUS -- check axes in F00/X07'}")


def _draw_head(ax, lw=0.6):
    t = np.linspace(0, 2 * np.pi, 200)
    ax.plot(HEAD_R * np.cos(t), HEAD_R * np.sin(t), "k", lw=lw)
    ax.plot([-0.09 * HEAD_R, 0, 0.09 * HEAD_R], [0.995 * HEAD_R, 1.09 * HEAD_R, 0.995 * HEAD_R], "k", lw=lw)
    for s in (-1, 1):
        ax.plot(s * HEAD_R * np.array([0.995, 1.04, 1.06, 1.04, 0.995]),
                HEAD_R * np.array([0.12, 0.10, 0.0, -0.12, -0.15]), "k", lw=lw)


# ---------------------------------------------------------------------------
# 3. TOPOMAP PRIMITIVE
# ---------------------------------------------------------------------------
def topo(ax, values, pos, vlim, cmap, strict=None, loose=None, outline=None, contours=3):
    """One head map. strict/loose: boolean masks for filled/open markers.
    outline: list of index arrays, each drawn as a dashed ROI outline."""
    import mne
    values = np.asarray(values, float)
    good = np.isfinite(values)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        im, _ = mne.viz.plot_topomap(values[good], pos[good], axes=ax, show=False, cmap=cmap,
                                     vlim=vlim, contours=contours, sensors=False, extrapolate="head",
                                     sphere=(0.0, 0.0, 0.0, HEAD_R), outlines="head",
                                     image_interp="cubic", res=128)
    ax.scatter(pos[:, 0], pos[:, 1], s=1.2, c="0.25", lw=0, zorder=5)
    if loose is not None and np.any(loose):
        ax.scatter(pos[loose, 0], pos[loose, 1], s=9, facecolors="none", edgecolors="k", lw=0.5, zorder=6)
    if strict is not None and np.any(strict):
        ax.scatter(pos[strict, 0], pos[strict, 1], s=9, c="k", lw=0, zorder=7)
    for idx in outline or []:
        _roi_outline(ax, pos[idx])
    return im


ROI_COLOR = "#009E73"


def _roi_outline(ax, pts, pad=0.011):
    """Rounded convex hull around the ROI electrodes (works for collinear rows too)."""
    from scipy.spatial import ConvexHull
    t = np.linspace(0, 2 * np.pi, 24, endpoint=False)
    cloud = (pts[:, None, :] + pad * np.stack([np.cos(t), np.sin(t)], 1)[None]).reshape(-1, 2)
    h = cloud[ConvexHull(cloud).vertices]
    h = np.vstack([h, h[:1]])
    ax.plot(h[:, 0], h[:, 1], ls="-", lw=1.0, c=ROI_COLOR, zorder=8)


def robust_limits(arrays, symmetric):
    v = np.concatenate([np.asarray(a, float).ravel() for a in arrays])
    v = v[np.isfinite(v)]
    if v.size == 0:
        return (-1, 1)
    if symmetric:
        m = np.percentile(np.abs(v), ROBUST_PCT)
        return (-m, m) if m > 0 else (-1, 1)
    lo, hi = np.percentile(v, [100 - ROBUST_PCT, ROBUST_PCT])
    return (lo, hi) if hi > lo else (lo - 1, hi + 1)


def add_cbar(fig, im, cax, label, extend="neither", orientation="horizontal"):
    cb = fig.colorbar(im, cax=cax, orientation=orientation, extend=extend)
    cb.set_label(label, labelpad=1)
    cb.ax.tick_params(length=2, pad=1)
    cb.outline.set_linewidth(0.4)
    return cb


def _cax(fig, spec, frac=0.7):
    """Colour-bar axis occupying the middle `frac` of a grid cell (vertical bars)."""
    ax = fig.add_subplot(spec)
    b = ax.get_position()
    ax.set_position([b.x0, b.y0 + b.height * (1 - frac) / 2, b.width, b.height * frac])
    return ax


def save(fig, folder: Path, stem: str):
    fig.savefig(folder / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(folder / f"{stem}.png", bbox_inches="tight", dpi=300)
    plt.close(fig)
    log(f"  wrote {folder.name}/{stem}.pdf/.png")


# ---------------------------------------------------------------------------
# 4. ELECTRODE-LEVEL STATISTICS (from X09)
# ---------------------------------------------------------------------------
class ElecStats:
    def __init__(self, paths: Paths, mode: str, labels):
        self.labels_uc = [l.upper() for l in labels]
        est = read_csv(need(paths.elec(mode, "electrode_model_estimates.csv"),
                            f"run X09 with RUN_MODE='{mode}'"))
        omn = read_csv(need(paths.elec(mode, "electrode_omnibus_results.csv"),
                            f"run X09 with RUN_MODE='{mode}'"))
        pw_path = paths.elec(mode, "electrode_pairwise_results.csv")
        pw = read_csv(pw_path) if pw_path.exists() else pd.DataFrame(
            columns=["Electrode", "Metric", "Band", "Contrast", "pHolm"])
        for d in (est, omn, pw):
            d["ElectrodeUC"] = uc(d["Electrode"])
        self.est, self.omn, self.pw = est, omn, pw
        missing = set(self.labels_uc) - set(est["ElectrodeUC"])
        if missing:
            sys.exit(f"ERROR: electrodes in channel_locations.csv but not in X09 estimates ({mode}): {sorted(missing)}")

    def _vec(self, df, band, metric, col):
        d = df[(df["Band"] == band) & (df["Metric"] == metric)].set_index("ElectrodeUC")
        return d.reindex(self.labels_uc)[col].to_numpy(float)

    def emm(self, band, metric, cond):
        return self._vec(self.est, band, metric, f"EMM_{cond}")

    def tmap(self, band, metric, ctag):
        return self._vec(self.est, band, metric, f"Est_{ctag}") / self._vec(self.est, band, metric, f"SE_{ctag}")

    def markers(self, band, metric, clabel, ctag):
        p_unc = self._vec(self.est, band, metric, f"p_{ctag}_uncorrected")
        fdr = self._vec(self.omn.assign(S=self.omn["SignificantFDR"].astype(str).str.lower()
                                        .isin(["1", "true"]).astype(float)), band, metric, "S") > 0
        d = self.pw[(self.pw["Band"] == band) & (self.pw["Metric"] == metric) & (self.pw["Contrast"] == clabel)]
        holm = d.set_index("ElectrodeUC").reindex(self.labels_uc)["pHolm"].to_numpy(float)
        strict = fdr & (holm < ALPHA)
        loose = (p_unc < ALPHA) & ~strict
        return strict, loose


# ---------------------------------------------------------------------------
# 5. FIGURES
# ---------------------------------------------------------------------------
def fig_topo_panel(st: ElecStats, pos, band, metric, mode, paths: Paths):
    """2 x 3 panel: EMMs (shared sequential bar) / t-maps (shared diverging bar)."""
    mname, unit = METRICS[metric]
    emms = [st.emm(band, metric, c) for c in CONDITIONS]
    ts = [st.tmap(band, metric, tag) for _, tag in CONTRASTS]
    if all(np.all(~np.isfinite(e)) for e in emms):
        log(f"  skip {band}/{metric}/{mode}: no estimates")
        return
    vlim_e = robust_limits(emms, symmetric=False)
    m = max(robust_limits(ts, symmetric=True)[1], T_MIN_LIMIT)
    vlim_t = (-m, m)

    fig = plt.figure(figsize=(FIG_W_CM * CM, 11.5 * CM))
    gs = GridSpec(2, 4, figure=fig, width_ratios=[1, 1, 1, 0.045], hspace=0.22, wspace=0.06)
    for k, c in enumerate(CONDITIONS):
        ax = fig.add_subplot(gs[0, k])
        im_e = topo(ax, emms[k], pos, vlim_e, CMAP_SEQ, contours=0)
        ax.set_title(f"{COND_NAMES[c]} ({c})", color=COND_COLORS[c], fontweight="bold")
    add_cbar(fig, im_e, _cax(fig, gs[0, 3]), f"Model-estimated mean\n{unit}", extend="both", orientation="vertical")
    for k, (clabel, tag) in enumerate(CONTRASTS):
        ax = fig.add_subplot(gs[1, k])
        strict, loose = st.markers(band, metric, clabel, tag)
        im_t = topo(ax, ts[k], pos, vlim_t, CMAP_DIV, strict, loose)
        ax.set_title(clabel.replace(" - ", " \u2212 "))
    add_cbar(fig, im_t, _cax(fig, gs[1, 3]), "t (estimate / SE)", extend="both", orientation="vertical")
    suffix = "  [exploratory]" if band in EXPLORATORY else ""
    fig.suptitle(f"{BAND_LABEL[band].rstrip('*')} — {mname}{suffix}\n{RUN_MODES[mode]}", y=0.995)
    fig.legend(handles=_marker_legend(), loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.0))
    save(fig, paths.panels, f"topo_{band}_{metric}_{mode}")


def _marker_legend():
    return [Line2D([], [], marker="o", ls="", mfc="k", mec="k", ms=3.5, label="FDR (electrodes) + Holm (3 contrasts)"),
            Line2D([], [], marker="o", ls="", mfc="none", mec="k", ms=3.5, label="p < 0.05 uncorrected only")]


def roi_indices(paths: Paths, labels):
    f = paths.roi2_def()
    if not f.exists():
        log("  (ROI v2 definitions not found -- ROI outlines and ROI-based panels skipped)")
        return {}
    d = read_csv(f)
    idx = {l.upper(): k for k, l in enumerate(labels)}
    out = {}
    for roi, g in d.groupby("ROI_v2"):
        if roi and roi.lower() not in ("nan", "none", ""):
            out[roi] = np.array([idx[l.upper()] for l in g["Label"] if l.upper() in idx])
    return out


def fig1_key_EN(stats, pos, rois, paths: Paths):
    """Fig 1: E-N t-maps for the key band/metric (both RUN_MODEs) + log Alpha/Theta, ROI outline."""
    panels = []
    for mode in RUN_MODES:
        if mode in stats:
            st = stats[mode]
            s, l = st.markers(KEY_BAND, KEY_METRIC, "E - N", "EminusN")
            panels.append((f"Relative theta, E − N\n{RUN_MODES[mode]}", st.tmap(KEY_BAND, KEY_METRIC, "EminusN"), s, l))
    ratio = _ratio_tmap(paths, "full_sample", stats["full_sample"].labels_uc)
    if ratio is not None:
        panels.append(("ln(Alpha/Theta), E − N\nFull sample (n = 30)",) + ratio)
    outline = [rois[r] for r in POSTERIOR_ROIS if r in rois]
    m = max(robust_limits([p[1] for p in panels], True)[1], T_MIN_LIMIT)
    fig = plt.figure(figsize=(FIG_W_CM * CM, 7.2 * CM))
    gs = GridSpec(2, len(panels), figure=fig, height_ratios=[1, 0.06], hspace=0.25, wspace=0.08)
    for k, (title, t, s, l) in enumerate(panels):
        ax = fig.add_subplot(gs[0, k])
        im = topo(ax, t, pos, (-m, m), CMAP_DIV, s, l, outline)
        ax.set_title(title)
    add_cbar(fig, im, fig.add_subplot(gs[1, len(panels) // 2] if len(panels) > 1 else gs[1, 0]),
             "t (estimate / SE)", extend="both")
    handles = _marker_legend() + [Line2D([], [], ls="-", c=ROI_COLOR, lw=1.0, label="ROI v2 parietal / occipital")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.06))
    save(fig, paths.main, "Fig1_key_EminusN_posterior_theta")


def _ratio_tmap(paths, mode, labels_uc):
    f = paths.ratio_pairwise(mode)
    if not f.exists():
        log("  (X12 electrode pairwise file not found -- Alpha/Theta panel skipped)")
        return None
    d = read_csv(f)
    unit = "Channel" if "Channel" in d.columns else d.columns[2]
    d = d[(d["Variant"] == "LogMeanBased") & (d["RatioType"] == "AlphaTheta") & (d["Contrast"] == "E - N")]
    if d.empty:
        return None
    d = d.assign(UC=uc(d[unit])).set_index("UC").reindex(labels_uc)
    t = (d["Estimate"] / d["SE"]).to_numpy(float)
    strict = (d["pHolm"] < ALPHA).to_numpy()          # rows exist only where the electrode passed FDR
    loose = ((d["pValue"] < ALPHA) & ~(d["pHolm"] < ALPHA)).to_numpy()
    return t, strict, loose


def _paired_values(bp_e, bp_g, rois):
    """Per-subject values for the E&N subjects: global relative alpha and posterior relative theta."""
    out = {}
    g = bp_g[(bp_g["Band"] == "Alpha")]
    out["Relative alpha, whole scalp"] = g.pivot_table(index="Subject", columns="Condition", values="RelativePower")
    chans = set()
    for r in POSTERIOR_ROIS:
        chans |= set(rois.get(r, []))
    if chans:
        e = bp_e[(bp_e["Band"] == "Theta") & bp_e["ChanIdx"].isin(chans)]
        out["Relative theta, parietal + occipital"] = e.groupby(["Subject", "Condition"])["RelativePower"].mean().unstack()
    return {k: v.dropna(subset=["E", "N"]) for k, v in out.items() if {"E", "N"} <= set(v.columns)}


def fig2_paired(bp_e, bp_g, rois, mode_excl, paths: Paths):
    """Fig 2: subjects who saw both Reels and Neutral (descriptive, uncorrected)."""
    data = _paired_values(bp_e, bp_g, rois)
    if not data:
        log("  skip Fig2: no E&N subjects found"); return
    fig, axes = plt.subplots(1, len(data), figsize=(FIG_W_CM * 0.6 * CM, 6.0 * CM))
    axes = np.atleast_1d(axes)
    for ax, (title, d) in zip(axes, data.items()):
        diff = d["E"] - d["N"]
        for s, row in d.iterrows():
            excluded = int(s) in mode_excl
            ax.plot([0, 1], [row["E"], row["N"]], c="0.55", lw=0.8, ls=":" if excluded else "-")
            ax.scatter([0, 1], [row["E"], row["N"]], s=8, c=[COND_COLORS["E"], COND_COLORS["N"]], zorder=3)
        for x, c in ((0, "E"), (1, "N")):
            mu, ci = d[c].mean(), 1.96 * d[c].std(ddof=1) / np.sqrt(len(d))
            ax.errorbar(x + (-0.12 if x == 0 else 0.12), mu, yerr=ci, fmt="s", ms=4, c="k", capsize=2, lw=0.8)
        k = int((diff < 0).sum()) if "alpha" in title else int((diff > 0).sum())
        word = "lower" if "alpha" in title else "higher"
        ax.set_title(f"{title}\n{k}/{len(d)} subjects {word} in Reels")
        ax.set_xticks([0, 1], [COND_NAMES["E"], COND_NAMES["N"]]); ax.set_xlim(-0.35, 1.35)
        ax.set_ylabel("Relative power")
    fig.text(0.5, -0.04, "Lines = individual subjects with both conditions; squares = mean ± 95% CI. "
             "Descriptive, not corrected for multiplicity.", ha="center", fontsize=6)
    fig.tight_layout()
    save(fig, paths.main, "Fig2_paired_subjects_EN")


def fig3_spectra(paths: Paths, labels, rois):
    """Fig 3: posterior spectra per condition and the within-subject E-N difference."""
    from scipy.io import loadmat
    f = paths.inputs / "psd_recordings.mat"
    if not f.exists():
        log("  skip Fig3: psd_recordings.mat missing (run F00)"); return
    m = loadmat(f, squeeze_me=True)
    psd, fr = np.asarray(m["PSD"], float), np.asarray(m["FreqHz"], float).ravel()
    chl = [str(x).strip().upper() for x in np.atleast_1d(m["ChannelLabels"])]
    subj = np.atleast_1d(m["RecSubject"]).astype(int)
    cond = np.array([str(x).strip() for x in np.atleast_1d(m["RecCondition"])])
    want = set()
    for r in POSTERIOR_ROIS:
        want |= {labels[i].upper() for i in rois.get(r, [])}
    ci = [k for k, l in enumerate(chl) if l in want] or list(range(len(chl)))
    lp = np.log10(psd[:, ci, :].mean(axis=1))                 # recording x freq
    keep = fr <= 45
    fr, lp = fr[keep], lp[:, keep]

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(FIG_W_CM * 0.6 * CM, 10 * CM), sharex=True,
                                 gridspec_kw={"height_ratios": [1.4, 1], "hspace": 0.32})
    for b, (lo, hi), col in (("theta", (4, 8), "#F0E442"), ("alpha", (8, 13), "#56B4E9")):
        for a in (a1, a2):
            a.axvspan(lo, hi, color=col, alpha=0.18, lw=0)
        a1.text((lo + hi) / 2, 1.0, b, transform=a1.get_xaxis_transform(), ha="center", va="bottom", fontsize=6)
    for c in CONDITIONS:
        x = lp[cond == c]
        if len(x) < 2:
            continue
        mu, se = x.mean(0), x.std(0, ddof=1) / np.sqrt(len(x))
        a1.plot(fr, mu, c=COND_COLORS[c], lw=1, label=f"{COND_NAMES[c]} (n = {len(x)})")
        a1.fill_between(fr, mu - 1.96 * se, mu + 1.96 * se, color=COND_COLORS[c], alpha=0.2, lw=0)
    a1.set_ylabel("log₁₀ power (µV²/Hz)")
    a1.legend(frameon=False, loc="upper right")
    a1.set_title("Parietal + occipital ROI (recording means, 95% CI)", pad=10)
    # within-subject difference, E&N subjects only
    d = []
    for s in np.unique(subj):
        ie, inn = np.where((subj == s) & (cond == "E"))[0], np.where((subj == s) & (cond == "N"))[0]
        if len(ie) == 1 and len(inn) == 1:
            d.append(lp[ie[0]] - lp[inn[0]])
    if len(d) >= 3:
        d = np.array(d)
        mu, se = d.mean(0), d.std(0, ddof=1) / np.sqrt(len(d))
        from scipy.stats import t as tdist
        q = tdist.ppf(0.975, len(d) - 1)
        a2.plot(fr, mu, c="k", lw=1)
        a2.fill_between(fr, mu - q * se, mu + q * se, color="k", alpha=0.18, lw=0)
        a2.axhline(0, c="0.5", lw=0.6)
        a2.set_title(f"Within-subject Reels − Neutral (n = {len(d)} subjects with both; 95% CI)")
    a2.set_ylabel("Δ log₁₀ power")
    a2.set_xlabel("Frequency (Hz)")
    a2.set_xlim(fr.min(), fr.max())
    save(fig, paths.main, "Fig3_posterior_spectra")


def _forest(ax, rows, title, xlabel):
    """rows: list of (label, {mode: (est, lo, hi, p)}, exploratory)."""
    y = np.arange(len(rows))[::-1]
    offs = {m: o for m, o in zip(RUN_MODES, (0.15, -0.15))}
    mk = {"full_sample": "o", "sensitivity_excl_noisy": "s"}
    for yi, (lab, vals, expl) in zip(y, rows):
        for mode, (e, lo, hi, p) in vals.items():
            col = "0.55" if expl else "k"
            fill = col if p < ALPHA else "white"
            ax.plot([lo, hi], [yi + offs[mode]] * 2, c=col, lw=0.8)
            ax.plot(e, yi + offs[mode], marker=mk[mode], ms=3.5, mfc=fill, mec=col, ls="")
    ax.axvline(0, c="0.6", lw=0.6)
    ax.set_yticks(y, [r[0] for r in rows])
    ax.set_title(title); ax.set_xlabel(xlabel)


def _forest_legend(fig):
    h = [Line2D([], [], marker="o", ls="", mfc="k", mec="k", ms=3.5, label=RUN_MODES["full_sample"]),
         Line2D([], [], marker="s", ls="", mfc="k", mec="k", ms=3.5, label=RUN_MODES["sensitivity_excl_noisy"]),
         Line2D([], [], marker="o", ls="", mfc="white", mec="k", ms=3.5, label="open = p ≥ 0.05 (model, uncorrected)"),
         Patch(fc="0.55", label="* exploratory band")]
    fig.legend(handles=h, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.08))


def fig4_forest_global(paths: Paths):
    """Fig 4: whole-scalp E-N estimates with 95% CI (X05), both RUN_MODEs."""
    tabs = {}
    for mode in RUN_MODES:
        f = paths.global_contrasts(mode)
        if f.exists():
            tabs[mode] = pd.read_excel(f)
    if not tabs:
        log("  skip Fig4: phase4_all_contrasts_ungated.xlsx not found"); return
    fig, axes = plt.subplots(1, 2, figsize=(FIG_W_CM * CM, 7 * CM), sharey=True)
    for ax, metric in zip(axes, METRICS):
        rows = []
        for b in BANDS:
            vals = {}
            for mode, t in tabs.items():
                r = t[(t["Metric"] == metric) & (t["Band"] == b) & (t["Contrast"] == "E - N")]
                if len(r):
                    r = r.iloc[0]
                    vals[mode] = (r["Estimate"], r["CI_Low"], r["CI_High"], r["pValue"])
            rows.append((BAND_LABEL[b], vals, b in EXPLORATORY))
        _forest(ax, rows, f"{METRICS[metric][0]}: Reels − Neutral", f"Estimate ({METRICS[metric][1]})")
    _forest_legend(fig); fig.tight_layout()
    save(fig, paths.main, "Fig4_forest_global_EminusN")


def fig5_forest_roi(paths: Paths):
    """Fig 5: ROI v2 E-N relative-power estimates (X25), both RUN_MODEs."""
    tabs = {m: read_csv(paths.roi2(m)) for m in RUN_MODES if paths.roi2(m).exists()}
    if not tabs:
        log("  skip Fig5: roi_v2_all_contrasts.csv not found"); return
    rois = list(dict.fromkeys(next(iter(tabs.values()))["ROI"]))
    bands = ["Delta", "Theta", "Alpha", "Beta", "Beta_low_exploratory", "Beta_high_exploratory", "Gamma_low_exploratory"]
    fig, axes = plt.subplots(1, len(rois), figsize=(FIG_W_CM * CM, 6.5 * CM), sharey=True)
    for ax, roi in zip(np.atleast_1d(axes), rois):
        rows = []
        for b in bands:
            vals = {}
            for mode, t in tabs.items():
                r = t[(t["ROI"] == roi) & (t["Metric"] == "RelativePower") & (t["Band"] == b) & (t["Contrast"] == "E - N")]
                if len(r):
                    r = r.iloc[0]
                    vals[mode] = (r["Estimate"], r["CI_Low"], r["CI_High"], r["pValue"])
            rows.append((BAND_LABEL[b], vals, b in EXPLORATORY))
        _forest(ax, rows, roi, "Δ relative power")
    _forest_legend(fig); fig.tight_layout()
    save(fig, paths.main, "Fig5_forest_ROIv2_EminusN")


def fig6_sensitivity(stats, pos, paths: Paths):
    """Fig 6: E-N t-maps, FullSample vs Sensitivity side by side, same scale."""
    sel = [("Theta", "RelativePower"), ("Alpha", "RelativePower"), ("Gamma_low_exploratory", "RelativePower")]
    modes = [m for m in RUN_MODES if m in stats]
    if len(modes) < 2:
        log("  skip Fig6: needs X09 output for both RUN_MODEs"); return
    fig = plt.figure(figsize=(9.5 * CM, 12.5 * CM))
    gs = GridSpec(len(sel) + 1, 3, figure=fig, width_ratios=[0.28, 1, 1],
                  height_ratios=[1] * len(sel) + [0.05], hspace=0.18, wspace=0.04)
    allt = [stats[m].tmap(b, me, "EminusN") for b, me in sel for m in modes]
    lim = max(robust_limits(allt, True)[1], T_MIN_LIMIT)
    for i, (b, me) in enumerate(sel):
        lab = fig.add_subplot(gs[i, 0]); lab.axis("off")
        lab.text(0.5, 0.5, BAND_LABEL[b].replace(" ", "\n", 1), rotation=90, ha="center", va="center", fontsize=7)
        for j, mode in enumerate(modes):
            st = stats[mode]
            s_, l_ = st.markers(b, me, "E - N", "EminusN")
            ax = fig.add_subplot(gs[i, j + 1])
            im = topo(ax, st.tmap(b, me, "EminusN"), pos, (-lim, lim), CMAP_DIV, s_, l_)
            if i == 0:
                ax.set_title(RUN_MODES[mode].replace(" (", "\n("))
    cax = fig.add_subplot(gs[-1, 1:])
    b0 = cax.get_position(); cax.set_position([b0.x0 + b0.width * 0.15, b0.y0, b0.width * 0.7, b0.height])
    add_cbar(fig, im, cax, "t, Reels \u2212 Neutral (relative power)", extend="both")
    fig.legend(handles=_marker_legend(), loc="lower center", ncol=1, frameon=False, bbox_to_anchor=(0.55, -0.06))
    save(fig, paths.main, "Fig6_EminusN_FullSample_vs_Sensitivity")


# ---------------------------------------------------------------------------
# 6. MAIN
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--root", default=DEFAULT_ROOT, help="ANALYSIS_ROOT of the intervention pipeline")
    ap.add_argument("--bands", default="all", help="comma list for the topo panels, or 'all'")
    ap.add_argument("--exclude", default="7,23", help="subjects of the sensitivity run (dotted lines in Fig 2)")
    a = ap.parse_args()
    import mne
    log(f"F01 | mne {mne.__version__} | root = {a.root}")
    paths = Paths(Path(a.root))
    labels, pos = load_positions(paths)
    check_montage(labels, pos, paths)
    rois = roi_indices(paths, labels)

    stats = {}
    for mode in RUN_MODES:
        if paths.elec(mode, "electrode_model_estimates.csv").exists():
            stats[mode] = ElecStats(paths, mode, labels)
        else:
            log(f"  (no X09 output for {mode} -- that RUN_MODE is skipped)")
    if "full_sample" not in stats:
        sys.exit("ERROR: X09 full_sample output is required.")

    log("Topo panels (2 x 3) ...")
    bands = BANDS if a.bands == "all" else [b.strip() for b in a.bands.split(",")]
    for mode, st in stats.items():
        for b in bands:
            for me in METRICS:
                fig_topo_panel(st, pos, b, me, mode, paths)

    log("Main figures ...")
    fig1_key_EN(stats, pos, rois, paths)
    bp_e = read_csv(need(paths.inputs / "bandpower_electrode.csv", "run F00_export_for_python.m"))
    bp_g = read_csv(need(paths.inputs / "bandpower_global.csv", "run F00_export_for_python.m"))
    for d in (bp_e, bp_g):
        d["Subject"] = d["Subject"].astype(str).str.replace(r"\D", "", regex=True).astype(int)
    idx = {l.upper(): k for k, l in enumerate(labels)}
    bp_e["ChanIdx"] = uc(bp_e["Channel"]).map(idx)
    excl = {int(x) for x in a.exclude.split(",") if x.strip()}
    fig2_paired(bp_e, bp_g, rois, excl, paths)
    fig3_spectra(paths, labels, rois)
    fig4_forest_global(paths)
    fig5_forest_roi(paths)
    fig6_sensitivity(stats, pos, paths)
    (paths.out / "figure_log.txt").write_text("\n".join(LOG) + "\n")
    log(f"Done. Figures in {paths.out}")


if __name__ == "__main__":
    main()
