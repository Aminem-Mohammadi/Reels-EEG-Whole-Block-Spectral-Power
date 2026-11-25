#!/usr/bin/env python3
"""
F01_publication_figures_cognitive.py  --  cognitive_task pipeline, publication figures

TASK-SIDE PORT OF F01_publication_figures.py (Intervention pipeline). Produces
the same family of figures (topo panels, forest plots, paired-subject plot,
posterior spectra, sensitivity comparison), run as THREE INDEPENDENT, PARALLEL
passes -- one per Task (Full / Nback / SART) -- never pooled, with each Task's
figures written to its own output folder, exactly mirroring how every MATLAB
statistics script in this pipeline (X05/X06/X08/X09/X24/X25) already treats
Task as an outer loop, never as a factor.

NO statistics are re-estimated here. Every model number (EMM, estimate, SE, p,
Holm, FDR) is read from the MATLAB outputs. The only things computed here are
descriptive (means and 95% CIs of spectra, paired differences) and the
t = estimate / SE used for colouring -- identical scope to the Intervention
script.

===========================================================================
WHAT IS DIFFERENT FROM THE INTERVENTION VERSION, AND WHY (full disclosure --
read this before interpreting any figure)
===========================================================================

1. TASK LOOPING. Every figure is produced once per Task (Full/Nback/SART),
   written to <ANALYSIS_ROOT>/<Task>/31_publication_figures/. This is not a
   stylistic choice -- it mirrors the fact that every one of this pipeline's
   own statistics scripts (X05, X08, X09, X09B, X12, X25) fits Task as three
   fully independent copies of the same model, never crossed with Condition
   and never pooled. channel_locations.csv, the ROI v2 definitions (X24),
   and the montage check are Task-invariant (same 62-electrode montage,
   same ROI v2 set for all three Tasks) and are therefore computed once and
   reused across all three Task passes, but every quantity that comes from
   a Task-looped statistics script is kept strictly separate per Task.

2. TWO COEXISTING, DELIBERATELY NOT-HARMONIZED CONTRAST CONVENTIONS. This
   pipeline (unlike Intervention, which uses one convention throughout) has
   two different pairwise-contrast label/sign conventions, carried over
   verbatim from the underlying MATLAB scripts and NEVER silently
   relabeled or sign-flipped here:
     - ELECTRODE-LEVEL (X09 / X09B): E-N / E-R / N-R (tags EminusN, EminusR,
       NminusR) -- this happens to be the SAME convention Intervention uses,
       so the ElecStats class, the topo panels, and Fig6 below are ported
       UNCHANGED (no sign handling needed).
     - GLOBAL-LEVEL (X05 / X05B) and ROI v2 (X24 / X25_cognitive) and the
       band-ratio indices (X12): N-E / R-E / R-N (the OPPOSITE sign on every
       one of the 3 contrasts). Fig4 (global forest) and Fig5 (ROI v2
       forest) below therefore filter on "N - E" (NOT "E - N" as in the
       Intervention script) and are titled using that literal label, so the
       plotted sign always matches the underlying CSV's own Estimate/CI
       columns exactly. This is a real, disclosed consequence of reusing
       each MATLAB script's own locked convention rather than harmonizing
       them for this Python port -- see X05/X09's own script headers, which
       already state the two conventions were never required to match.
   DISPLAY CONVENTION: Fig1/Fig4/Fig5
   DISPLAY the global/ROI-level values under the electrode-level's own E-N/
   E-R/N-R label too, so every figure in this script reads the same
   direction. This is implemented as a pure sign flip of the ALREADY-READ
   Estimate/CI (E-N := -(N-E), CI swapped-and-negated) done only in this
   Python script at render time -- see CONTRAST_DISPLAY_FLIP/flip_sign()
   above. The underlying MATLAB CSVs are read/filtered using their own
   native N-E/R-E/R-N label exactly as before (nothing upstream changed,
   nothing reran); only what gets drawn on the page is negated, and every
   such figure says so in a footnote. Fig6/topo panels use the native X09/
   X09B output directly (already E-N/E-R/N-R, see item above), so they need
   no transform.

3. "FIG1" IS REDESIGNED, NOT PORTED 1:1. Intervention's own Fig1 was a
   bespoke figure built around ONE specific, already-identified finding
   (posterior Theta, E-N, chosen after seeing Intervention's own results --
   see that script's KEY_BAND/KEY_METRIC constants). This Task pipeline does
   not have an equivalent single "key finding" pick; what it does have is
   its own pre-specified ROI v2 confirmatory layer (X24/X25_cognitive): H1
   (Occipital, Alpha, R - E) and H2 (FrontalMidline, Theta, N - E) -- see
   X24_roi_v2_definitions_cognitive.m's own disclosure that these were
   chosen post-hoc but are literature-grounded and deliberately parallel to
   Intervention's own H1/H2 for cross-dataset comparability. Fig1 here is
   therefore a small forest-style summary of H1 and H2's own already-
   computed, correctly-signed Estimate + 95% CI (both RUN_MODEs), read
   directly from X25_cognitive's roi_v2_confirmatory.csv -- NOT an
   electrode-level topomap.
   DISPLAY CONVENTION: H1/H2's own native contrast label (R-E/N-E) is
   sign-flipped for DISPLAY to E-R/E-N, the same transform as Fig4/Fig5 --
   see item 2's update note above. This remains a verbatim-sourced,
   already-computed number (just negated for display); no electrode-level
   map is built here, and no quantity is estimated or re-derived.

4. X12's band-ratio Alpha/Theta panel (Intervention's Fig1 third panel) is
   NOT reproduced here. The only targeted Alpha/Theta analysis in this
   pipeline is X21_targeted_alphatheta_confirmatory_cognitive.m, which that
   script's own extensive DISCLOSURE_TEXT states was chosen because it was
   the closest-to-significant result in an earlier run and explicitly
   warns its p-values must NOT be read as independent confirmatory
   evidence. Promoting that specific result to a polished publication
   figure here would contradict its own disclosed selection-bias caveat,
   so it is deliberately omitted. (X12's own, non-targeted, pre-specified
   ratio results are not currently used by any figure in this script either
   -- Intervention's own Fig1 only used X12 for the one targeted ratio.)

5. Source file/column names: 'InterventionNumber' does not exist in this
   pipeline; 'TaskType' + 'BlockNumber' are used wherever Intervention
   referenced InterventionNumber (bandpower CSVs, psd_recordings.mat). Band
   list is unchanged (same 7 bands, same tier structure) -- verified
   identical PARAMS3.BandDefinitionsHz between pipelines.

6. COND_NAMES / COND_COLORS are REUSED UNCHANGED. X01_audit_cognitive_data.m
   itself documents the Condition letters identically to Intervention: "...
   Instagram Reels [E] / Neutral Images [N] / Rest [R]...", verifying this
   is the SAME E/N/R manipulation, recorded during cognitive-task
   performance (Full/Nback/SART) rather than during a dedicated
   intervention-only session. This is not a guess -- it is read verbatim
   from this pipeline's own script header.

RUN ORDER
  MATLAB: X03, X04, X07, X09 (both RUN_MODEs, all computed inside its own
          Task loop), X09B (both RUN_MODEs; companion to X09), X05 (both
          RUN_MODEs), X05B (both RUN_MODEs; companion to X05), X24, X25
          (both RUN_MODEs), then F00_export_for_python_cognitive.m
  Python: python3 F01_publication_figures_cognitive.py  [--root <ANALYSIS_ROOT>]

REQUIREMENTS  python >= 3.9, mne >= 1.4, matplotlib, numpy, pandas, scipy,
              openpyxl  (pip install mne matplotlib pandas scipy openpyxl)

OUTPUT  <ANALYSIS_ROOT>/<Full|Nback|SART>/31_publication_figures/
          main/          Fig1, Fig2, Fig3, Fig4, Fig5, Fig6 (PDF + PNG)
          topo_panels/   one 2x3 panel per band x metric x RUN_MODE
          check_montage.png   electrode positions -- CHECK THIS FIRST
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
DEFAULT_ROOT = "results"  # pass --root to point at your own ANALYSIS_ROOT

TASKS = {"full": "Full", "nback": "Nback", "sart": "SART"}

# NOTE: the subject numbers named in the "Sensitivity" label below are a
# display label only; EXCLUSION_CONFIG_cognitive.m's EXCLUDE_NOISY_SUBJECTS
# is this pipeline's single source of truth for which subjects are actually
# excluded in sensitivity_excl_noisy mode. Verify this label text against
# that file's current value before publication -- it is not re-derived from
# it automatically.
RUN_MODES = {"full_sample": "Full sample", "sensitivity_excl_noisy": "Sensitivity (excl. 7, 23)"}
CONDITIONS = ["E", "N", "R"]
# Reused UNCHANGED from Intervention -- verified identical manipulation; see
# disclosure item 6 above.
COND_NAMES = {"E": "Reels", "N": "Neutral", "R": "Rest"}
COND_COLORS = {"E": "#D55E00", "N": "#0072B2", "R": "#6E6E6E"}   # Okabe-Ito, colour-blind safe

# Electrode-level (X09/X09B) native tags -- SAME convention as Intervention.
# UNCHANGED/UNTOUCHED by the display-flip below: X09/X09B's own
# output is already E-N/E-R/N-R, so Fig6 and the topo panels (which read
# ONLY from X09/X09B) need no transform at all.
CONTRASTS = [("E - N", "EminusN"), ("E - R", "EminusR"), ("N - R", "NminusR")]
# Global-level (X05/X05B) and ROI v2 (X25_cognitive) NATIVE contrast label as
# actually stored on disk -- OPPOSITE sign convention from the electrode
# level (see disclosure item 2). Used ONLY to filter/read the CSV; see
# CONTRAST_DISPLAY_FLIP below for what is actually shown in Fig1/Fig4/Fig5.
GLOBAL_ROI_CONTRAST_LABEL = "N - E"

# ---------------------------------------------------------------------------
# DISPLAY-ONLY SIGN FLIP --
# same logic as F02/F03's own CONTRAST_DISPLAY_FLIP/flip_sign(). Fig1/Fig4/
# Fig5 below read X05B/X25_cognitive CSVs, which are natively N-E/R-E/R-N
# (GLOBAL_ROI_CONTRAST_LABEL). This pipeline
# displays those figures under the X09/X09B-native E-N/E-R/N-R label too, so
# every figure in this script (and F02/F03) reads the same direction:
# E-N := -(N-E), CI swapped-and-negated. This is a pure sign flip of the
# ALREADY-COMPUTED Estimate, applied only in this Python script; MATLAB
# output is never touched/rerun, p-values/significance are unaffected.
# Fig6/topo panels are NOT touched -- they already read X09/X09B's native
# E-N/E-R/N-R output directly (see CONTRASTS above).
CONTRAST_DISPLAY_FLIP = {"N - E": "E - N", "R - E": "E - R", "R - N": "N - R"}
GLOBAL_ROI_DISPLAY_LABEL = CONTRAST_DISPLAY_FLIP[GLOBAL_ROI_CONTRAST_LABEL]  # "E - N"


def flip_sign(native_label: str, est: float, ci_low: float, ci_high: float):
    """Returns (display_label, display_estimate, display_ci_low, display_ci_high)."""
    disp = CONTRAST_DISPLAY_FLIP.get(native_label, native_label)
    if disp == native_label:
        return disp, est, ci_low, ci_high
    return disp, -est, -ci_high, -ci_low

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

POSTERIOR_ROIS = ["Parietal", "Occipital"]        # ROI v2 names used for spectra (Fig3)
# H1/H2 straight from X24_roi_v2_definitions_cognitive.m -- NOT re-chosen here.
CONFIRMATORY_HYPOTHESES = [
    ("H1", "Occipital", "Alpha", "R - E"),
    ("H2", "FrontalMidline", "Theta", "N - E"),
]

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
    """Task-aware path resolver. channel_locations.csv, psd_recordings.mat,
    bandpower_*.csv and the ROI v2 DEFINITION file (X24) are Task-invariant
    (shared); every statistics output below them is Task-looped on disk, so
    every method that reaches into a Task's own statistics output takes a
    `task` key ('full'/'nback'/'sart')."""

    def __init__(self, root: Path):
        self.root = root
        self.inputs = root / "30_python_inputs"

    def task_label(self, task: str) -> str:
        return TASKS[task]

    def out(self, task: str) -> Path:
        d = self.root / self.task_label(task) / "31_publication_figures"
        d.mkdir(parents=True, exist_ok=True)
        return d

    def main(self, task: str) -> Path:
        d = self.out(task) / "main"
        d.mkdir(parents=True, exist_ok=True)
        return d

    def panels(self, task: str) -> Path:
        d = self.out(task) / "topo_panels"
        d.mkdir(parents=True, exist_ok=True)
        return d

    def elec(self, task: str, mode: str, name: str) -> Path:
        return self.root / self.task_label(task) / "09_electrode_stats" / mode / name

    def roi2(self, task: str, mode: str) -> Path:
        return self.root / self.task_label(task) / "25_roi_v2_stats" / mode / "roi_v2_all_contrasts.csv"

    def roi2_confirmatory(self, task: str, mode: str) -> Path:
        return self.root / self.task_label(task) / "25_roi_v2_stats" / mode / "roi_v2_confirmatory.csv"

    def roi2_def(self) -> Path:
        # Shared, Task-invariant -- X24 is a single, non-Task-looped execution.
        return self.root / "24_roi_v2_definitions" / "roi_v2_definitions.csv"

    def global_contrasts(self, task: str, mode: str) -> Path:
        return self.root / self.task_label(task) / "05_primary_stats" / mode / "phase4_all_contrasts_ungated.xlsx"


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
# 2. ELECTRODE POSITIONS (Task-invariant -- same montage for all three Tasks)
#    EEGLAB / X07 convention: +X nose, +Y left ear, +Z up.  We project the 3D
#    positions azimuthally (distance from the vertex = colatitude), the same
#    projection EEGLAB topoplot uses: an electrode on the equator (Fpz, T7,
#    Oz ...) lands on the head circle, lower ones (P9, P10) slightly outside.
# ---------------------------------------------------------------------------
HEAD_R = 0.095


def load_positions(paths: Paths):
    loc = read_csv(need(paths.inputs / "channel_locations.csv", "run F00_export_for_python_cognitive.m in MATLAB"))
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


def check_montage(labels, pos, out_dir: Path):
    fig, ax = plt.subplots(figsize=(9 * CM, 9 * CM))
    _draw_head(ax)
    ax.scatter(pos[:, 0], pos[:, 1], s=6, c="k")
    for (x, y), lab in zip(pos, labels):
        ax.text(x, y + 0.004, lab, ha="center", va="bottom", fontsize=5)
    ax.set_title("Montage check: nose up, C3 on the LEFT, O1/O2 at the back")
    ax.set_aspect("equal"); ax.axis("off")
    lim = 1.25 * HEAD_R
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    f = out_dir / "check_montage.png"
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
# 3. TOPOMAP PRIMITIVE (unchanged from Intervention)
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
# 4. ELECTRODE-LEVEL STATISTICS (from X09 + X09B), per Task
# ---------------------------------------------------------------------------
class ElecStats:
    def __init__(self, paths: Paths, task: str, mode: str, labels):
        self.labels_uc = [l.upper() for l in labels]
        est = read_csv(need(paths.elec(task, mode, "electrode_model_estimates.csv"),
                            f"run X09B with RUN_MODE='{mode}' for Task={paths.task_label(task)}"))
        omn = read_csv(need(paths.elec(task, mode, "electrode_omnibus_results.csv"),
                            f"run X09 with RUN_MODE='{mode}' for Task={paths.task_label(task)}"))
        pw_path = paths.elec(task, mode, "electrode_pairwise_results.csv")
        pw = read_csv(pw_path) if pw_path.exists() else pd.DataFrame(
            columns=["Electrode", "Metric", "Band", "Contrast", "pHolm"])
        for d in (est, omn, pw):
            d["ElectrodeUC"] = uc(d["Electrode"])
        self.est, self.omn, self.pw = est, omn, pw
        missing = set(self.labels_uc) - set(est["ElectrodeUC"])
        if missing:
            sys.exit(f"ERROR: electrodes in channel_locations.csv but not in X09B estimates "
                      f"(Task={paths.task_label(task)}, {mode}): {sorted(missing)}")

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
def fig_topo_panel(st: ElecStats, pos, band, metric, mode, task_label, out_panels: Path):
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
        ax.set_title(clabel.replace(" - ", " − "))
    add_cbar(fig, im_t, _cax(fig, gs[1, 3]), "t (estimate / SE)", extend="both", orientation="vertical")
    suffix = "  [exploratory]" if band in EXPLORATORY else ""
    fig.suptitle(f"{BAND_LABEL[band].rstrip('*')} — {mname}{suffix}\nTask={task_label} | {RUN_MODES[mode]}", y=0.995)
    fig.legend(handles=_marker_legend(), loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.0))
    save(fig, out_panels, f"topo_{band}_{metric}_{mode}")


def _marker_legend():
    return [Line2D([], [], marker="o", ls="", mfc="k", mec="k", ms=3.5, label="FDR (electrodes) + Holm (3 contrasts)"),
            Line2D([], [], marker="o", ls="", mfc="none", mec="k", ms=3.5, label="p < 0.05 uncorrected only")]


def roi_indices(paths: Paths, labels):
    """Shared, Task-invariant -- same ROI v2 set for all three Tasks (X24)."""
    f = paths.roi2_def()
    if not f.exists():
        log("  (ROI v2 definitions not found -- ROI-based panels skipped)")
        return {}
    d = read_csv(f)
    idx = {l.upper(): k for k, l in enumerate(labels)}
    out = {}
    for roi, g in d.groupby("ROI_v2"):
        if roi and roi.lower() not in ("nan", "none", ""):
            out[roi] = np.array([idx[l.upper()] for l in g["Label"] if l.upper() in idx])
    return out


def fig1_confirmatory(paths: Paths, task: str, task_label: str):
    """Fig 1 (REDESIGNED -- see disclosure item 3 in the module docstring):
    H1 and H2's own already-computed Estimate + 95% CI (both RUN_MODEs), read
    verbatim from X25_cognitive's roi_v2_confirmatory.csv -- no electrode-level
    map. H1/H2's native contrast ('R - E'/'N - E') is sign-flipped for DISPLAY
    only to 'E - R'/'E - N' via flip_sign()
    above, same transform as Fig4/Fig5; no new statistic is computed here."""
    tabs = {}
    for mode in RUN_MODES:
        f = paths.roi2_confirmatory(task, mode)
        if f.exists():
            tabs[mode] = read_csv(f)
    if not tabs:
        log(f"  skip Fig1 ({task_label}): roi_v2_confirmatory.csv not found for any RUN_MODE"); return

    fig, axes = plt.subplots(1, len(CONFIRMATORY_HYPOTHESES), figsize=(FIG_W_CM * 0.65 * CM, 6.0 * CM), sharex=False)
    axes = np.atleast_1d(axes)
    for ax, (hyp, roi, band, contrast) in zip(axes, CONFIRMATORY_HYPOTHESES):
        disp_contrast = CONTRAST_DISPLAY_FLIP.get(contrast, contrast)
        rows = []
        for m in ("LogAbsolutePower", "RelativePower"):
            vals = {}
            for mode, t in tabs.items():
                r = t[(t["Hypothesis"] == hyp) & (t["Metric"] == m)]
                if len(r):
                    r = r.iloc[0]
                    _, e, lo, hi = flip_sign(contrast, r["Estimate"], r["CI_Low"] if "CI_Low" in r else np.nan,
                                              r["CI_High"] if "CI_High" in r else np.nan)
                    vals[mode] = (e, lo, hi, r["pValue"])
            rows.append((METRICS[m][0], vals, False))
        _forest(ax, rows, f"{hyp}: {roi}, {BAND_LABEL.get(band, band)}\ncontrast = {disp_contrast}", "Estimate")
    _forest_legend(fig); fig.suptitle(f"Confirmatory ROI v2 hypotheses (X24/X25) — Task={task_label}", y=1.02)
    fig.text(0.0, -0.16, "Sign flipped for display to E−R/E−N (this pipeline's native MATLAB output, "
              "R−E/N−E, is unchanged on disk).", fontsize=5.5, color="0.4")
    fig.tight_layout()
    save(fig, paths.main(task), "Fig1_confirmatory_H1_H2")


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


def fig2_paired(bp_e, bp_g, rois, mode_excl, task_label, out_main: Path):
    """Fig 2: subjects who saw both Reels and Neutral for this Task (descriptive, uncorrected)."""
    data = _paired_values(bp_e, bp_g, rois)
    if not data:
        log(f"  skip Fig2 ({task_label}): no E&N subjects found"); return
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
    fig.text(0.5, -0.04, f"Task={task_label}. Lines = individual subjects with both conditions; squares = mean "
             "± 95% CI. Descriptive, not corrected for multiplicity.", ha="center", fontsize=6)
    fig.tight_layout()
    save(fig, out_main, "Fig2_paired_subjects_EN")


def fig3_spectra(paths: Paths, task: str, task_label: str, labels, rois):
    """Fig 3: posterior spectra per condition and the within-subject E-N difference, this Task only."""
    from scipy.io import loadmat
    f = paths.inputs / "psd_recordings.mat"
    if not f.exists():
        log(f"  skip Fig3 ({task_label}): psd_recordings.mat missing (run F00_export_for_python_cognitive.m)"); return
    m = loadmat(f, squeeze_me=True)
    psd, fr = np.asarray(m["PSD"], float), np.asarray(m["FreqHz"], float).ravel()
    chl = [str(x).strip().upper() for x in np.atleast_1d(m["ChannelLabels"])]
    subj = np.atleast_1d(m["RecSubject"]).astype(int)
    cond = np.array([str(x).strip() for x in np.atleast_1d(m["RecCondition"])])
    tasktype = np.array([str(x).strip().lower() for x in np.atleast_1d(m["RecTaskType"])])
    keepTask = tasktype == task
    psd, subj, cond = psd[keepTask], subj[keepTask], cond[keepTask]
    if psd.shape[0] == 0:
        log(f"  skip Fig3 ({task_label}): no recordings for this Task in psd_recordings.mat"); return

    want = set()
    for r in POSTERIOR_ROIS:
        want |= {labels[i].upper() for i in rois.get(r, [])}
    ci = [k for k, l in enumerate(chl) if l in want] or list(range(len(chl)))
    lp = np.log10(psd[:, ci, :].mean(axis=1))                 # recording x freq
    keep = fr <= 45
    fr_k, lp = fr[keep], lp[:, keep]

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
        a1.plot(fr_k, mu, c=COND_COLORS[c], lw=1, label=f"{COND_NAMES[c]} (n = {len(x)})")
        a1.fill_between(fr_k, mu - 1.96 * se, mu + 1.96 * se, color=COND_COLORS[c], alpha=0.2, lw=0)
    a1.set_ylabel("log₁₀ power (µV²/Hz)")
    a1.legend(frameon=False, loc="upper right")
    a1.set_title(f"Task={task_label}. Parietal + occipital ROI (recording means, 95% CI)", pad=10)
    # within-subject difference, E&N subjects only, this Task only
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
        a2.plot(fr_k, mu, c="k", lw=1)
        a2.fill_between(fr_k, mu - q * se, mu + q * se, color="k", alpha=0.18, lw=0)
        a2.axhline(0, c="0.5", lw=0.6)
        a2.set_title(f"Within-subject Reels − Neutral (n = {len(d)} subjects with both; 95% CI)")
    else:
        log(f"  Fig3 ({task_label}): fewer than 3 subjects with both E and N for this Task -- difference panel thin/empty.")
    a2.set_ylabel("Δ log₁₀ power")
    a2.set_xlabel("Frequency (Hz)")
    a2.set_xlim(fr_k.min(), fr_k.max())
    save(fig, paths.main(task), "Fig3_posterior_spectra")


def _forest(ax, rows, title, xlabel):
    """rows: list of (label, {mode: (est, lo, hi, p)}, exploratory)."""
    y = np.arange(len(rows))[::-1]
    offs = {m: o for m, o in zip(RUN_MODES, (0.15, -0.15))}
    mk = {"full_sample": "o", "sensitivity_excl_noisy": "s"}
    for yi, (lab, vals, expl) in zip(y, rows):
        for mode, (e, lo, hi, p) in vals.items():
            col = "0.55" if expl else "k"
            fill = col if (p is not None and p < ALPHA) else "white"
            if np.isfinite(lo) and np.isfinite(hi):
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


def fig4_forest_global(paths: Paths, task: str, task_label: str):
    """Fig 4: whole-scalp estimates with 95% CI (X05B), both RUN_MODEs.
    CSV is read/filtered using the NATIVE 'N - E' label (GLOBAL_ROI_CONTRAST_
    LABEL -- that's what phase4_all_contrasts_ungated actually contains).
    By display convention, the PLOTTED Estimate/CI is sign-
    flipped to display as 'E - N' (GLOBAL_ROI_DISPLAY_LABEL), matching X09/
    X09B's native direction used elsewhere in this pipeline -- see
    flip_sign()/CONTRAST_DISPLAY_FLIP above. Nothing is re-derived; this is
    -(N-E), the same already-computed number, negated for display only."""
    tabs = {}
    for mode in RUN_MODES:
        f = paths.global_contrasts(task, mode)
        if f.exists():
            tabs[mode] = pd.read_excel(f)
    if not tabs:
        log(f"  skip Fig4 ({task_label}): phase4_all_contrasts_ungated.xlsx not found (run X05B)"); return
    fig, axes = plt.subplots(1, 2, figsize=(FIG_W_CM * CM, 7 * CM), sharey=True)
    for ax, metric in zip(axes, METRICS):
        rows = []
        for b in BANDS:
            vals = {}
            for mode, t in tabs.items():
                r = t[(t["Metric"] == metric) & (t["Band"] == b) & (t["Contrast"] == GLOBAL_ROI_CONTRAST_LABEL)]
                if len(r):
                    r = r.iloc[0]
                    _, e, lo, hi = flip_sign(GLOBAL_ROI_CONTRAST_LABEL, r["Estimate"], r["CI_Low"], r["CI_High"])
                    vals[mode] = (e, lo, hi, r["pValue"])
            rows.append((BAND_LABEL[b], vals, b in EXPLORATORY))
        _forest(ax, rows, f"{METRICS[metric][0]}: {GLOBAL_ROI_DISPLAY_LABEL}", f"Estimate ({METRICS[metric][1]})")
    _forest_legend(fig)
    fig.suptitle(f"Task={task_label}", y=1.02)
    fig.text(0.0, -0.16, f"Sign flipped for display: {GLOBAL_ROI_DISPLAY_LABEL} shown here = "
              f"−({GLOBAL_ROI_CONTRAST_LABEL}), this pipeline's native MATLAB output (unchanged on disk).",
              fontsize=5.5, color="0.4")
    fig.tight_layout()
    save(fig, paths.main(task), "Fig4_forest_global")


def fig5_forest_roi(paths: Paths, task: str, task_label: str):
    """Fig 5: ROI v2 relative-power estimates (X25_cognitive), both RUN_MODEs.
    Same native-read/display-flip note as Fig4 -- see flip_sign() above."""
    tabs = {m: read_csv(paths.roi2(task, m)) for m in RUN_MODES if paths.roi2(task, m).exists()}
    if not tabs:
        log(f"  skip Fig5 ({task_label}): roi_v2_all_contrasts.csv not found"); return
    rois = list(dict.fromkeys(next(iter(tabs.values()))["ROI"]))
    bands = ["Delta", "Theta", "Alpha", "Beta", "Beta_low_exploratory", "Beta_high_exploratory", "Gamma_low_exploratory"]
    fig, axes = plt.subplots(1, len(rois), figsize=(FIG_W_CM * CM, 6.5 * CM), sharey=True)
    for ax, roi in zip(np.atleast_1d(axes), rois):
        rows = []
        for b in bands:
            vals = {}
            for mode, t in tabs.items():
                r = t[(t["ROI"] == roi) & (t["Metric"] == "RelativePower") & (t["Band"] == b) & (t["Contrast"] == GLOBAL_ROI_CONTRAST_LABEL)]
                if len(r):
                    r = r.iloc[0]
                    _, e, lo, hi = flip_sign(GLOBAL_ROI_CONTRAST_LABEL, r["Estimate"], r["CI_Low"], r["CI_High"])
                    vals[mode] = (e, lo, hi, r["pValue"])
            rows.append((BAND_LABEL[b], vals, b in EXPLORATORY))
        _forest(ax, rows, roi, "Δ relative power")
    _forest_legend(fig)
    fig.suptitle(f"Task={task_label} | contrast = {GLOBAL_ROI_DISPLAY_LABEL}", y=1.02)
    fig.text(0.0, -0.16, f"Sign flipped for display: {GLOBAL_ROI_DISPLAY_LABEL} shown here = "
              f"−({GLOBAL_ROI_CONTRAST_LABEL}), this pipeline's native MATLAB output (unchanged on disk).",
              fontsize=5.5, color="0.4")
    fig.tight_layout()
    save(fig, paths.main(task), "Fig5_forest_ROIv2")


def fig6_sensitivity(stats, pos, task_label: str, out_main: Path):
    """Fig 6: E-N t-maps (native electrode-level convention), FullSample vs
    Sensitivity side by side, same scale."""
    sel = [("Theta", "RelativePower"), ("Alpha", "RelativePower"), ("Gamma_low_exploratory", "RelativePower")]
    modes = [m for m in RUN_MODES if m in stats]
    if len(modes) < 2:
        log(f"  skip Fig6 ({task_label}): needs X09/X09B output for both RUN_MODEs"); return
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
    add_cbar(fig, im, cax, "t, Reels − Neutral (relative power)", extend="both")
    fig.suptitle(f"Task={task_label}", y=1.0)
    fig.legend(handles=_marker_legend(), loc="lower center", ncol=1, frameon=False, bbox_to_anchor=(0.55, -0.06))
    save(fig, out_main, "Fig6_EminusN_FullSample_vs_Sensitivity")


# ---------------------------------------------------------------------------
# 6. PER-TASK DRIVER
# ---------------------------------------------------------------------------
def run_task(task: str, paths: Paths, labels, pos, rois, bp_e_all, bp_g_all, excl, bands_arg):
    task_label = paths.task_label(task)
    log(f"\n########################################################\n### TASK = {task_label}\n########################################################\n")

    stats = {}
    for mode in RUN_MODES:
        if paths.elec(task, mode, "electrode_model_estimates.csv").exists():
            stats[mode] = ElecStats(paths, task, mode, labels)
        else:
            log(f"  (no X09B output for Task={task_label}, {mode} -- that RUN_MODE is skipped)")
    if "full_sample" not in stats:
        log(f"  *** Task={task_label}: X09B full_sample output is required and missing -- skipping this Task entirely. ***")
        return

    log(f"Topo panels (2 x 3), Task={task_label} ...")
    bands = BANDS if bands_arg == "all" else [b.strip() for b in bands_arg.split(",")]
    for mode, st in stats.items():
        for b in bands:
            for me in METRICS:
                fig_topo_panel(st, pos, b, me, mode, task_label, paths.panels(task))

    log(f"Main figures, Task={task_label} ...")
    fig1_confirmatory(paths, task, task_label)

    bp_e = bp_e_all[bp_e_all["TaskType"].str.lower() == task].copy()
    bp_g = bp_g_all[bp_g_all["TaskType"].str.lower() == task].copy()
    if bp_e.empty or bp_g.empty:
        log(f"  Task={task_label}: no rows in bandpower_electrode/global.csv for this Task -- Fig2/Fig3 may be empty.")
    fig2_paired(bp_e, bp_g, rois, excl, task_label, paths.main(task))
    fig3_spectra(paths, task, task_label, labels, rois)
    fig4_forest_global(paths, task, task_label)
    fig5_forest_roi(paths, task, task_label)
    fig6_sensitivity(stats, pos, task_label, paths.main(task))
    check_montage(labels, pos, paths.out(task))


# ---------------------------------------------------------------------------
# 7. MAIN
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--root", default=DEFAULT_ROOT, help="ANALYSIS_ROOT of the cognitive_task pipeline")
    ap.add_argument("--bands", default="all", help="comma list for the topo panels, or 'all'")
    ap.add_argument("--exclude", default="7,23", help="subjects of the sensitivity run (dotted lines in Fig 2)")
    ap.add_argument("--tasks", default="all", help="comma list of full,nback,sart or 'all'")
    a = ap.parse_args()
    import mne
    log(f"F01 (cognitive_task) | mne {mne.__version__} | root = {a.root}")
    paths = Paths(Path(a.root))
    labels, pos = load_positions(paths)
    rois = roi_indices(paths, labels)

    bp_e_all = read_csv(need(paths.inputs / "bandpower_electrode.csv", "run F00_export_for_python_cognitive.m"))
    bp_g_all = read_csv(need(paths.inputs / "bandpower_global.csv", "run F00_export_for_python_cognitive.m"))
    for d in (bp_e_all, bp_g_all):
        d["Subject"] = d["Subject"].astype(str).str.replace(r"\D", "", regex=True).astype(int)
        d["TaskType"] = d["TaskType"].astype(str).str.lower().str.strip()
    idx = {l.upper(): k for k, l in enumerate(labels)}
    bp_e_all["ChanIdx"] = uc(bp_e_all["Channel"]).map(idx)
    excl = {int(x) for x in a.exclude.split(",") if x.strip()}

    tasks = list(TASKS) if a.tasks == "all" else [t.strip() for t in a.tasks.split(",")]
    for task in tasks:
        if task not in TASKS:
            sys.exit(f"ERROR: unknown task '{task}', expected one of {list(TASKS)}")
        run_task(task, paths, labels, pos, rois, bp_e_all, bp_g_all, excl, a.bands)

    log_path = paths.root / "31_publication_figures_cognitive_figure_log.txt"
    log_path.write_text("\n".join(LOG) + "\n")
    log(f"\nDone. Figures written under each Task's own <ANALYSIS_ROOT>/<Task>/31_publication_figures/. "
        f"Combined log: {log_path}")


if __name__ == "__main__":
    main()
