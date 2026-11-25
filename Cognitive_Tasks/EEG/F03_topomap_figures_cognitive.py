#!/usr/bin/env python3
"""
F03_topomap_figures_cognitive.py  --  cognitive_task pipeline, full-scalp
topographic figures for the X11/X12/X17/X18 "exploratory add-on" families.

WHY THIS SCRIPT EXISTS: the forest-plot
figures produced by F02_followup_family_figures_cognitive.py are accurate but
harder to read at a glance than a scalp map for a paper. F02's own CSVs
(the real, BH-FDR-gated pairwise outputs of X11/X12/X17/X18) only contain the
electrodes that SURVIVED FDR -- they cannot feed a continuous full-scalp heat
map, because a map needs a value at EVERY electrode, not just the survivors.
That is exactly why the X11B/X12B/X17B/X18B companion scripts were built: each
refits the SAME locked model at every one of the ~62 electrodes (ungated,
cross-checked against its own parent script's saved omnibus result) and
exports one electrode_all_contrasts_ungated.csv per family. THIS script reads
those four new CSVs (plus the parent scripts' own already-existing, FDR-gated
omnibus/pairwise CSVs, used ONLY to decide which electrodes get a filled vs.
open significance marker) and draws one scalp map per electrode, instead of
one forest-plot row per electrode.

NO new statistics are computed here. Every Estimate/SE/t/p plotted is copied
verbatim from a MATLAB-produced CSV. The only things computed in this script
are: (a) t = Estimate / SE for colour scaling (identical in spirit to F01's
own t-maps), and (b) which marker (filled/open/none) an electrode gets, from
already-computed p-values -- never a new p-value, never a new FDR correction.

===========================================================================
WHAT IS DIFFERENT FROM F01's OWN TOPO PANELS (X09/X09B) -- FULL DISCLOSURE
===========================================================================

1. NO continuous-condition (EMM) maps. F01's Fig6/topo panels show both the
   per-condition estimated mean (top row) AND the pairwise t-map (bottom row)
   because X09B's own saved table carries EMM_<cond> columns. X11B/X12B/
   X17B/X18B do NOT save per-condition EMMs -- their parent scripts (X11/
   X12/X17/X18) only ever fit and report CONTRAST estimates (differences),
   never raw condition means, at the electrode level. So every map in this
   script is a contrast/difference map (Estimate, in the metric's own
   units), coloured with a diverging colormap, exactly like F01's own
   t-map row -- there is no sequential-colormap "mean level" row here.
   This is a direct, disclosed consequence of what the companion scripts
   could legitimately re-derive from already-saved data, not a stylistic
   choice.

2. CONTRAST SIGN CONVENTION: all four families (X11,
   X12, X17, X18) still store N-E / R-E / R-N natively on disk (REF_LEVEL='E',
   verified by reading each script's own CONTRASTS literal directly: X11
   line ~290, X12 line ~281, X17 line ~264, X18 line ~738, all "{'N','E','N -
   E'; 'R','E','R - E'; 'R','N','R - N'}") -- the GLOBAL/ROI convention,
   opposite sign from X09/X09B's native E-N/E-R/N-R used by F01's own topo
   panels (see F01's own header item 2). By design,
   this script DISPLAYS every value under the E-N/E-R/N-R
   label instead, so every figure across this whole pipeline reads the same
   direction: E-N := -(N-E), CI swapped-and-negated -- see
   CONTRAST_DISPLAY_FLIP / flip_sign_for_display() below. This is a pure sign
   flip of the ALREADY-COMPUTED Estimate, done only in this Python script;
   nothing in MATLAB is touched, nothing is re-derived, p-values/markers are
   unaffected. CSV lookups (which electrode gets a filled/open marker) always
   query the on-disk file using the NATIVE label, since that is what the CSV
   actually contains -- only the plotted color/title/__rows.csv *_display
   columns use the flipped label+sign. "Second - First" / "TimeHalf=... -
   TimeHalf=..." (X18 Position / X17 TimeHalf marginal axes) are NOT
   Condition contrasts and are never flipped.

3. SIGNIFICANCE MARKERS are necessarily two-tier and weaker than F01's own,
   because the companion CSVs are explicitly "ungated" / non-inferential
   (see each XxxB script's own header):
     - FILLED marker ("passed the real BH-FDR + Holm gate"): the electrode
       appears in the ORIGINAL (not companion) gated pairwise CSV for that
       exact Band/Metric/Contrast(/ContrastFamily), AND that CSV's own pHolm
       < ALPHA. This is the one and only truly inferential marker in this
       script -- identical evidentiary standard to F02's forest plots.
     - OPEN marker ("nominal, uncorrected only"): the electrode's own
       companion-CSV pValue < ALPHA but it is not FILLED (either because the
       omnibus/gated step did not run yet, or it did run and this specific
       electrode did not survive). This mirrors F01's own "loose" marker,
       EXCEPT F01 also requires the electrode's own omnibus to be
       FDR-significant for the "loose" tier; here there is no cheap way to
       re-derive that without an omnibus CSV lookup, so "open" here is purely
       "ungated p < 0.05", clearly weaker than F01's -- labelled as such in
       every legend and in the figure caption.
     - NO marker: pValue >= ALPHA, or gated CSV not found/not matched.
   If the ORIGINAL gated pairwise CSV for a given family/Task/RUN_MODE/
   contrast is simply missing from disk, this is NOT an error -- it means
   nothing survived BH-FDR+Holm for that combination (same convention as
   F02), and every marker in that panel falls back to open-or-none only.
   This is logged explicitly in the manifest, never silently.

4. MONTAGE / PROJECTION / HEAD-DRAWING CODE is copied VERBATIM from
   F01_publication_figures_cognitive.py (load_positions, check_montage,
   _draw_head, topo, robust_limits, add_cbar, _cax, save) so that every
   figure in this pipeline -- F01's own topo panels and this script's new
   ones -- uses the exact same projection, same head outline, same
   orientation convention, same font/colour settings. Nothing about the
   montage math was re-derived or changed.

5. X12 is at the CHANNEL (not Electrode) level in its own source column
   naming (band-ratio indices per electrode), but it is the SAME 62-channel
   montage as the other three families (verified: X12's own ELECTRODES list
   is Task-invariant, same montage source as X09/X11/X17/X18). Its own group
   dimensions are Variant (MeanBased/MedianBased) x RatioType (ThetaBeta/
   AlphaTheta/AlphaBeta/EngagementIndex) instead of Band x Metric.

RUN ORDER
  MATLAB (once per Task x RUN_MODE, already built/verified this round):
    X11B_electrode_model_estimates_cognitive.m
    X12B_electrode_ratio_estimates_cognitive.m
    X17B_electrode_all_contrasts_ungated_cognitive.m   (Nback/SART only)
    X18B_electrode_all_contrasts_ungated_cognitive.m   (Nback/SART only)
  (no rerun of X11/X12/X17/X18 or anything upstream is required -- see the
  companion scripts' own headers)
  Python:
    python3 F00_export_for_python_cognitive.m   [MATLAB -- for channel_locations.csv, Task-invariant, already required by F01]
    python3 F03_topomap_figures_cognitive.py [--root <ANALYSIS_ROOT>] [--tasks all] [--modes all]

OUTPUT   <ANALYSIS_ROOT>/<Task>/33_topomap_figures/<X11|X12|X17|X18>/<RUN_MODE>/
           one <group>__<ContrastFamily-if-any>.png/.pdf per Band(or
           Variant/RatioType) x Metric combination, each with a companion
           __rows.csv (the exact rows plotted, for one-by-one traceability,
           same convention as F02's __rows.csv)
         33_topomap_figures_coverage_manifest.txt  (one line per
           Task x RUN_MODE x family x group: FOUND / none survived FDR here
           (no filled markers possible) / ungated data not found (run the
           XxxB companion) / OUT OF SCOPE (Full Task, X17/X18))
         33_topomap_figures_figure_log.txt
"""

from __future__ import annotations

import argparse
import re
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

# ---------------------------------------------------------------------------
# 0. SETTINGS
# ---------------------------------------------------------------------------
DEFAULT_ROOT = "results"  # pass --root to point at your own ANALYSIS_ROOT

TASKS = {"full": "Full", "nback": "Nback", "sart": "SART"}
RUN_MODES_ORDER = ["sensitivity_excl_noisy", "full_sample"]   # sensitivity listed first, by design
# NOTE: the subject numbers named in the "Sensitivity" label below are a
# display label only; EXCLUSION_CONFIG_cognitive.m's EXCLUDE_NOISY_SUBJECTS
# is this pipeline's single source of truth for which subjects are actually
# excluded in sensitivity_excl_noisy mode. Verify this label text against
# that file's current value before publication -- it is not re-derived from
# it automatically.
RUN_MODES = {"full_sample": "Full sample", "sensitivity_excl_noisy": "Sensitivity (excl. 7, 23)"}

BAND_ORDER = ["Delta", "Theta", "Alpha", "Beta", "Beta_low_exploratory", "Beta_high_exploratory",
              "Gamma_low_exploratory", "Gamma_high1_exploratory", "Gamma_high2_exploratory"]
BAND_LABEL = {"Delta": "Delta 1-4 Hz", "Theta": "Theta 4-8 Hz", "Alpha": "Alpha 8-13 Hz",
              "Beta": "Beta 13-30 Hz", "Beta_low_exploratory": "Beta low 13-20 Hz*",
              "Beta_high_exploratory": "Beta high 20-30 Hz*",
              "Gamma_low_exploratory": "Gamma 30-45 Hz*", "Gamma_high1_exploratory": "Gamma 55-80 Hz*",
              "Gamma_high2_exploratory": "Gamma 80-100 Hz*"}

# CONTRAST_ORDER stays NATIVE (N-E/R-E/R-N, this pipeline's on-disk convention
# for X11/X12/X17/X18) -- it is only used to pick a stable panel order and to
# drive the CSV lookups in build_marker_lookup(), which must match the
# Contrast strings actually on disk. CONTRAST_COLOR below is DISPLAY-keyed
# (E-N/E-R/N-R) -- see CONTRAST_DISPLAY_FLIP / flip_sign_for_display().
CONTRAST_ORDER = ["N - E", "R - E", "R - N"]
# Okabe-Ito, colour-blind safe -- same palette family as F01/F02
CONTRAST_COLOR = {"E - N": "#D55E00", "E - R": "#0072B2", "N - R": "#009E73"}

# ---------------------------------------------------------------------------
# DISPLAY-ONLY SIGN FLIP --
# identical logic/rationale to F02_followup_family_figures_cognitive.py's own
# CONTRAST_DISPLAY_FLIP/flip_sign(). The ungated companion CSVs (X11B/X12B/
# X17B/X18B) are natively N-E/R-E/R-N, same as their parent scripts. By
# design, every topomap here is displayed under the X09/X09B-
# native E-N/E-R/N-R label instead: E-N := -(N-E), CI swapped-and-negated.
# NOTHING in MATLAB is touched -- this negates the already-computed Estimate
# purely for the plotted colour/marker and the printed title, in Python only.
# CSV lookups in build_marker_lookup() ALWAYS use the native label (clabel),
# since that is what is actually written in the on-disk CSVs -- only the
# final display (panel title, grid row/col labels, plotted Estimate sign,
# and the extra *_display columns in each __rows.csv) is flipped.
CONTRAST_DISPLAY_FLIP = {"N - E": "E - N", "R - E": "E - R", "R - N": "N - R"}
CONTRAST_ORDER_DISPLAY = ["E - N", "E - R", "N - R"]


def display_contrast(native_label: str) -> str:
    """Flips only the base of a '<base> @ <Dim>=<level>' label (X17/X18
    InteractionSimple); labels with no flip entry (TimeHalf/Position
    marginal contrasts) are returned unchanged."""
    if native_label in CONTRAST_DISPLAY_FLIP:
        return CONTRAST_DISPLAY_FLIP[native_label]
    if " @ " in str(native_label):
        base, suffix = native_label.split(" @ ", 1)
        if base in CONTRAST_DISPLAY_FLIP:
            return f"{CONTRAST_DISPLAY_FLIP[base]} @ {suffix}"
    return native_label


def flip_sign_for_display(native_label: str, vals):
    """vals: array of Estimates for this contrast, one per electrode.
    Returns (display_label, display_vals) -- negated iff native_label (or its
    '<base> @ ...' base) is in CONTRAST_DISPLAY_FLIP."""
    disp = display_contrast(native_label)
    if disp == native_label:
        return disp, vals
    return disp, -vals


ROBUST_PCT = 98
T_MIN_LIMIT = 2.0
ALPHA = 0.05
FIG_W_CM = 17.8
FONT = "DejaVu Sans"
CMAP_DIV = "RdBu_r"
HEAD_R = 0.095
CM = 1 / 2.54

plt.rcParams.update({
    "font.family": FONT, "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
    "ps.fonttype": 42, "savefig.dpi": 300, "figure.dpi": 100})

MANIFEST = []
FIGLOG = []


def log_manifest(line: str) -> None:
    print(line)
    MANIFEST.append(line)


def log_fig(line: str) -> None:
    print(line)
    FIGLOG.append(line)


def read_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    for c in df.columns:
        if df[c].dtype == object:
            df[c] = df[c].astype(str).str.strip()
    return df


def uc(s):
    return s.astype(str).str.upper().str.strip()


# ---------------------------------------------------------------------------
# 1. MONTAGE / TOPOMAP PRIMITIVES -- copied VERBATIM from
#    F01_publication_figures_cognitive.py (same projection, same head
#    outline, same orientation convention). See disclosure item 4 above.
# ---------------------------------------------------------------------------
def load_positions(inputs_dir: Path):
    f = inputs_dir / "channel_locations.csv"
    if not f.exists():
        sys.exit(f"ERROR: missing input {f}\n  -> run F00_export_for_python_cognitive.m in MATLAB")
    loc = read_csv(f)
    xyz = loc[["X", "Y", "Z"]].to_numpy(float)
    xyz = xyz / np.linalg.norm(xyz, axis=1, keepdims=True)
    colat = np.arccos(np.clip(xyz[:, 2], -1, 1))
    az = np.arctan2(xyz[:, 1], xyz[:, 0])
    r = colat / (np.pi / 2)
    r = r / max(r.max() / 0.92, 1.0) * HEAD_R
    pos = np.column_stack([-r * np.sin(az), r * np.cos(az)])
    labels = loc["Label"].astype(str).str.strip().tolist()
    return labels, pos


def _draw_head(ax, lw=0.6):
    t = np.linspace(0, 2 * np.pi, 200)
    ax.plot(HEAD_R * np.cos(t), HEAD_R * np.sin(t), "k", lw=lw)
    ax.plot([-0.09 * HEAD_R, 0, 0.09 * HEAD_R], [0.995 * HEAD_R, 1.09 * HEAD_R, 0.995 * HEAD_R], "k", lw=lw)
    for s in (-1, 1):
        ax.plot(s * HEAD_R * np.array([0.995, 1.04, 1.06, 1.04, 0.995]),
                HEAD_R * np.array([0.12, 0.10, 0.0, -0.12, -0.15]), "k", lw=lw)


def topo(ax, values, pos, vlim, cmap, strict=None, loose=None, contours=3):
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
    return im


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
    ax = fig.add_subplot(spec)
    b = ax.get_position()
    ax.set_position([b.x0, b.y0 + b.height * (1 - frac) / 2, b.width, b.height * frac])
    return ax


def save(fig, folder: Path, stem: str):
    folder.mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(folder / f"{stem}.png", bbox_inches="tight", dpi=300)
    plt.close(fig)
    log_fig(f"  wrote {folder}/{stem}.pdf/.png")


def _marker_legend():
    return [Line2D([], [], marker="o", ls="", mfc="k", mec="k", ms=3.5,
                    label="Survived BH-FDR (omnibus) + Holm (gated pairwise CSV)"),
            Line2D([], [], marker="o", ls="", mfc="none", mec="k", ms=3.5,
                    label="Ungated p < 0.05 only (companion script; non-inferential)")]


# ---------------------------------------------------------------------------
# 2. PATHS
# ---------------------------------------------------------------------------
FAMILY_FOLDER = {"X11": "11_median_robust_family", "X12": "12_band_ratio_indices",
                  "X17": "17_half_split_condition_time", "X18": "18_position_pos1_pos2"}


class Paths:
    def __init__(self, root: Path):
        self.root = root
        self.inputs = root / "30_python_inputs"

    def task_label(self, task: str) -> str:
        return TASKS[task]

    def family_dir(self, task: str, family: str, mode: str) -> Path:
        return self.root / self.task_label(task) / FAMILY_FOLDER[family] / mode

    def out(self, task: str, family: str, mode: str) -> Path:
        d = self.root / self.task_label(task) / "33_topomap_figures" / family / mode
        d.mkdir(parents=True, exist_ok=True)
        return d


# ---------------------------------------------------------------------------
# 3. FAMILY SPECS -- one entry per XxxB companion output. Each spec says
#    where to find the ungated (all-electrode) CSV, where to find the
#    ORIGINAL gated omnibus/pairwise CSVs (for markers only), which Tasks are
#    in scope, and which columns define one figure ("group_dims").
# ---------------------------------------------------------------------------
# orig_id_col: the electrode-identifying column name used by the ORIGINAL (not
# companion) omnibus/pairwise CSVs -- verified by reading each parent script's
# own writetable()'d column-name literal directly. This is NOT always the same
# name the companion (XxxB) script uses for its own output (id_col below):
#   X11/X17 use the pipeline's generic Level/Unit scheme (Unit = electrode name
#     at Level=='Electrode', verified: X11 line 588/740 "omnibusBandCols =
#     {'Level','Unit',...}" / "pairwiseBandCols = {'Level','Unit',...}"; X17
#     line 596 "omnibusCols = {'Level','Unit',...}" / line 760 "pairwiseCols =
#     {'Level','Unit',...}"), while X11B/X17B's OWN output already says
#     "Electrode" (clearer, single-level-only column name).
#   X12 uses 'Channel' in both the original ElecOmnibusT/ElecPairwiseT and its
#     own companion's output -- same name, no mapping needed.
#   X18 uses 'Electrode' in both the original OmnibusElectrodeT/pairwise
#     *_electrode.csv files and its own companion's output -- same name.
# id_col and orig_id_col must NOT be assumed equal: X11's omnibus/pairwise
# CSVs (phase11_omnibus_band_results.csv / phase11_pairwise_band_contrasts.csv)
# use the column name "Unit", not "Electrode", so the lookup below reads
# orig_id_col explicitly per family rather than hard-coding "Electrode".
FAMILY_SPECS = [
    dict(family="X11", id_col="Electrode", orig_id_col="Unit", group_dims=["Metric", "Band"],
         tasks=("full", "nback", "sart"),
         ungated_fname="electrode_all_contrasts_ungated.csv",
         omnibus_fname="phase11_omnibus_band_results.csv",
         omnibus_level_col="Level", omnibus_level_val="Electrode",
         omnibus_group_cols=["Metric", "Band"],
         pairwise_fname="phase11_pairwise_band_contrasts.csv",
         pairwise_level_col="Level", pairwise_level_val="Electrode",
         pairwise_group_cols=["Metric", "Band"],
         family_tag_col=None, metric_list=["AbsMedian", "RelMedian", "CV_Absolute", "CV_Relative"]),
    dict(family="X12", id_col="Channel", orig_id_col="Channel", group_dims=["Variant", "RatioType"],
         tasks=("full", "nback", "sart"),
         ungated_fname="electrode_all_contrasts_ungated.csv",
         omnibus_fname="phase12_band_ratio_electrode_omnibus.csv",
         omnibus_level_col=None, omnibus_level_val=None,
         omnibus_group_cols=["Variant", "RatioType"],
         pairwise_fname="phase12_band_ratio_electrode_pairwise.csv",
         pairwise_level_col=None, pairwise_level_val=None,
         pairwise_group_cols=["Variant", "RatioType"],
         family_tag_col=None, metric_list=None),
    dict(family="X17", id_col="Electrode", orig_id_col="Unit", group_dims=["Metric", "Band"],
         tasks=("nback", "sart"),
         ungated_fname="electrode_all_contrasts_ungated.csv",
         omnibus_fname="phase17_omnibus_results.csv",
         omnibus_level_col="Level", omnibus_level_val="Electrode",
         omnibus_group_cols=["Metric", "Band"],
         pairwise_fname="phase17_pairwise_contrasts.csv",
         pairwise_level_col="Level", pairwise_level_val="Electrode",
         pairwise_group_cols=["Metric", "Band"],
         family_tag_col="ContrastFamily",
         # companion's own ContrastFamily tag -> parent script's own omnibus "Term" value
         # (NOT the same strings -- verified by reading X17's own TERM_KEYS vs X17B's own
         # allRows tags directly; see module docstring item 3)
         omnibus_term_col="Term",
         family_tag_to_term={"ConditionMarginal": "ConditionMain", "TimeHalfMarginal": "TimeHalfMain",
                              "InteractionSimple": "Interaction"},
         pairwise_term_col="SurvivingTerm",
         family_tag_to_pairwise_term={"ConditionMarginal": "ConditionMain", "TimeHalfMarginal": "TimeHalfMain",
                                       "InteractionSimple": "Interaction"},
         metric_list=None),
    dict(family="X18", id_col="Electrode", orig_id_col="Electrode", group_dims=["Metric", "Band"],
         tasks=("nback", "sart"),
         ungated_fname="electrode_all_contrasts_ungated.csv",
         omnibus_fname="phase18_omnibus_electrode.csv",
         omnibus_level_col=None, omnibus_level_val=None,
         omnibus_group_cols=["Metric", "Band"],
         pairwise_fname=None,   # three separate files, chosen per family_tag below
         pairwise_level_col=None, pairwise_level_val=None,
         pairwise_group_cols=["Metric", "Band"],
         family_tag_col="ContrastFamily",
         omnibus_term_col="Term",
         family_tag_to_term={"ConditionMarginal": "Condition", "PositionMarginal": "Position",
                              "InteractionSimple": "Condition:Position"},
         pairwise_fname_by_tag={"ConditionMarginal": "phase18_pairwise_condition_marginal_electrode.csv",
                                 "PositionMarginal": "phase18_pairwise_position_marginal_electrode.csv",
                                 "InteractionSimple": "phase18_pairwise_interaction_simple_electrode.csv"},
         metric_list=None),
]


# ---------------------------------------------------------------------------
# 4. MARKER LOOKUP
# ---------------------------------------------------------------------------
def build_marker_lookup(spec, fdir: Path, group_key: dict, contrast_label: str, family_tag):
    """Returns (survived_fdr: set[ElectrodeUC] or None, holm_df: DataFrame or None).
    survived_fdr is None if the omnibus CSV for this family isn't on disk (no way to
    know which electrodes even passed the omnibus step); holm_df is None if the
    gated pairwise CSV for this exact group/contrast isn't on disk (no filled
    markers possible for this panel -- logged, not an error)."""
    omn_path = fdir / spec["omnibus_fname"]
    survived = None
    if omn_path.exists():
        omn = read_csv(omn_path)
        mask = pd.Series(True, index=omn.index)
        if spec.get("omnibus_level_col"):
            mask &= omn[spec["omnibus_level_col"]] == spec["omnibus_level_val"]
        for c in spec["omnibus_group_cols"]:
            mask &= omn[c].astype(str) == str(group_key[c])
        if spec.get("omnibus_term_col") and family_tag is not None:
            term = spec["family_tag_to_term"].get(family_tag)
            if term is not None:
                mask &= omn[spec["omnibus_term_col"]].astype(str) == term
        sub = omn[mask]
        idc = spec["orig_id_col"]
        if idc in sub.columns and "SignificantFDR" in sub.columns:
            fdr_col = sub["SignificantFDR"]
            if fdr_col.dtype == object:
                fdr_bool = fdr_col.astype(str).str.lower().isin(["1", "true"])
            else:
                fdr_bool = fdr_col.astype(bool)
            survived = set(uc(sub.loc[fdr_bool, idc]))
        else:
            survived = set()  # no per-electrode id / no FDR column in this omnibus row shape -- treat as "none known"

    # gated pairwise CSV -> electrodes that have pHolm < ALPHA for this exact contrast
    if spec["family"] == "X18":
        pw_fname = spec["pairwise_fname_by_tag"].get(family_tag)
    else:
        pw_fname = spec["pairwise_fname"]
    holm_elecs = None
    if pw_fname is not None:
        pw_path = fdir / pw_fname
        if pw_path.exists():
            pw = read_csv(pw_path)
            mask = pd.Series(True, index=pw.index)
            if spec.get("pairwise_level_col"):
                mask &= pw[spec["pairwise_level_col"]] == spec["pairwise_level_val"]
            for c in spec["pairwise_group_cols"]:
                mask &= pw[c].astype(str) == str(group_key[c])
            if spec.get("pairwise_term_col") and family_tag is not None:
                pterm = spec["family_tag_to_pairwise_term"].get(family_tag)
                if pterm is not None:
                    mask &= pw[spec["pairwise_term_col"]].astype(str) == pterm
            mask &= pw["Contrast"].astype(str) == contrast_label
            sub = pw[mask]
            idc = spec["orig_id_col"]
            holm_elecs = set(uc(sub[sub["pHolm"] < ALPHA][idc])) if idc in sub.columns else set()
    return survived, holm_elecs


# ---------------------------------------------------------------------------
# 5. ONE FAMILY, ONE TASK, ONE RUN_MODE
# ---------------------------------------------------------------------------
def process_family(paths: Paths, spec: dict, task: str, mode: str, labels, pos):
    family = spec["family"]
    task_label = paths.task_label(task)
    tag = f"{family} | Task={task_label} | {RUN_MODES[mode]}"
    if task not in spec["tasks"]:
        log_manifest(f"{tag}: OUT OF SCOPE (this family's own TASK_TYPES does not include {task_label})")
        return
    fdir = paths.family_dir(task, family, mode)
    ung_path = fdir / spec["ungated_fname"]
    if not ung_path.exists():
        log_manifest(f"{tag}: ungated data not found at {ung_path} "
                      f"-- run the X{family[1:]}B companion script for this Task/RUN_MODE first")
        return
    df = read_csv(ung_path)
    df["IdUC"] = uc(df[spec["id_col"]])
    labels_uc = [l.upper() for l in labels]
    missing = set(labels_uc) - set(df["IdUC"])
    if missing:
        log_manifest(f"{tag}: ERROR -- {len(missing)} montage electrodes absent from the companion's own "
                      f"output ({sorted(missing)[:5]}{'...' if len(missing) > 5 else ''}); skipping this family "
                      f"for this Task/RUN_MODE rather than plotting an incomplete scalp map.")
        return

    group_cols = spec["group_dims"]
    tag_col = spec["family_tag_col"]
    groupers = group_cols + ([tag_col] if tag_col else [])
    n_figs = 0
    for key_vals, g in df.groupby(groupers, dropna=False):
        key_vals = key_vals if isinstance(key_vals, tuple) else (key_vals,)
        group_key = dict(zip(group_cols, key_vals[:len(group_cols)]))
        family_tag = key_vals[len(group_cols)] if tag_col else None

        contrasts_present = [c for c in g["Contrast"].unique()]
        # keep a stable, readable order: the 3 canonical ones first if present, then the rest as-is
        ordered = [c for c in CONTRAST_ORDER if c in contrasts_present] + \
                  [c for c in contrasts_present if c not in CONTRAST_ORDER]
        if not ordered:
            continue

        panel_data = []  # (contrast_label, values[62], strict_mask, loose_mask, holm_found: bool)
        for clabel in ordered:
            gc = g[g["Contrast"] == clabel].set_index("IdUC").reindex(labels_uc)
            vals = gc["Estimate"].to_numpy(float)
            if np.all(~np.isfinite(vals)):
                continue
            survived, holm_elecs = build_marker_lookup(spec, fdir, group_key, clabel, family_tag)
            p_unc = gc["pValue"].to_numpy(float)
            if holm_elecs is not None:
                strict = np.array([l in holm_elecs for l in labels_uc])
                if survived is not None:
                    strict = strict & np.array([l in survived for l in labels_uc])
            else:
                strict = np.zeros(len(labels_uc), dtype=bool)
            loose = (p_unc < ALPHA) & ~strict
            # CSV lookups above used the NATIVE clabel (matches what's on disk);
            # from here on only the DISPLAY label/sign are used (plot title,
            # grid row/col text, plotted Estimate colour) -- see
            # CONTRAST_DISPLAY_FLIP / flip_sign_for_display().
            disp_label, disp_vals = flip_sign_for_display(clabel, vals)
            panel_data.append((disp_label, disp_vals, strict, loose, holm_elecs is not None, clabel))

        if not panel_data:
            continue

        gated_missing_any = any(not found for *_, found, _nat in panel_data)
        group_name = "_".join(str(v) for v in group_key.values())
        fam_suffix = f"__{family_tag}" if family_tag else ""
        stem = f"{group_name}{fam_suffix}"

        # Layout: a "simple-effects" family (6 contrasts of the form
        # "<base contrast> @ <SplitDim>=<level>", e.g. X17's InteractionSimple
        # "N - E @ TimeHalf=TimeHalf1" or X18's "N - E @ Position=First") is
        # laid out as a 2D grid -- rows = split level, cols = base contrast --
        # instead of cramming 6 panels into one row, which made titles overlap
        # into unreadable text. Every other family (<=3 flat contrasts) keeps
        # the original single-row layout.
        pat = re.compile(r"^(.*) @ (\w+)=(.+)$")
        parsed = [pat.match(clabel) for clabel, *_ in panel_data]
        is_grid = len(panel_data) > 3 and all(parsed)
        m = max(robust_limits([v for _, v, *_ in panel_data], symmetric=True)[1], 1e-6)
        vlim = (-m, m)

        if is_grid:
            bases, split_dim, levels = [], parsed[0].group(2), []
            for p in parsed:
                b, _, lvl = p.group(1), p.group(2), p.group(3)
                if b not in bases:
                    bases.append(b)
                if lvl not in levels:
                    levels.append(lvl)
            bases_ordered = [b for b in CONTRAST_ORDER_DISPLAY if b in bases] + \
                             [b for b in bases if b not in CONTRAST_ORDER_DISPLAY]
            grid = {(p.group(1), p.group(3)): pd for p, pd in zip(parsed, panel_data)}
            nrows, ncols = len(levels), len(bases_ordered)
            fig = plt.figure(figsize=(FIG_W_CM * CM, (3.6 * nrows + 1.0) * CM))
            gs = GridSpec(nrows, ncols + 1, figure=fig, width_ratios=[1] * ncols + [0.06],
                          wspace=0.10, hspace=0.38)
            im = None
            for ri, lvl in enumerate(levels):
                for ci, b in enumerate(bases_ordered):
                    clabel, vals, strict, loose, found, _nat = grid[(b, lvl)]
                    ax = fig.add_subplot(gs[ri, ci])
                    im = topo(ax, vals, pos, vlim, CMAP_DIV, strict, loose)
                    if ri == 0:
                        ax.set_title(b, color=CONTRAST_COLOR.get(b, "k"), fontsize=6.5)
                    if ci == 0:
                        ax.text(-1.35, 0.5, f"{split_dim}={lvl}", transform=ax.transAxes, fontsize=6.5,
                                ha="center", va="center", rotation=90)
            add_cbar(fig, im, _cax(fig, gs[:, ncols]), "Estimate\n(model-based difference)",
                     extend="both", orientation="vertical")
            suptitle_y = 1.0 + 0.10 * (2 / max(nrows, 1))
        else:
            ncols = len(panel_data)
            fig = plt.figure(figsize=(FIG_W_CM * CM, 7.0 * CM))
            gs = GridSpec(1, ncols + 1, figure=fig, width_ratios=[1] * ncols + [0.06], wspace=0.08)
            im = None
            for k, (clabel, vals, strict, loose, found, _nat) in enumerate(panel_data):
                ax = fig.add_subplot(gs[0, k])
                im = topo(ax, vals, pos, vlim, CMAP_DIV, strict, loose)
                suffix = "" if found else "  [no gated CSV]"
                ax.set_title(clabel + suffix, color=CONTRAST_COLOR.get(clabel, "k"), fontsize=6.5)
            add_cbar(fig, im, _cax(fig, gs[0, ncols]), "Estimate\n(model-based difference)",
                     extend="both", orientation="vertical")
            suptitle_y = 1.08

        group_label = ", ".join(f"{k}={v}" for k, v in group_key.items())
        if family_tag:
            group_label += f" [{family_tag}]"
        warn = "  (gated pairwise CSV missing for >=1 panel -- open markers only there)" \
            if gated_missing_any else ""
        fig.suptitle(f"{family} -- {group_label}\nTask={task_label} | {RUN_MODES[mode]}{warn}",
                     y=suptitle_y, fontsize=7.5)
        legend_y = -0.06 if is_grid else -0.14
        fig.legend(handles=_marker_legend(), loc="lower center", ncol=1, frameon=False,
                   bbox_to_anchor=(0.5, legend_y), fontsize=6)
        fig.text(0.0, legend_y - 0.09, "Sign flipped for display: E−N/E−R/N−R shown here = −(N−E)/−(R−E)/−(R−N),"
                  " this pipeline's native companion-script output (unchanged on disk). See build-record.",
                  fontsize=5, color="0.4")

        out_dir = paths.out(task, family, mode)
        save(fig, out_dir, stem)

        rows_out = g[g["Contrast"].isin(ordered)].copy()
        rows_out.drop(columns=["IdUC"], inplace=True, errors="ignore")
        rows_out["Contrast_display"] = rows_out["Contrast"].map(display_contrast)
        flipped_mask = rows_out["Contrast_display"] != rows_out["Contrast"]
        rows_out["Estimate_display"] = np.where(flipped_mask, -rows_out["Estimate"], rows_out["Estimate"])
        rows_out["CI_Low_display"] = np.where(flipped_mask, -rows_out["CI_High"], rows_out["CI_Low"])
        rows_out["CI_High_display"] = np.where(flipped_mask, -rows_out["CI_Low"], rows_out["CI_High"])
        rows_out.to_csv(out_dir / f"{stem}__rows.csv", index=False)
        n_figs += 1

    if n_figs == 0:
        log_manifest(f"{tag}: ungated data found but every group had all-NaN estimates -- nothing plotted")
    else:
        log_manifest(f"{tag}: FOUND -- {n_figs} figure(s) written to {paths.out(task, family, mode).parent}")


# ---------------------------------------------------------------------------
# 6. MAIN
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=DEFAULT_ROOT)
    ap.add_argument("--tasks", default="all", help="comma list of full,nback,sart or 'all'")
    ap.add_argument("--modes", default="all", help="comma list of RUN_MODE labels or 'all'")
    a = ap.parse_args()

    import mne
    print(f"F03 (cognitive_task) | mne {mne.__version__} | root = {a.root}")

    root = Path(a.root)
    paths = Paths(root)
    labels, pos = load_positions(paths.inputs)
    print(f"Montage: {len(labels)} electrodes loaded from {paths.inputs / 'channel_locations.csv'}")

    tasks = list(TASKS) if a.tasks == "all" else [t.strip() for t in a.tasks.split(",")]
    modes = RUN_MODES_ORDER if a.modes == "all" else [m.strip() for m in a.modes.split(",")]
    modes = [m for m in RUN_MODES_ORDER if m in modes]  # keep sensitivity-first order

    for mode in modes:
        for task in tasks:
            for spec in FAMILY_SPECS:
                process_family(paths, spec, task, mode, labels, pos)

    manifest_path = root / "33_topomap_figures_coverage_manifest.txt"
    manifest_path.write_text("\n".join(MANIFEST) + "\n")
    figlog_path = root / "33_topomap_figures_figure_log.txt"
    figlog_path.write_text("\n".join(FIGLOG) + "\n")
    print(f"\nWrote {manifest_path}")
    print(f"Wrote {figlog_path}")


if __name__ == "__main__":
    main()
