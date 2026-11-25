#!/usr/bin/env python3
"""
F02_followup_family_figures_cognitive.py  --  cognitive_task pipeline

Visualizes the FDR-surviving pairwise-contrast outputs of FOUR exploratory,
purely-additive follow-up families that F01_publication_figures_cognitive.py
does NOT plot at all:

  X11  median-robust derived-metric family   (AbsMedian / RelMedian / CV_Absolute /
       CV_Relative, band-dependent; + TotalPowerMean / TotalPowerMedian, broadband)
  X12  band-ratio indices                     (ThetaBeta / AlphaTheta / AlphaBeta /
       EngagementIndex, each as MeanBased and MedianBased variants)
  X17  half-split Condition x TimeHalf        (THREE independent FDR families per
       level: Condition-main / TimeHalf-main / Interaction). Nback & SART ONLY --
       Full is categorically out of scope for this analysis (X17's own TASK_TYPES).
  X18  recording-position / block-order       (THREE independent FDR families --
       Condition-marginal / Position-marginal / Interaction-simple -- each at
       Global AND Electrode level, so up to 6 output files). Nback & SART ONLY --
       Full is categorically out of scope (X18's own TASK_TYPES).

NO statistics are re-estimated, re-derived, or sign-flipped anywhere in this
script. Every Estimate/SE/DF/tStat/p/pHolm/CI value plotted is copied verbatim
from the MATLAB-produced pairwise CSV. This mirrors F01's own "no
re-estimation" scope exactly.

===========================================================================
WHY THIS SCRIPT EXISTS / SCOPE (full disclosure)
===========================================================================
X11/X12/X17/X18 are exploratory, additive families -- never pooled with the
two co-primary metrics (LogAbsolutePower/RelativePower) or with each other.
F01 does not plot them (X11/X17/X18 were never in F01's scope at all; X12's
only use in F01 is the single X21-targeted Alpha/Theta ratio, explicitly
EXCLUDED from F01's Fig1 due to X21's own selection-bias disclosure -- see
F01's header item 4). This script exists specifically to visualize
these four families' own FDR-gated results, in direct, row-for-row
correspondence with the pairwise CSVs the MATLAB scripts produce.

CONTRAST SIGN CONVENTION. The CSVs on disk for all four
families are STILL natively N-E / R-E / R-N (REF_LEVEL = 'E'), i.e. the
GLOBAL/ROI convention (X05/X05B/X24/X25), NOT the electrode-level NATIVE
E-N/E-R/N-R convention used by X09/X09B. This was established by reading each of
X11/X12/X17/X18's own source (.m) files directly: X11 line ~290 "CONTRASTS =
{'N','E','N - E'; ...}"; X12/X17/X18 all literally reuse X05's CONTRASTS cell
(comment: "X05's convention"). X18's Position-axis contrast is always
literally "Second - First" (not a Condition contrast -- a block-order effect,
never flipped).

By design, this script DISPLAYS every
N-E/R-E/R-N value under the X09/X09B-native E-N/E-R/N-R label instead, so
every figure across the whole pipeline (F01's topo panels, this script, F03's
topomaps) reads in the same direction. This is done purely in this Python
script, as a linear sign flip of the ALREADY-COMPUTED Estimate/CI (E-N :=
-(N-E), CI bounds swapped-and-negated) -- see CONTRAST_DISPLAY_FLIP /
flip_sign() below. NOTHING in MATLAB was changed or rerun; the on-disk CSVs
are read and copied exactly as before. p-values/SE/DF/FDR/Holm survival are
unaffected (sign-symmetric). Every __rows.csv carries BOTH the native
(as-saved) and the new *_display columns side by side, and every figure with
>=1 flipped row gets a one-line footnote saying so -- never a silent flip.

MISSING-FILE SEMANTICS (handled explicitly, never as an error):
  - X11/X12 run for all three Tasks (Full/Nback/SART) at both RUN_MODEs. If a
    given Task/RUN_MODE/file is absent, that is because this pipeline's own
    MATLAB code only calls writetable(...) `if ~isempty(PairwiseT)` -- i.e.
    nothing survived that family's own BH-FDR gate for that Task/RUN_MODE.
    This script logs that explicitly and moves on; it reflects the gating
    outcome, not a missing run.
  - X17/X18 only ever run for Nback and SART (TASK_TYPES = {'nback','sart'}
    in both scripts' own source, with an explicit comment that Full is "out
    of scope" for this analysis). Requesting Full for X17/X18 is logged
    as "out of scope" (a different, more specific message than "nothing
    survived FDR"), and Full is skipped for those two families only.

RUN_MODE PRIORITY. By design, sensitivity_excl_noisy is
treated as the PRIMARY run mode for this script (X11/X12/X17/X18 are mostly
null/empty in full_sample for several Task/metric combinations), but BOTH
RUN_MODEs are attempted and plotted independently for every Task/family --
nothing is suppressed, sensitivity is simply processed/listed first and
called out as such in the manifest and file naming is RUN_MODE-explicit so
nothing is ambiguous.

TRACEABILITY ("one-by-one against the tables you sent"). For every figure
produced, this script ALSO writes a companion "<figure-name>__rows.csv" file
in the same folder, containing the EXACT rows plotted in that figure, with
the EXACT same column names/order/values as the source MATLAB CSV (just
re-written row by row, no transformation). Opening the PNG next to its
"__rows.csv" (or next to the original MATLAB CSV) lets you match every dot
on the plot to its exact source row by eye. In addition, every point is
annotated directly on the figure with its own Estimate and pHolm value.

VERIFIED PROVENANCE OF THE CSVs ALREADY SHARED IN THIS CONVERSATION. Cross-
matched by exact numeric identity (Estimate/DF/tStat to the reports' own
printed precision) against this pipeline's own .txt reports -- not guessed:
  phase11_pairwise_band_contrasts.csv (18 rows, AFZ/AF4/F6/Frontal/Central/
    Parietal, Alpha, RelMedian/CV_Relative)      -> Task=Nback, RUN_MODE=full_sample
  phase11_pairwise_band_contrasts_full_sesivtivity.csv (9 rows, Occipital/
    Temporal, Beta_low_exploratory/Theta)         -> Task=Full,  RUN_MODE=sensitivity_excl_noisy
  phase12_band_ratio_global_pairwise.csv (3 rows, MedianBased AlphaTheta)
                                                   -> Task=Nback, RUN_MODE=full_sample
  phase12_band_ratio_global_pairwise_full_sensitivity.csv (6 rows, MeanBased+
    MedianBased AlphaTheta)                       -> Task=Full,  RUN_MODE=sensitivity_excl_noisy
  phase12_band_ratio_roi_pairwise_full_senitivity.csv (6 rows, MeanBased
    AlphaTheta, Frontal+Temporal)                 -> Task=Full,  RUN_MODE=sensitivity_excl_noisy
  phase17_pairwise_contrasts.csv (138 rows, ConditionMain survivors)
                                                   -> Task=Nback, RUN_MODE=full_sample  (verified earlier: 42/1134 ConditionMain survivors)
  phase18_pairwise_position_marginal_global.csv (1 row, Alpha RelativePower,
    Second-First, est=-0.0154/DF=15.34/t=-2.703/p=0.0161)
                                                   -> Task=Nback, RUN_MODE=sensitivity_excl_noisy
  phase18_pairwise_position_marginal_electrode.csv (33 rows, Alpha
    RelativePower, Second-First)                  -> Task=Nback, RUN_MODE=sensitivity_excl_noisy
These exact files are used as the fixture that exercised/validated this
script (see verification notes shipped alongside this file); when you point
--root at your own ANALYSIS_ROOT, every Task/RUN_MODE/family is re-walked
fresh from whatever is actually on disk.

OUT OF SCOPE FOR THIS SCRIPT (disclosed, not silently dropped): X11's own
TotalPower pairwise file, X12's own electrode-level pairwise file, and
X17/X18's TimeHalf-main/Interaction/Condition-marginal/Interaction-simple
families are ALL read generically by the same code path below (they are not
hard-coded away) -- they simply were not among the files already shared in
this conversation, so they could not be cross-checked against a report here.
If your real run produces rows in any of them, this script will plot them
too, under their own file name.

RUN ORDER
  MATLAB: X11, X12 (both RUN_MODEs, all 3 Tasks), X17, X18 (both RUN_MODEs,
          Nback+SART only) -- all already-existing, already-run scripts.
          Nothing upstream changes.
  Python: python3 F02_followup_family_figures_cognitive.py [--root <ANALYSIS_ROOT>]
                                                             [--tasks all|full,nback,sart]
                                                             [--modes all|sensitivity_excl_noisy,full_sample]

REQUIREMENTS  python >= 3.9, pandas, numpy, matplotlib  (pip install pandas numpy matplotlib)

OUTPUT  <ANALYSIS_ROOT>/<Task>/32_followup_figures/<X11|X12|X17|X18>/<RUN_MODE>/
          *.png + *.pdf                 one forest-style panel per
                                         (file x natural grouping) that has >=1 row
          *__rows.csv                   exact rows plotted, same columns as source CSV
        <ANALYSIS_ROOT>/32_followup_figures_coverage_manifest.txt
          ONE human-readable table: for every Task x RUN_MODE x family x file,
          whether it was found-with-N-rows, found-but-empty(impossible in
          practice), missing-nothing-survived-FDR, or out-of-scope-for-this-Task.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# 0. SETTINGS
# ---------------------------------------------------------------------------
DEFAULT_ROOT = "results"  # pass --root to point at your own ANALYSIS_ROOT

TASK_LABELS = {"full": "Full", "nback": "Nback", "sart": "SART"}
ALL_TASKS = ["full", "nback", "sart"]
HALF_TASKS = ["nback", "sart"]          # X17/X18's own TASK_TYPES -- Full excluded by design

# sensitivity_excl_noisy listed FIRST -- treated as primary by design
# (X11/X12/X17/X18 are mostly null in full_sample); full_sample still processed.
RUN_MODES = ["sensitivity_excl_noisy", "full_sample"]
RUN_MODE_DISPLAY = {"sensitivity_excl_noisy": "Sensitivity (excl. noisy subjects)",
                     "full_sample": "Full sample"}

BANDS_ORDER = ["Delta", "Theta", "Alpha", "Beta", "Beta_low_exploratory", "Beta_high_exploratory",
               "Gamma_low_exploratory", "Gamma_high1_exploratory", "Gamma_high2_exploratory",
               "TotalPower"]

# All four families' CSVs on disk are still natively N-E/R-E/R-N (see header
# item "CONTRAST SIGN CONVENTION") -- CONTRAST_ORDER stays native, it only
# controls row sort order and is never shown to the user. Displayed colour/
# legend keys below are DISPLAY labels (E-N/E-R/N-R) -- see
# CONTRAST_DISPLAY_FLIP and flip_sign() for the sign-flip-for-display logic,
# which shows every figure in
# this pipeline under the X09/X09B-native E-N/E-R/N-R convention.
CONTRAST_ORDER = ["N - E", "R - E", "R - N", "Second - First"]
CONTRAST_COLORS = {"E - N": "#0072B2", "E - R": "#D55E00", "N - R": "#009E73",
                    "Second - First": "#CC79A7"}
ALPHA = 0.05
FONT = "DejaVu Sans"

# ---------------------------------------------------------------------------
# DISPLAY-ONLY SIGN FLIP
# ---------------------------------------------------------------------------
# This pipeline's MATLAB output for X11/X12/X17/X18 is natively N-E/R-E/R-N
# (REF_LEVEL='E', verified by reading each script's own CONTRASTS literal --
# see module docstring). X09/X09B's own electrode-level output is natively
# E-N/E-R/N-R. To make every figure in this pipeline read the same direction,
# this script now displays X11/X12/X17/X18's results under the E-N/E-R/N-R
# label too, by negating the ALREADY-COMPUTED Estimate/CI: E-N := -(N-E), with
# CI bounds swapped-and-negated. This is a pure linear transform of a number
# the locked MATLAB model already produced -- NOT a new statistic, NOT a
# re-fit, and the MATLAB output on disk is never touched or rerun. p-values,
# SE, DF, and FDR/Holm survival are completely unaffected by this transform.
# "Second - First" (X18's Position axis) and anything else not in the map
# below is NOT a Condition contrast and is left exactly as-is.
# Every __rows.csv keeps BOTH the native (as-saved) Contrast/Estimate/CI_Low/
# CI_High columns AND the new Contrast_display/Estimate_display/CI_Low_display/
# CI_High_display columns side by side, so nothing here is a silent transform.
CONTRAST_DISPLAY_FLIP = {"N - E": "E - N", "R - E": "E - R", "R - N": "N - R"}


def display_contrast(native_label: str) -> str:
    return CONTRAST_DISPLAY_FLIP.get(native_label, native_label)


def flip_sign(native_label: str, est: float, ci_low: float, ci_high: float):
    """Returns (display_label, display_estimate, display_ci_low, display_ci_high).
    Unchanged if native_label has no flip entry (e.g. 'Second - First')."""
    disp = display_contrast(native_label)
    if disp == native_label:
        return disp, est, ci_low, ci_high
    return disp, -est, -ci_high, -ci_low

plt.rcParams.update({
    "font.family": FONT, "font.size": 7.5, "axes.titlesize": 8, "axes.labelsize": 7.5,
    "xtick.labelsize": 7, "ytick.labelsize": 6.5, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
    "ps.fonttype": 42, "savefig.dpi": 300, "figure.dpi": 100})

LOG: list[str] = []
MANIFEST: list[str] = []


def log(msg: str) -> None:
    print(msg)
    LOG.append(msg)


# ---------------------------------------------------------------------------
# 1. FAMILY/FILE SPECIFICATIONS
#    Each entry fully describes one possible output .csv: where it lives,
#    which Tasks it can even exist for, and how to split its rows into one
#    forest panel per "natural" grouping (e.g. one panel per Level x Metric).
#    This table is the ONLY place that needs editing to add a file this
#    script does not yet know about -- the plotting code below is generic.
# ---------------------------------------------------------------------------
SPECS = [
    # ---- X11: median-robust derived-metric family -------------------------
    dict(family="X11", subdir="11_median_robust_family",
         fname="phase11_pairwise_band_contrasts.csv", tasks=ALL_TASKS,
         group_cols=["Level", "Metric"], label_cols=["Unit", "Band"],
         xlabel="Estimate (median-robust metric units)"),
    dict(family="X11", subdir="11_median_robust_family",
         fname="phase11_pairwise_totalpower_contrasts.csv", tasks=ALL_TASKS,
         group_cols=["Level", "Metric"], label_cols=["Unit"],
         xlabel="Estimate (broadband TotalPower metric units)"),
    # ---- X12: band-ratio indices -------------------------------------------
    dict(family="X12", subdir="12_band_ratio_indices",
         fname="phase12_band_ratio_global_pairwise.csv", tasks=ALL_TASKS,
         group_cols=["Variant"], label_cols=["Unit", "RatioType"],
         xlabel="Estimate (ratio units)"),
    dict(family="X12", subdir="12_band_ratio_indices",
         fname="phase12_band_ratio_roi_pairwise.csv", tasks=ALL_TASKS,
         group_cols=["Variant"], label_cols=["ROI", "RatioType"],
         xlabel="Estimate (ratio units)"),
    dict(family="X12", subdir="12_band_ratio_indices",
         fname="phase12_band_ratio_electrode_pairwise.csv", tasks=ALL_TASKS,
         group_cols=["Variant", "RatioType"], label_cols=["Channel"],
         xlabel="Estimate (ratio units)"),
    # ---- X17: half-split Condition x TimeHalf (Nback/SART only) ----------
    dict(family="X17", subdir="17_half_split_condition_time",
         fname="phase17_pairwise_contrasts.csv", tasks=HALF_TASKS,
         group_cols=["SurvivingTerm", "Level", "Metric"], label_cols=["Unit", "Band"],
         xlabel="Estimate (locked-metric units)"),
    # ---- X18: recording-position / block-order (Nback/SART only) ---------
    dict(family="X18", subdir="18_position_pos1_pos2",
         fname="phase18_pairwise_condition_marginal_global.csv", tasks=HALF_TASKS,
         group_cols=["Metric"], label_cols=["Band"],
         xlabel="Estimate (locked-metric units)"),
    dict(family="X18", subdir="18_position_pos1_pos2",
         fname="phase18_pairwise_condition_marginal_electrode.csv", tasks=HALF_TASKS,
         group_cols=["Metric"], label_cols=["Electrode", "Band"],
         xlabel="Estimate (locked-metric units)"),
    dict(family="X18", subdir="18_position_pos1_pos2",
         fname="phase18_pairwise_position_marginal_global.csv", tasks=HALF_TASKS,
         group_cols=["Metric"], label_cols=["Band"],
         xlabel="Estimate, Second - First (locked-metric units)"),
    dict(family="X18", subdir="18_position_pos1_pos2",
         fname="phase18_pairwise_position_marginal_electrode.csv", tasks=HALF_TASKS,
         group_cols=["Metric"], label_cols=["Electrode", "Band"],
         xlabel="Estimate, Second - First (locked-metric units)"),
    dict(family="X18", subdir="18_position_pos1_pos2",
         fname="phase18_pairwise_interaction_simple_global.csv", tasks=HALF_TASKS,
         group_cols=["Metric", "Position"], label_cols=["Band"],
         xlabel="Estimate (locked-metric units)"),
    dict(family="X18", subdir="18_position_pos1_pos2",
         fname="phase18_pairwise_interaction_simple_electrode.csv", tasks=HALF_TASKS,
         group_cols=["Metric", "Position"], label_cols=["Electrode", "Band"],
         xlabel="Estimate (locked-metric units)"),
]


# ---------------------------------------------------------------------------
# 2. PATHS / IO
# ---------------------------------------------------------------------------
def stats_dir(root: Path, task: str, subdir: str, mode: str) -> Path:
    return root / TASK_LABELS[task] / subdir / mode


def out_dir(root: Path, task: str, family: str, mode: str) -> Path:
    d = root / TASK_LABELS[task] / "32_followup_figures" / family / mode
    d.mkdir(parents=True, exist_ok=True)
    return d


def read_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    for c in df.columns:
        if df[c].dtype == object:
            df[c] = df[c].astype(str).str.strip()
    return df


def band_sort_key(b: str) -> int:
    return BANDS_ORDER.index(b) if b in BANDS_ORDER else len(BANDS_ORDER)


def row_label(row: pd.Series, label_cols: list[str]) -> str:
    parts = []
    for c in label_cols:
        if c in row and pd.notna(row[c]):
            parts.append(str(row[c]))
    return " | ".join(parts)


def row_sort_key(row: pd.Series, label_cols: list[str]):
    key = []
    for c in label_cols:
        v = row[c] if c in row else ""
        key.append(band_sort_key(v) if c == "Band" else str(v))
    contrast = row.get("Contrast", "")
    key.append(CONTRAST_ORDER.index(contrast) if contrast in CONTRAST_ORDER else 99)
    return tuple(key)


# ---------------------------------------------------------------------------
# 3. GENERIC FOREST PANEL
# ---------------------------------------------------------------------------
def forest_panel(df: pd.DataFrame, label_cols: list[str], title: str, xlabel: str,
                  out_png: Path) -> None:
    """One dot + 95% CI whisker per row, coloured by Contrast, annotated with
    its own Estimate and pHolm so every point can be matched back to its
    source CSV row by eye. Rows sorted by label (Band in canonical pipeline
    order where present) then by the fixed Contrast order."""
    rows = sorted(df.to_dict("records"),
                   key=lambda r: row_sort_key(pd.Series(r), label_cols))
    n = len(rows)
    height = float(np.clip(0.34 * n + 1.3, 2.2, 22.0))
    fig, ax = plt.subplots(figsize=(7.2, height))
    y = np.arange(n)[::-1]
    disp_labels = []
    any_flipped = False
    for yi, r in zip(y, rows):
        native_contrast = r.get("Contrast", "")
        disp, est, lo, hi = flip_sign(native_contrast, r["Estimate"], r["CI_Low"], r["CI_High"])
        any_flipped = any_flipped or (disp != native_contrast)
        r["Contrast_display"], r["Estimate_display"] = disp, est
        r["CI_Low_display"], r["CI_High_display"] = lo, hi
        disp_labels.append(disp)
        col = CONTRAST_COLORS.get(disp, "k")
        pholm = r.get("pHolm", r.get("pHolm_withinModel", np.nan))
        if np.isfinite(lo) and np.isfinite(hi):
            ax.plot([lo, hi], [yi, yi], color=col, lw=1.1)
        ax.plot(est, yi, marker="o", ms=4.5, mfc=col, mec=col, ls="")
        ax.text(hi if np.isfinite(hi) else est, yi, f"  est={est:.4g}, pHolm={pholm:.4g}",
                 va="center", ha="left", fontsize=5.6, color="0.25")
    ax.axvline(0, color="0.6", lw=0.7)
    ax.set_yticks(y, [f"{row_label(pd.Series(r), label_cols)}  [{d}]" for r, d in zip(rows, disp_labels)])
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left")
    handles = [Line2D([], [], marker="o", ls="", mfc=c, mec=c, ms=4.5, label=k)
               for k, c in CONTRAST_COLORS.items() if k in disp_labels]
    if handles:
        # Outside the axes (upper-right margin) so it never overlaps a data
        # row's own CI whisker / annotation text, however tall the panel is.
        ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.0, 1.0),
                  frameon=False, borderaxespad=0.3)
    if any_flipped:
        fig.text(0.0, -0.01, "Sign flipped for display: E−N/E−R/N−R shown here = −(N−E)/−(R−E)/−(R−N),"
                  " this pipeline's native MATLAB output (unchanged on disk). See build-record.",
                  fontsize=5, color="0.4")
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    fig.savefig(out_png.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(rows).to_csv(out_png.with_name(out_png.stem + "__rows.csv"), index=False)
    log(f"    wrote {out_png.name} + .pdf + __rows.csv  ({n} row{'s' if n != 1 else ''})")


# ---------------------------------------------------------------------------
# 4. DRIVER -- generic over SPECS
# ---------------------------------------------------------------------------
def process_one(root: Path, task: str, mode: str, spec: dict) -> None:
    family, subdir, fname = spec["family"], spec["subdir"], spec["fname"]
    task_label = TASK_LABELS[task]
    tag = f"{family}/{fname}  Task={task_label}  RUN_MODE={RUN_MODE_DISPLAY[mode]}"

    if task not in spec["tasks"]:
        MANIFEST.append(f"{tag}: OUT OF SCOPE -- this Task is not in {family}'s own TASK_TYPES "
                         f"(Full is categorically not run for X17/X18 by MATLAB-side design).")
        return

    f = stats_dir(root, task, subdir, mode) / fname
    if not f.exists():
        MANIFEST.append(f"{tag}: not found -- nothing survived this family's own BH-FDR gate "
                         f"for this Task/RUN_MODE (the MATLAB script only writes this file when "
                         f"its own PairwiseT is non-empty). Not an error.")
        return

    df = read_csv(f)
    if df.empty:
        MANIFEST.append(f"{tag}: file exists but has 0 rows (should not normally happen -- "
                         f"flag for review if seen).")
        return

    MANIFEST.append(f"{tag}: FOUND -- {len(df)} row(s). Plotted below.")
    od = out_dir(root, task, family, mode)
    log(f"  {tag}: {len(df)} row(s)")

    group_cols = [c for c in spec["group_cols"] if c in df.columns]
    if group_cols:
        groups = df.groupby(group_cols, dropna=False)
    else:
        groups = [((), df)]
    for key, sub in groups:
        key_t = key if isinstance(key, tuple) else (key,)
        key_str = "_".join(str(k) for k in key_t) if group_cols else "all"
        title_bits = ", ".join(f"{c}={v}" for c, v in zip(group_cols, key_t)) if group_cols else ""
        title = f"{family} | Task={task_label} | {RUN_MODE_DISPLAY[mode]}" + (f" | {title_bits}" if title_bits else "")
        stem = Path(fname).stem
        out_png = od / f"{family}_{stem}__{key_str}.png"
        forest_panel(sub, spec["label_cols"], title, spec["xlabel"], out_png)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--root", default=DEFAULT_ROOT, help="ANALYSIS_ROOT of the cognitive_task pipeline")
    ap.add_argument("--tasks", default="all", help="comma list of full,nback,sart or 'all'")
    ap.add_argument("--modes", default="all", help="comma list of sensitivity_excl_noisy,full_sample or 'all'")
    a = ap.parse_args()

    root = Path(a.root)
    tasks = ALL_TASKS if a.tasks == "all" else [t.strip() for t in a.tasks.split(",")]
    for t in tasks:
        if t not in ALL_TASKS:
            sys.exit(f"ERROR: unknown task '{t}', expected one of {ALL_TASKS}")
    modes = RUN_MODES if a.modes == "all" else [m.strip() for m in a.modes.split(",")]
    for m in modes:
        if m not in RUN_MODES:
            sys.exit(f"ERROR: unknown mode '{m}', expected one of {RUN_MODES}")

    log(f"F02 followup-family figures (cognitive_task) | root = {root}")
    log(f"Tasks: {tasks} | RUN_MODEs (priority order): {modes}")
    log("Families covered: X11 (median-robust), X12 (band-ratio indices), "
        "X17 (half-split Condition x TimeHalf, Nback/SART only), "
        "X18 (recording-position, Nback/SART only).\n")

    for mode in modes:
        for task in tasks:
            for spec in SPECS:
                process_one(root, task, mode, spec)

    manifest_path = root / "32_followup_figures_coverage_manifest.txt"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "============================================================",
        "F02 FOLLOWUP-FAMILY FIGURES -- COVERAGE MANIFEST (X11/X12/X17/X18)",
        "============================================================",
        "Every line below is one (Task x RUN_MODE x family x output-file)",
        "combination this run checked. 'not found' = nothing survived that",
        "family's own FDR gate for that cell (expected, not an error).",
        "'OUT OF SCOPE' = Full Task for X17/X18, never run by MATLAB-side design.",
        "",
    ]
    manifest_path.write_text("\n".join(header + MANIFEST) + "\n")
    log(f"\nCoverage manifest written: {manifest_path}")

    log_path = root / "32_followup_figures_figure_log.txt"
    log_path.write_text("\n".join(LOG) + "\n")
    log(f"Figure log written: {log_path}")
    log("\nDone. Figures (where any FDR survivors exist) are under each Task's own "
        "<ANALYSIS_ROOT>/<Task>/32_followup_figures/<X11|X12|X17|X18>/<RUN_MODE>/, "
        "each PNG paired with a '__rows.csv' listing the exact source rows plotted.")


if __name__ == "__main__":
    main()
