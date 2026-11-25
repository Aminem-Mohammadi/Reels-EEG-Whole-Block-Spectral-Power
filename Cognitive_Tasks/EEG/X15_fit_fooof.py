#!/usr/bin/env python3
"""
================================================================
X15_fit_fooof.py

Standalone Python companion script for
X15_aperiodic_1f_family_cognitive.m (cognitive-task EEG pipeline,
exploratory-extension family X11-X15, purely additive -- see that
script's header for full context; this file is never imported by
anything else in the pipeline, it is invoked as a plain OS-level
subprocess).

PURPOSE
  Batch-fits the FOOOF ("Fitting Oscillations & One Over F") model to
  many individual power spectra -- one row per Subject x TaskType x
  BlockNumber x Condition x Channel, restricted to the CONTINUOUS
  1-45 Hz retained band (never the 55-100 Hz portion, which is
  separated from 1-45 Hz by a 45-55 Hz notch guard band and would
  break FOOOF's continuous-spectrum assumption) -- and writes back
  ONLY the APERIODIC parameters (offset, exponent) plus a
  goodness-of-fit (R^2) and a per-row fit status.

  This is the PRIMARY method for aperiodic-exponent/offset estimation
  in X15 because it performs a proper periodic/aperiodic decomposition
  (it fits and removes putative oscillatory peaks before estimating
  the 1/f background), which a plain log-log linear fit cannot do --
  a strong alpha or beta peak can otherwise bias a naive linear fit's
  slope estimate. MATLAB's own robustfit()/polyfit() log-log linear
  fit (computed independently, entirely inside the .m script, no
  round trip through this file) is the MATLAB-native sensitivity /
  fallback check, in the same methodological spirit as the
  gamma-robustness jackknife already used elsewhere in this pipeline
  (X06_gamma_robustness_cognitive_check.m) -- two independent methods,
  neither privileged as "ground truth", compared for agreement.

INSTALLATION
    pip install fooof

  NOTE ON PACKAGE NAMING: 'fooof' is the PyPI package name current as
  of when this script was written. Its author-designated
  successor/rename is 'specparam' (same maintainers, same underlying
  algorithm, newer name/API). This script targets 'fooof' as primary
  because that is the package installable
  on the machine that will run this. If 'pip install fooof' ever stops
  resolving, the documented migration is: `pip install specparam` and
  change the `from fooof import FOOOF` line below to
  `from specparam import SpectralModel as FOOOF` (specparam's fit
  API is a drop-in match for aperiodic_params_/r_squared_). This
  script does NOT attempt to auto-detect or silently fall back to
  specparam -- that substitution is a disclosed manual step, not an
  automatic one.

WHY 'fixed' APERIODIC MODE (not 'knee')
  FOOOF's aperiodic component can be fit in two modes:
    - 'fixed' : log10(power) = offset - exponent * log10(freq)
                a single straight line in log-log space.
    - 'knee'  : adds a bend/knee parameter for spectra that flatten
                out at low frequencies (common in intracranial/LFP
                recordings, or very broad ranges that extend down
                toward 0 Hz).
  This script's fit range is 1-45 Hz: a comparatively narrow,
  broadband scalp-EEG range, with no a priori physiological or visual
  evidence of a low-frequency knee/bend within that range (per this
  pipeline's own diagnostic PSD plots, X03's Diagnostic 1). Adding an
  unneeded knee parameter would only inflate the model's degrees of
  freedom and destabilize the exponent/offset estimates on individual
  (not grand-averaged) recording spectra, which are already noisier
  than a grand average. 'fixed' aperiodic_mode is therefore the
  standard, recommended default for a range like this and is used
  uniformly for every row, every Task. Testing whether a knee
  parameter improves fit is a natural follow-on analysis but is
  explicitly OUT OF SCOPE for this exploratory extension.

USAGE (invoked by X15_aperiodic_1f_family_cognitive.m via MATLAB's
plain system() call -- NOT MATLAB's py.* interop, so that the
round-trip stays a simple OS-level call MATLAB can check the exit
status of):

    python3 X15_fit_fooof.py \
        --input  <fooof_input_TaskLabel_VariantName.csv> \
        --output <fooof_output_TaskLabel_VariantName.csv> \
        --freq-range 1 45

INPUT CSV FORMAT (written by the MATLAB script; WIDE format -- one row
per Subject x TaskType x BlockNumber x Condition x Channel spectrum;
see the MATLAB script's header for why wide, not long, was chosen):
  Identifier columns : Subject, TaskType, BlockNumber, Condition,
                        Channel, RowID. RowID is a unique integer key
                        used ONLY so the MATLAB script can merge this
                        script's output back onto its own table
                        row-for-row; it carries no scientific meaning
                        and must be echoed back unchanged.
  Frequency columns  : one column per frequency bin, named exactly
                        "F_<freqHz>" with 4 decimal places, e.g.
                        "F_1.0000", "F_1.2500", ..., "F_45.0000".
                        Every column matching this pattern is treated
                        as a power value (uV^2/Hz, MATLAB pwelch()
                        units, unchanged from X03) at that frequency.
                        Columns are re-sorted by their embedded
                        frequency value before fitting -- their order
                        in the CSV is never assumed to already be
                        ascending.

OUTPUT CSV FORMAT: the same identifier columns (Subject, TaskType,
BlockNumber, Condition, Channel, RowID), in the same row order as the
input, plus four new columns:
    AperiodicExponent_FOOOF, AperiodicOffset_FOOOF, FOOOF_RSquared,
    FitStatus   ('OK' or 'FAILED')
A row whose fit fails (any exception during FOOOF fitting, or
non-finite/non-positive input power) gets FitStatus='FAILED' and NaN
for the three numeric columns -- it is NEVER silently dropped, so
output.csv always has exactly the same number of rows as input.csv.
This lets MATLAB count and log FOOOF fit failures per Task without
guessing which rows are missing.
================================================================
"""

import argparse
import re
import sys
import time

import numpy as np
import pandas as pd

try:
    from fooof import FOOOF
except ImportError:
    sys.stderr.write(
        "ERROR: the 'fooof' package is not installed in this Python "
        "environment.\n"
        "  Fix with:  pip install fooof\n"
        "  (If 'fooof' is ever removed from PyPI, its documented "
        "successor is 'specparam' -- pip install specparam, then "
        "edit the import at the top of this script; see this "
        "script's header for the exact one-line change.)\n"
    )
    sys.exit(2)

# Matches column headers of the exact form "F_<digits>.<digits>",
# e.g. "F_1.0000", "F_45.0000". Group 1 captures the frequency value.
FREQ_COL_PATTERN = re.compile(r'^F_([0-9]+\.[0-9]+)$')

# Standard, disclosed default FOOOF settings for this 1-45 Hz broadband
# scalp-EEG range (see header for the 'fixed' mode rationale).
FOOOF_PEAK_WIDTH_LIMITS = [1, 8]
FOOOF_MAX_N_PEAKS = 6
FOOOF_APERIODIC_MODE = 'fixed'


def parse_args():
    p = argparse.ArgumentParser(
        description='Batch FOOOF aperiodic (offset/exponent) fitting for X15 '
                    '(cognitive-task EEG pipeline, exploratory extension).')
    p.add_argument('--input', required=True,
                    help='Input CSV path (wide format, written by the MATLAB script).')
    p.add_argument('--output', required=True,
                    help='Output CSV path to write (same row order/count as input).')
    p.add_argument('--freq-range', nargs=2, type=float, required=True,
                    metavar=('FMIN', 'FMAX'),
                    help='Frequency range in Hz to restrict fitting to, e.g. --freq-range 1 45. '
                         'Frequency columns outside this range (if any are present in the CSV) '
                         'are ignored, not deleted.')
    p.add_argument('--progress-every', type=int, default=500,
                    help='Print a progress line every N rows (default 500), so a long batch '
                         'has visible stdout while MATLAB''s system() call is blocked on it.')
    return p.parse_args()


def main():
    args = parse_args()
    freq_range = [min(args.freq_range), max(args.freq_range)]

    print(f'[X15_fit_fooof] fooof/aperiodic_mode={FOOOF_APERIODIC_MODE}, '
          f'peak_width_limits={FOOOF_PEAK_WIDTH_LIMITS}, max_n_peaks={FOOOF_MAX_N_PEAKS}, '
          f'freq_range={freq_range}', flush=True)
    print(f'[X15_fit_fooof] Reading input CSV: {args.input}', flush=True)
    df = pd.read_csv(args.input)
    n = len(df)
    print(f'[X15_fit_fooof] Loaded {n} row(s), {len(df.columns)} column(s).', flush=True)

    id_cols = [c for c in df.columns if not FREQ_COL_PATTERN.match(c)]
    freq_cols_all = [c for c in df.columns if FREQ_COL_PATTERN.match(c)]
    if not freq_cols_all:
        sys.stderr.write(
            "ERROR: no frequency columns matching 'F_<freqHz>' (e.g. 'F_1.0000') found in "
            f"{args.input}. Columns present: {list(df.columns)}\n")
        sys.exit(3)

    # Sort frequency columns by their embedded numeric value -- the CSV's
    # on-disk column order is never assumed to already be ascending.
    freq_vals_all = np.array([float(FREQ_COL_PATTERN.match(c).group(1)) for c in freq_cols_all])
    order = np.argsort(freq_vals_all)
    freq_cols_all = [freq_cols_all[i] for i in order]
    freq_vals_all = freq_vals_all[order]

    in_range = (freq_vals_all >= freq_range[0]) & (freq_vals_all <= freq_range[1])
    if not np.any(in_range):
        sys.stderr.write(
            f"ERROR: none of the {len(freq_cols_all)} frequency column(s) in {args.input} "
            f"fall within --freq-range {freq_range}.\n")
        sys.exit(4)
    freq_cols = [c for c, keep in zip(freq_cols_all, in_range) if keep]
    freq_vals = freq_vals_all[in_range]

    print(f'[X15_fit_fooof] Fitting {n} row(s) over {len(freq_cols)} frequency bin(s) '
          f'({freq_vals.min():.4f}-{freq_vals.max():.4f} Hz).', flush=True)

    power_block = df[freq_cols].to_numpy(dtype=float)  # n x nFreq

    exponent = np.full(n, np.nan)
    offset = np.full(n, np.nan)
    r2 = np.full(n, np.nan)
    fit_status = np.array(['FAILED'] * n, dtype=object)

    n_failed = 0
    n_printed_failures = 0
    t0 = time.time()
    for i in range(n):
        power_row = power_block[i, :]
        try:
            if np.any(~np.isfinite(power_row)) or np.any(power_row <= 0):
                raise ValueError('non-finite or non-positive power value(s) in this row')
            fm = FOOOF(peak_width_limits=FOOOF_PEAK_WIDTH_LIMITS,
                       max_n_peaks=FOOOF_MAX_N_PEAKS,
                       aperiodic_mode=FOOOF_APERIODIC_MODE,
                       verbose=False)
            fm.fit(freq_vals, power_row, freq_range)
            ap = fm.aperiodic_params_  # ['fixed' mode -> [offset, exponent]
            offset[i] = ap[0]
            exponent[i] = ap[1]
            r2[i] = fm.r_squared_
            fit_status[i] = 'OK'
        except Exception as exc:  # noqa: BLE001 -- deliberately broad: a single
            # bad row (degenerate spectrum, FOOOF internal convergence
            # failure, etc.) must never crash the whole batch; per spec,
            # mark it and move on.
            n_failed += 1
            fit_status[i] = 'FAILED'
            if n_printed_failures < 20:
                print(f'[X15_fit_fooof] Row {i} (RowID='
                      f'{df["RowID"].iloc[i] if "RowID" in df.columns else "?"}) FAILED: {exc}',
                      flush=True)
                n_printed_failures += 1
            elif n_printed_failures == 20:
                print('[X15_fit_fooof] (further per-row failure messages suppressed; '
                      'final failure count reported at the end)', flush=True)
                n_printed_failures += 1

        if (i + 1) % args.progress_every == 0 or (i + 1) == n:
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else float('nan')
            print(f'[X15_fit_fooof] {i + 1}/{n} rows fit ({n_failed} failed so far), '
                  f'{elapsed:.1f}s elapsed ({rate:.1f} rows/s).', flush=True)

    out = df[id_cols].copy()
    out['AperiodicExponent_FOOOF'] = exponent
    out['AperiodicOffset_FOOOF'] = offset
    out['FOOOF_RSquared'] = r2
    out['FitStatus'] = fit_status

    out.to_csv(args.output, index=False)
    n_ok = n - n_failed
    print(f'[X15_fit_fooof] DONE. {n_ok}/{n} row(s) fit OK, {n_failed}/{n} row(s) FAILED. '
          f'Wrote: {args.output}', flush=True)


if __name__ == '__main__':
    main()
