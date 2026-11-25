#!/usr/bin/env python3
"""
================================================================
X15_fit_fooof.py

Helper for X15_aperiodic_1f_family_intervention.m, called by it as a subprocess
(python3 X15_fit_fooof.py --input <csv> --output <csv> --freq-range 1 45).

Fits the FOOOF model to each power spectrum in the input CSV (one row per
recording x channel, columns F_<frequency>) over the continuous 1-45 Hz range and
writes the aperiodic parameters (offset, exponent), the goodness of fit (R^2) and a
per-row fit status. Settings: aperiodic_mode = 'fixed', peak_width_limits = [1, 8],
max_n_peaks = 6. The MATLAB script additionally computes a robust log-log linear
fit for comparison.

INSTALLATION
    pip install fooof
  ('specparam' is the renamed successor package; to use it, replace
   `from fooof import FOOOF` with `from specparam import SpectralModel as FOOOF`.)
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

# FOOOF settings for the 1-45 Hz scalp-EEG range ('fixed' aperiodic mode).
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
