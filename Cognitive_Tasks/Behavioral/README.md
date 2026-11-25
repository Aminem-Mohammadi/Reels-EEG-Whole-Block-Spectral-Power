# Behavioral analyses: SART and N-back

Scores the SART and N-back task log files and runs the behavioral statistics.

The scoring scripts read the raw task files, one `.mat` file per participant and testing block. They compute the task indices for every block and form change scores: the post-intervention block minus the block immediately before that intervention. The order of the two tasks in each block (SART first, SN; N-back first, NS) is taken from the start times in the file names. The statistics scripts then fit linear mixed models on the change scores:

```
change ~ −1 + E + N + R + Position + OrderPre + OrderPost + (1 | Subject)
```

The models are fitted with `fitlme` and Satterthwaite degrees of freedom. The scripts report the omnibus tests (FDR-corrected within the primary and within the secondary indices), the pairwise contrasts (Holm-corrected) with Cohen's d, the adjusted means, and the covariate effects.

## Setup

Set `DATA_DIR` in `config.m` and run the scripts from this folder. Outputs go to `results/`.

## Run order

Steps 1 and 3 are independent; step 2 needs step 1, and step 4 needs step 3.

| Order | Script | What it does | Manuscript item |
|---|---|---|---|
| 1 | `sart_analysis.m` | Scores every SART block (commission and omission rates, d′, c, Go RT indices, commission RT, ex-Gaussian parameters, pre-error speeding, post-error slowing) and forms the change scores | Methods (SART) |
| 2 | `sart_stats.m` | Change-score models for the SART indices | Fig. 2c; SART results; Supplementary Table S1 |
| 3 | `nback_analysis_full.m` | Scores every N-back block (hit and false-alarm rates, accuracy, d′, c, hit RT indices, working-memory composites over 1- to 3-back, load slopes) and forms the change scores | Methods (N-back) |
| 4 | `nback_stats.m` | Change-score models for the N-back indices (section A); load effect on raw scores (section D) | Fig. 2a; N-back results; Supplementary Table S1 |

Sections B and C of `nback_stats.m` (per-load and time-course models), the time-on-task and font-size breakdowns in `sart_analysis.m`, and the additional hit-RT indices are descriptive and are not reported in the manuscript.

## Input data

```
DATA_DIR/
    sart_behavioral/      one file per participant and block
    n_back_behavioral/    one file per participant and block
```

- The participant number, block and start time are read from the file names:
  - SART: `S<N>_SART-data-block<B>-…-<Ddd-Mmm-D-YYYY-H-M-S>.mat`
  - N-back: `s<N>_Nback-data-pre<B>-…-<Ddd-Mmm-D-YYYY-H-M-S>.mat`
- Block 1 is the baseline, block 2 follows the first intervention, and block 3 follows the second.
- Each file holds one variable, `data`:
  - **SART:** fields `trial`, `digit`, `fontSize`, `RT` (s from digit onset) and `response`; 225 trials; digit 3 = No-Go.
  - **N-back:** a 1 × 240 struct array with fields `block`, `nBack`, `trial`, `digit`, `isTarget`, `responded`, `responseTime`, `stimOnset` and `correct`.
- The assignment of conditions and task orders to participants is written in the scripts as the `design` table (also given in Supplementary Table S5).

## Notes

- `sart_analysis.m` and `nback_analysis_full.m` need base MATLAB only. `sart_stats.m` and `nback_stats.m` need the Statistics and Machine Learning Toolbox.
- The ex-Gaussian fit (`fminsearch`) is deterministic. The only random numbers are the horizontal jitter of points in the figures, which does not affect any statistic.
- Commission and omission rates are stored as proportions (0–1); the paper reports percentages.
