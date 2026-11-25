# Cognitive-task EEG

EEG analyses of the testing blocks that followed each intervention. Three contexts are analysed: the N-back, the SART, and both tasks together ("Full"). The pre-intervention baseline block (block 1) is not used in these models.

Each statistics script fits `Metric ~ Condition + (1|Subject)` by REML with Satterthwaite degrees of freedom. It corrects the omnibus tests with Benjamini-Hochberg FDR within each family and applies Holm correction to the three pairwise contrasts (E − N, E − R, N − R) of every unit that survives FDR. Absolute and relative power always form separate families.

## Setup

1. Set `DATA_DIR` in `config.m`. Outputs go to `results/` in this folder, in a `<Full|Nback|SART>/<NN_name>/<RUN_MODE>/` structure.
2. Start MATLAB inside this folder. Start EEGLAB before the scripts that read `.set` files.
3. Each statistics script has a `RUN_MODE` line near the top. Run it once with `'full_sample'` (N = 30) and once with `'sensitivity_excl_noisy'` (cleaned sample, N = 28). The excluded participants are listed in `EXCLUSION_CONFIG_cognitive.m`.
4. Run the figure scripts from this folder with Python 3 (`numpy`, `pandas`, `matplotlib`, `mne`). `X15` also needs `fooof`. Use `--root` if the outputs are elsewhere.

## Run order

| Order | Script | What it does | Manuscript item |
|---|---|---|---|
| 1 | `X00_noise_audit_cognitive.m` | Composite noise score per recording and participant | Methods |
| 2 | `X01_audit_cognitive_data.m` | Inventory of the preprocessed recordings | Methods |
| 3 | `X02_build_cognitive_metadata.m` | Recording-level metadata and preprocessing QC values | Methods |
| 4 | `X03_compute_cognitive_psd.m` | 4-s windows with 50% overlap, window rejection, power spectra | Methods |
| 5 | `X04_aggregate_cognitive_subject_level.m` | Band power per recording, electrode and band | input to all statistics |
| 6 | `X07_roi_electrode_definitions_cognitive.m` | Channel table; five predefined regions | Methods |
| 7 | `X24_roi_v2_definitions_cognitive.m` | Four literature-based regions and the planned contrasts | Methods |
| 8 | `F00_export_for_python_cognitive.m` | Exports channel positions and band power for Python | input to F01–F03 |
| 9 | `X05_primary_cognitive_statistics.m` ×2 | Whole-scalp band power | Tables S3, S4; Figs. S4–S6 |
| 10 | `X05B_all_contrasts_ungated_cognitive.m` ×2 | All whole-scalp contrasts without FDR gating (descriptive) | input to figures |
| 11 | `X08_roi_level_statistics_cognitive.m` ×2 | Band power in the five predefined regions | Tables S3, S4 |
| 12 | `X09_electrode_level_statistics_cognitive.m` ×2 | Electrode-level band power | Tables S3, S4; Fig. 5; Figs. S16–S19 |
| 13 | `X09B_electrode_model_estimates_cognitive.m` ×2 | All electrode-level contrasts without gating (descriptive) | input to figures |
| 14 | `X11_median_robust_family_cognitive.m` ×2 | Median-based band power; CV across windows | Tables S3, S4; Fig. S21; Figs. S23–S25 |
| 15 | `X11B_electrode_model_estimates_cognitive.m` ×2 | Ungated electrode-level contrasts for `X11` | input to F03 |
| 16 | `X12_band_ratio_indices_cognitive.m` ×2 | θ/β, α/θ, α/β, engagement index | Tables S3, S4; Figs. S23–S25 |
| 17 | `X12B_electrode_ratio_estimates_cognitive.m` ×2 | Ungated electrode-level contrasts for `X12` | input to F03 |
| 18 | `X13_asymmetry_indices_cognitive.m` ×2 | Alpha asymmetry | Tables S3, S4 |
| 19 | `X14_spectral_shape_features_cognitive.m` | Individual alpha frequency and SEF95 | Tables S3, S4 |
| 20 | `X15_aperiodic_1f_family_cognitive.m` | Aperiodic exponent and offset (FOOOF; calls `X15_fit_fooof.py`) | Tables S3, S4 |
| 21 | `X16_full_reconstructed_cognitive.m` ×2 | Full context rebuilt from the N-back and SART windows | Tables S3, S4 |
| 22 | `X17_half_split_condition_time_cognitive.m` ×2 | Condition × recording-half model | Tables S3, S4; Fig. S20; Figs. S23–S25 |
| 23 | `X17B_electrode_all_contrasts_ungated_cognitive.m` ×2 | Ungated electrode-level contrasts for `X17` | input to F03 |
| 24 | `X18_position_pos1_pos2_cognitive.m` ×2 | Condition × intervention-position model | Tables S3, S4; Fig. S22; Figs. S23–S25 |
| 25 | `X18B_electrode_all_contrasts_ungated_cognitive.m` ×2 | Ungated electrode-level contrasts for `X18` | input to F03 |
| 26 | `X20_individualized_alpha_cognitive.m` ×2 | Power in IAF ± 2 Hz | Tables S3, S4 |
| 27 | `X25_roi_v2_statistics_cognitive.m` ×2 | Planned contrasts C-T1 and C-T2 with order-adjusted, within-person and leave-one-out checks; regional band power | Table 2; Tables S3, S4; Figs. S7–S9 |
| 28 | `F01_publication_figures_cognitive.py` | Main and supplementary figures | Fig. 5; Figs. S4–S19 |
| 29 | `F02_followup_family_figures_cognitive.py` | Forest plots for the follow-up families | Fig. S21; Figs. S23–S25 |
| 30 | `F03_topomap_figures_cognitive.py` | Electrode topographies for the follow-up families | Figs. S20, S22 |

Helper files: `config.m`, `EXCLUSION_CONFIG_cognitive.m`, `X15_fit_fooof.py`.

## Input data

- `DATA_DIR` holds the preprocessed recordings, one `.set`/`.fdt` pair per participant, task and block, named `sub<N>_<code>_<block>_<full|nback|sart>_pp.set` (N = 1–30; block 1 = baseline, not used in the models).
- Each recording has 62 scalp channels at 512 Hz, preprocessed as described in the Methods.
- The preprocessing QC files are read by `X02`.
- `channel_locations_template.mat` (in `DATA_DIR`) is required by `X07`.

## Notes

- No analysis step uses random numbers, so no seed is needed.
- The slowest steps are `X03` and `X15`.
- In some output tables the native contrast label is N − E, R − E or R − N, the opposite sign of the paper. The figure scripts flip the sign for display, so the figures use the paper's convention.
