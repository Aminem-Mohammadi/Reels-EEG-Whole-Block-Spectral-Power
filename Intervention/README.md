# Intervention EEG

EEG analyses of the three 10-min interventions: Instagram Reels (E), neutral images (N), and eyes-closed rest (R).

Starting from the preprocessed EEGLAB recordings, the code:

1. inventories the recordings and builds a metadata table (`X01`, `X02`);
2. splits each recording into 4-s windows with 50% overlap, rejects artefactual windows, and computes Welch power spectra (`X03`);
3. derives absolute, ln-absolute and relative band power for each recording, electrode and band (`X04`);
4. computes the participant-level composite noise score that defines the cleaned sample (`X00`);
5. fits linear mixed models, `Metric ~ Condition + (1|Subject)`, with Benjamini-Hochberg FDR correction of the omnibus tests and Holm correction of the three pairwise contrasts, for every analysis family (`X05`–`X25`);
6. exports the inputs for the figures and draws them (`F00`, `F01`).

## Setup

1. Set `DATA_DIR` in `config.m`. Outputs go to `results/` in this folder.
2. Start MATLAB inside this folder, so that `config.m` and the helper functions are on the path.
3. Start EEGLAB once before `X01`–`X03`, which use `pop_loadset`.
4. Each statistics script has a line `RUN_MODE = 'full_sample';` near the top. Run the script, change the line to `'sensitivity_excl_noisy'`, and run it again. The excluded participants are listed in `EXCLUSION_CONFIG_intervention.m`.
5. Install the Python packages in `requirements_figures.txt` and `fooof`. `X15` calls Python through `system()`.

## Run order

"×2" means the script is run once for each `RUN_MODE`.

| Order | Script | What it does | Manuscript item |
|---|---|---|---|
| 1 | `X01_audit_intervention_data.m` | Inventory of the recordings, header checks, condition-pair coverage | Methods |
| 2 | `X02_build_intervention_metadata.m` | Recording-level metadata and preprocessing QC values | Methods |
| 3 | `X03_compute_intervention_psd.m` | 4-s windows, window rejection, Welch spectra | Sample and data quality |
| 4 | `X04_aggregate_subject_level.m` | Band power per recording, electrode and band | input to all statistics |
| 5 | `X00_noise_audit_intervention.m` | Composite noise score per participant | Sample and data quality; cleaned sample |
| 6 | `X07_roi_electrode_definitions.m` | Channel table; five predefined regions | Methods |
| 7 | `X05_primary_statistics.m` ×2 | Whole-scalp band power; order-adjusted and within-person checks | Table 3; Table S2; Fig. S1 |
| 8 | `X06_gamma_robustness_check.m` ×2 | Power above 30 Hz without recordings flagged for high-frequency artefact | Results; Table 3 |
| 9 | `X08_roi_level_statistics.m` ×2 | Band power in the five predefined regions | Table 3 |
| 10 | `X09_electrode_level_statistics.m` ×2 | Electrode-level band power (FDR across 62 electrodes) | Table 3; Fig. 3; Fig. S3; Table S2 |
| 11 | `X24_roi_v2_definitions.m` | Four literature-based regions and the planned contrasts | Methods |
| 12 | `X25_roi_v2_statistics.m` ×2 | Planned contrasts C-I1 and C-I2 with order-adjusted, within-person and leave-one-out checks; regional band power | Table 2; Table 3; Fig. S2 |
| 13 | `X11_median_robust_family_intervention.m` ×2 | Median-based band power; CV across windows | Table S2 |
| 14 | `X12_band_ratio_indices_intervention.m` ×2 | θ/β, α/θ, α/β, engagement index; mean-based, median-based and ln variants | Table 3; Table S2; Supplementary Note; Fig. 3 |
| 15 | `X13_asymmetry_indices_intervention.m` ×2 | Alpha asymmetry, ln(right) − ln(left) | Table S2 |
| 16 | `X14_spectral_shape_features_intervention.m` ×2 | Individual alpha frequency and SEF95 | Table 3; Table S2 |
| 17 | `X15_aperiodic_1f_family_intervention.m` ×2 | Aperiodic exponent and offset (FOOOF; calls `X15_fit_fooof.py`) | Table 3; Table S2 |
| 18 | `X17_half_split_condition_time_intervention.m` ×2 | Condition × recording-half model | Table S2 |
| 19 | `X18_position_int1_int2_intervention.m` ×2 | Condition × intervention-position model | Table S2 |
| 20 | `X20_individualized_alpha_intervention.m` ×2 | Power in IAF ± 2 Hz | Table 3; Table S2 |
| 21 | `F00_export_for_python.m` | Exports channel positions, band power and spectra for Python | input to F01 |
| 22 | `F01_publication_figures.py` | Draws the figures (no statistics) | see below |

Helper files: `config.m`, `EXCLUSION_CONFIG_intervention.m`, `apply_run_mode_filter.m`, `expected_condition_counts.m`, `intervention_load.m`, `X15_fit_fooof.py`.

### Figures drawn by `F01_publication_figures.py`

| Output in `results/31_publication_figures/` | Manuscript |
|---|---|
| `main/Fig1_key_EminusN_posterior_theta` | Fig. 3 |
| `main/Fig3_posterior_spectra` | Fig. 4 |
| `main/Fig4_forest_global_EminusN` | Supplementary Fig. S1 |
| `main/Fig5_forest_ROIv2_EminusN` | Supplementary Fig. S2 |
| `topo_panels/topo_<band>_RelativePower_<mode>` (theta, alpha, 30–45 Hz) | Supplementary Fig. S3 |

Other outputs of `F01` are not used in the manuscript.

## Input data

- `DATA_DIR` holds the preprocessed recordings, searched recursively for `*.set` with the matching `.fdt`.
- File names follow `sub<N>_in<K>_<COND>_pp.set`. N is the participant number (1–30), K is the intervention (1 or 2), and COND is `i` (Reels), `n` (neutral images) or `r` (rest), in either case. Example: `sub2_in1_i_pp.set`.
- Each recording has 62 scalp channels at 512 Hz, referenced to the average of A1 and A2, with `boundary` events where data were removed.
- `Preprocessing_QC_<K>.xlsx` (optional) holds the preprocessing QC values used by `X00`.
- `channel_locations_template.mat` (in `DATA_DIR`) holds the channel locations and is required by `X07`.

## Notes

- No analysis step uses random numbers, so no seed is needed.
- Approximate run time on a laptop: `X01`–`X04` 10–30 min; most statistics scripts 2–20 min per `RUN_MODE`; `X15`, `X17` and `X18` 15–40 min per `RUN_MODE`.
- Text reports in `results/09_reports/` refer to the steps as phases: Phase 0 = `X01`, 1 = `X02`, 2 = `X03`, 3 = `X04`, 4 = `X05`/`X06`, 5D = `X07`, 5A = `X08`, 5B = `X09`, 6a = `X24`, 6b = `X25`.
- The literature-based regions (`X24`) were defined after the electrode-level results had been inspected, as stated in the manuscript.
