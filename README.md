# Instagram Reels, EEG, attention and working memory

Analysis code for the study:

> Ramezani, Z., Mohammadi, A. M., Hajiha, Z. & Dolatshahi, B. *Brain rhythms during Instagram Reels viewing and subsequent performance on sustained attention and working memory: an EEG study.* Manuscript submitted to *Scientific Reports* (2026).

The data are available on the Open Science Framework: **[OSF link / DOI]**.

## Study in brief

Thirty master's students each completed two of three 10-min interventions: Instagram Reels on their own account (E), a silent video of neutral images (N), or eyes-closed rest (R). Each participant did three testing blocks, before and after each intervention. Each block contained a Sustained Attention to Response Task (SART), a visual N-back task (0- to 3-back), and a fatigue rating. EEG (62 scalp channels, 512 Hz) was recorded throughout. EEG measures and change scores of task performance were compared between conditions with linear mixed models. Every EEG analysis was run in the full sample (N = 30) and in a cleaned sample (N = 28) that excludes the two participants with the highest EEG noise scores.

## Repository structure

```
Preprocessing/            EEG preprocessing (EEGLAB, HAPPE-based pipeline with ICLabel)
Task/                     SART and N-back presentation code (MATLAB, Psychtoolbox)
Intervention/             EEG during the three interventions
Cognitive_Tasks/
    EEG/                  EEG during the SART, the N-back, and both tasks together ("Full")
    Behavioral/           SART and N-back performance
```

Each folder has its own README with the run order, the inputs and outputs, and the manuscript item that each script produces.

```
raw EEG ──► Preprocessing ──► Intervention          (EEG during the interventions)
                         └──► Cognitive_Tasks/EEG  (EEG during the tasks)
task log files ──────────────► Cognitive_Tasks/Behavioral
```

## Requirements

| Component | Version | Used in |
|---|---|---|
| MATLAB | 25.2.0.2998904 (R2025b) | all `.m` files |
| Statistics and Machine Learning Toolbox | 25.2 | all statistics scripts (`fitlme`) |
| Signal Processing Toolbox | 25.2 | `Intervention/X03` (`pwelch`) |
| EEGLAB (with ICLabel) | 2025.1.0 (ICLabel 1.7) | Preprocessing; reading `.set` files |
| Psychtoolbox | Psychtoolbox - 3.0.19 | Task |
| Python | ≥ 3.9 | figure scripts (`F01`–`F03`), FOOOF fits |
| Python packages | see `Intervention/requirements_figures.txt`, plus `fooof` | |

## How to run

1. Download the data from OSF.
2. In each analysis folder, open `config.m` and set `DATA_DIR` to the matching data folder. This is the only edit needed. Outputs are written to a `results/` folder inside the analysis folder.
3. Start MATLAB from inside the analysis folder and run the scripts in the order given in that folder's README.
4. EEG statistics scripts have a `RUN_MODE` setting near the top. Run each of them twice: once with `'full_sample'` (N = 30) and once with `'sensitivity_excl_noisy'` (cleaned sample, N = 28).

## Where each result comes from

| Manuscript item | Folder and script |
|---|---|
| Sample and data quality (window retention, noise scores, cleaned sample) | `Intervention/X03`, `X00` |
| Fig. 2a, 2c; Supplementary Table S1; behavioral results | `Cognitive_Tasks/Behavioral/nback_stats.m`, `sart_stats.m` |
| Fig. 2b; fatigue results | **[fatigue analysis: to be added]** |
| Fig. 3, Fig. 4 | `Intervention/F01_publication_figures.py` |
| Fig. 5 | `Cognitive_Tasks/EEG/F01_publication_figures_cognitive.py` |
| Table 2 (planned contrasts) | `Intervention/X25` (C-I1, C-I2); `Cognitive_Tasks/EEG/X25` (C-T1, C-T2) |
| Table 3 | `Intervention/X05`, `X06`, `X09`, `X12`, `X14`, `X15`, `X20`, `X25` |
| Supplementary Table S2 | `Intervention/` (all statistics scripts) |
| Supplementary Tables S3, S4 | `Cognitive_Tasks/EEG/` (all statistics scripts) |
| Supplementary Note | `Intervention/X12` |
| Supplementary Figs. S1–S3 | `Intervention/F01_publication_figures.py` |
| Supplementary Figs. S4–S19 | `Cognitive_Tasks/EEG/F01_publication_figures_cognitive.py` |
| Supplementary Figs. S20, S22 | `Cognitive_Tasks/EEG/F03_topomap_figures_cognitive.py` (from `X17`, `X18`) |
| Supplementary Figs. S21, S23–S25 | `Cognitive_Tasks/EEG/F02_followup_family_figures_cognitive.py` (from `X11`, `X12`, `X17`, `X18`) |

Figure titles and labels written by the code differ from the published figures, which were relabelled and assembled separately. The main differences are:

| Code | Paper |
|---|---|
| `sensitivity_excl_noisy`, "Sensitivity (excl. 7, 23)" | cleaned sample (N = 28) |
| `full_sample` | full sample (N = 30) |
| H1, H2, H-I1, H-I2, H-T1, H-T2 | planned contrasts C-I1, C-I2, C-T1, C-T2 |
| `i` (in file names), E | Instagram Reels |
| `Nback`, `Full` | N-back, Full context |
| X11, X12, X17, X18 (family codes) | median-based, band-ratio, recording-half and intervention-position families |
| `Beta_low_exploratory`, `Beta_high_exploratory`, `Gamma_low_exploratory`, `Gamma_high1_exploratory`, `Gamma_high2_exploratory` | 13–20, 20–30, 30–45, 55–80 and 80–100 Hz |
| `RelativePower`, `LogAbsolutePower` | relative power (share of 1–45 Hz power), absolute power (ln µV²) |
| N − E, R − E, R − N (native sign in some output tables) | E − N, E − R, N − R (figures already use the paper's sign) |

## License

The code is released under the MIT License (see `LICENSE`).

## Citation

If you use this code, please cite the article above. Citation details will be updated on publication.

## Contact

- **Code, analysis methods, and EEG:** Amin M. Mohammadi, School of Electrical and Computer Engineering, College of Engineering, University of Tehran, Tehran, Iran (aminem.mohamadi@gmail.com)

- **Data availability, study design, and cognitive and social-science questions:** Zahra Ramezani, Department of Clinical Psychology, University of Social Welfare and Rehabilitation Sciences, Tehran, Iran (zahra.ramezani.2079@gmail.com)

- **The research group and its work:** Behrooz Dolatshahi, Department of Psychology, University of Social Welfare and Rehabilitation Sciences, Tehran, Iran (be.dolatshahi@uswr.ac.ir)
