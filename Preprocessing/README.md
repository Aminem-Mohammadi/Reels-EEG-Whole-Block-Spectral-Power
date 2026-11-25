# Preprocessing

EEG preprocessing in MATLAB with EEGLAB, adapted from the HAPPE pipeline (Gabard-Durnam et al., 2018), with ICLabel classification and visual inspection in place of wavelet-enhanced ICA. The steps and parameters are described in the Methods of the article: re-reference to the average of A1 and A2, 1-100 Hz FIR band-pass, 47-53 Hz notch, removal of the first and last 5 s, bad-channel detection, extended ICA on a 1-45 Hz copy, ICLabel-based component removal, and spherical interpolation.

The outputs of this folder (`*_pp.set` / `.fdt`) are the inputs of `Intervention/` and `Cognitive_Tasks/EEG/`.

[Add the preprocessing scripts and their run order here.]
