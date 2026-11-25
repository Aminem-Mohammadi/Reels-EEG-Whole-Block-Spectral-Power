%% ========================================================================
% EXCLUSION_CONFIG_cognitive.m
%
% Shared, single-source-of-truth exclusion list for this analysis folder.
% Every statistics script in this folder calls this as a script (not a
% function -- it must run with `run('EXCLUSION_CONFIG_cognitive.m')` or be
% on the MATLAB path) near its own Section 0, then uses the
% EXCLUDE_NOISY_SUBJECTS vector it defines below to build BOTH a
% full-sample run and a sensitivity run (subjects in the vector removed)
% in the SAME execution -- see any script's Section 0 for exactly how it
% is consumed.
%
% TO ADD OR REMOVE A SUBJECT: edit ONLY the vector below. Nothing else in
% this file, or in any script that reads it, needs to change. Every
% subject in the vector should have a one-line rationale comment here --
% never add a number without saying why, so the exclusion list stays
% auditable on its own, independent of any script's output.
%
% CURRENT LIST AND RATIONALE:
%   7  -- X00 noise-audit composite z-score = 1.75 (elevated on M2
%         high-frequency-ratio AND M3 broadband-power axes). Below the
%         original pre-specified z>2.0 threshold; included after lowering
%         the monitoring threshold to ~1.75 to also catch this subject
%         and Subject 23 (both independently elevated, not cherry-picked
%         by name) -- see noise_audit_report.txt.
%   23 -- X00 noise-audit composite z-score = 1.85 (elevated on M1
%         window-rejection-rate, 37.4% vs ~5% dataset median). Same
%         threshold-lowering rationale as Subject 7.
%
% MONITORING CEILING: no more than ~3-4 subjects should be excluded in
% total. Current count = 2 (7, 23) -- within the ceiling.
% Any further additions should be reviewed against this ceiling before
% being added below.
% ========================================================================

EXCLUDE_NOISY_SUBJECTS = [7, 23];

EXCLUSION_RATIONALE_cognitive = struct( ...
    's7',  'X00 noise-audit composite z=1.75 (M2+M3 elevated); included via threshold lowering to ~1.75.', ...
    's23', 'X00 noise-audit composite z=1.85 (M1 window-rejection-rate elevated, 37.4%); included via the same threshold lowering as Subject 7.');

RUN_MODES_cognitive = {'full_sample', 'sensitivity_excl_noisy'};
RUN_MODES_DISPLAY_cognitive = struct('full_sample', 'FullSample', 'sensitivity_excl_noisy', 'SensitivityExclNoisy');
