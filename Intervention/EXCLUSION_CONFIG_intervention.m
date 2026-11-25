%% ========================================================================
% EXCLUSION_CONFIG_intervention.m
%
% Sample definition shared by all statistics scripts. EXCLUDE_NOISY_SUBJECTS lists the
% participants removed in the cleaned-sample run (RUN_MODE = 'sensitivity_excl_noisy'):
% the two participants with the highest composite noise scores (see X00 and the Methods).
% The full-sample run (RUN_MODE = 'full_sample') never excludes anyone.
% ========================================================================

EXCLUDE_NOISY_SUBJECTS = [7, 23];

EXCLUSION_RATIONALE_intervention = struct( ...
    's7',  'Intervention X00 CompositeZ=1.93 (borderline) + task X00 z=1.75; Int1 lost 106/141 windows to boundary events.', ...
    's23', 'Task X00 CompositeZ=1.85 (37% window rejection); excluded here too so both pipelines drop the same subjects.');

EXPECTED_N_SUBJECTS_TOTAL = 30;   % full-sample design: 30 subjects, 2 conditions each

RUN_MODES_intervention = {'full_sample', 'sensitivity_excl_noisy'};
RUN_MODES_DISPLAY_intervention = struct('full_sample', 'FullSample', ...
                                        'sensitivity_excl_noisy', 'SensitivityExclNoisy');
