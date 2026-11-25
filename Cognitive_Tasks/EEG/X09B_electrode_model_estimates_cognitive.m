%% ========================================================================
% X09B_electrode_model_estimates_cognitive.m
%
% Purpose: companion to X09_electrode_level_statistics_cognitive.m. X09
% saves full pairwise estimates only for the electrode x Band x Metric
% combinations whose FDR family already survived; this script adds
% model-based estimated marginal means (EMMs) and all three pairwise
% contrast estimates/SEs/uncorrected p-values for every one of the 868
% electrode x band x metric models, matching the intervention pipeline's
% own equivalent output and supporting model-based EMM topomaps and
% full-scalp t-maps. It changes nothing about X09's or X10's own saved
% outputs, models, FDR, or Holm logic; it only adds a new, read-only-
% derived file, computed by refitting the identical 868 LMMs X09 already
% fit. The p_..._uncorrected columns are descriptive only and are not a
% substitute for X09's own FDR-gated, Holm-corrected ElectrodePairwiseT,
% which remains the only confirmatory/inferential pairwise output in this
% pipeline. (X10, this pipeline's descriptive topography script, plots
% raw descriptive condition means at all 62 electrodes by design; that is
% unrelated to, and unchanged by, this script.)
%
% Method: for each RUN_MODE x Task (Full/Nback/SART; 6 independent
% passes, never pooled): (1) loads that Task+RUN_MODE's already-saved
% electrode_statistics.mat (ElectrodeBandPowerT, ElectrodeOmnibusT,
% produced by X09); errors out if X09 has not already been run for that
% combination; (2) re-fits the same 868 models (<Metric> ~ Condition +
% (1|Subject), REML, REF_LEVEL='E'); (3) cross-checks every successful
% refit's omnibus F/DF1/DF2/p against X09's own saved ElectrodeOmnibusT
% row for that exact Electrode/Metric/Band, within floating-point
% tolerance, and stops with an error on any mismatch beyond tolerance;
% (4) extracts EMM_E/EMM_N/EMM_R (model-based estimated marginal means),
% Est_/SE_/p_uncorrected for all three pairwise contrasts, for all
% electrodes regardless of FDR/Holm outcome; (5) appends
% ElectrodeModelEstimatesT into the same electrode_statistics.mat X09
% already wrote (ElectrodeOmnibusT, ElectrodePairwiseT,
% ElectrodeBandPowerT and PARAMS5B are never rewritten), and writes
% electrode_model_estimates.csv to the same folder, using the same column
% names as the intervention pipeline's equivalent file.
%
% Requirements: MATLAB with Statistics and Machine Learning Toolbox
% (fitlme, anova(LinearMixedModel), coefTest).
% X09_electrode_level_statistics_cognitive.m must already have been run
% for every (Task, RUN_MODE) combination this script attempts.
%
% Inputs (read-only, not modified):
%   <ANALYSIS_ROOT>/04_subject_level/phase3_subject_level_bandpower.mat
%     (PARAMS3 only -- band-tier/band-name list, shared, Task-invariant)
%   <ANALYSIS_ROOT>/07_roi_electrode_definitions/
%     phase5_channel_roi_definitions.mat (ChannelDefT -- shared electrode order)
%   <ANALYSIS_ROOT>/<Full|Nback|SART>/09_electrode_stats/<RUN_MODE>/
%     electrode_statistics.mat (ElectrodeBandPowerT, ElectrodeOmnibusT --
%     produced by X09, read-only)
%
% Outputs (added, per Task per RUN_MODE; nothing pre-existing is overwritten):
%   <ANALYSIS_ROOT>/<Full|Nback|SART>/09_electrode_stats/<RUN_MODE>/
%       electrode_model_estimates.csv
%       electrode_statistics.mat        (ElectrodeModelEstimatesT appended)
%   <ANALYSIS_ROOT>/09_reports/
%       phase5B_electrode_model_estimates_report_<Full|Nback|SART>_<RunModeLabel>.txt
%
% RUN_MODE: full_sample / sensitivity_excl_noisy; run once per mode.
% ========================================================================

clear; clc;

%% ------------------------------------------------------------------
%  0. SETTINGS, PATHS, EXISTENCE CHECKS
% ------------------------------------------------------------------

config;
DATA_ROOT_FOR_SIBLING = DATA_DIR;
[DATA_ROOT_PARENT, ~] = fileparts(DATA_ROOT_FOR_SIBLING);
ANALYSIS_ROOT = OUT_DIR;

SUBJ_DIR   = fullfile(ANALYSIS_ROOT, '04_subject_level');
DEF_DIR    = fullfile(ANALYSIS_ROOT, '07_roi_electrode_definitions');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');

bandPowerFile = fullfile(SUBJ_DIR, 'phase3_subject_level_bandpower.mat');
roiDefFile    = fullfile(DEF_DIR, 'phase5_channel_roi_definitions.mat');

if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

run(fullfile(fileparts(mfilename('fullpath')), 'EXCLUSION_CONFIG_cognitive.m'));

fprintf('========================================================\n');
fprintf('X09B -- ELECTRODE-LEVEL MODEL-BASED ESTIMATES (EMM + contrasts), COGNITIVE-TASK DATA\n');
fprintf('Companion to X09. Re-fits the SAME 868 models X09 already fit, per Task per RUN_MODE,\n');
fprintf('purely to extract and save EMM_E/N/R + Est_/SE_/p_uncorrected for ALL electrodes.\n');
fprintf('========================================================\n');
fprintf('ANALYSIS_ROOT   : %s\n', ANALYSIS_ROOT);
fprintf('--------------------------------------------------------\n\n');

if ~exist(bandPowerFile, 'file')
    error('Phase5B_X09B:MissingInput', 'X04 band-power file not found: %s', bandPowerFile);
end
if ~exist(roiDefFile, 'file')
    error('Phase5B_X09B:MissingInput', 'X07 electrode/ROI definition file not found: %s', roiDefFile);
end
if isempty(ver('stats'))
    error('Phase5B_X09B:MissingToolbox', ...
        'Statistics and Machine Learning Toolbox not detected. This script requires fitlme, anova(LinearMixedModel), and coefTest.');
end

%% ------------------------------------------------------------------
%  1. SHARED (Task- and RUN_MODE-invariant) REFERENCE DATA
% ------------------------------------------------------------------

S3 = load(bandPowerFile, 'PARAMS3');
PARAMS3 = S3.PARAMS3;
clear S3;
ALL_BANDS = fieldnames(PARAMS3.BandDefinitionsHz);
nBands = numel(ALL_BANDS);

S5D = load(roiDefFile, 'ChannelDefT', 'PARAMS5D');
ChannelDefT = S5D.ChannelDefT;
clear S5D;
ELECTRODES = ChannelDefT.Label;
nElec = numel(ELECTRODES);
if nElec ~= 62
    error('Phase5B_X09B:ElecCount', 'Expected 62 electrodes from X07, found %d.', nElec);
end

METRICS = {'LogAbsolutePower', 'RelativePower'};
nMetrics = numel(METRICS);
CONDITION_ORDER = {'E', 'N', 'R'};
REF_LEVEL = CONDITION_ORDER{1};
% This script's own convention is reproduced VERBATIM from X09's own
% contrast labels (E-N, E-R, N-R) -- the SAME convention Intervention's
% own electrode_model_estimates.csv uses for its Est_Eminus*/SE_Eminus*/
% p_*_uncorrected column names, so no renaming/relabeling is needed on
% either side of this port.
CONTRASTS = {'E','N','E - N'; 'E','R','E - R'; 'N','R','N - R'};

TASK_TYPES = {'full', 'nback', 'sart'};
TASK_DISPLAY = struct('full', 'Full', 'nback', 'Nback', 'sart', 'SART');

FIT_TOL_F = 1e-4;     % relative tolerance for cross-checking refit FStat vs X09's saved FStat
FIT_TOL_P = 1e-6;     % absolute tolerance for cross-checking refit pValue vs X09's saved pValue

grandCrossCheckFails = 0;

%% ========================================================================
%  MAIN LOOP -- RUN_MODE x TASK, each a fully independent pass
% ========================================================================
for rmI = 1:numel(RUN_MODES_cognitive)
    RUN_MODE = RUN_MODES_cognitive{rmI};
    RUN_MODE_LABEL = RUN_MODES_DISPLAY_cognitive.(RUN_MODE);

    for taskI = 1:numel(TASK_TYPES)
        TASK_KEY = TASK_TYPES{taskI};
        TASK_LABEL = TASK_DISPLAY.(TASK_KEY);

        ELEC_DIR = fullfile(ANALYSIS_ROOT, TASK_LABEL, '09_electrode_stats', RUN_MODE);
        matFile = fullfile(ELEC_DIR, 'electrode_statistics.mat');

        fprintf('\n########################################################\n');
        fprintf('### TASK = %s | RUN_MODE = %s\n', TASK_LABEL, RUN_MODE_LABEL);
        fprintf('########################################################\n\n');

        if ~exist(matFile, 'file')
            error('Phase5B_X09B:MissingX09Output', ...
                ['X09''s saved output not found for Task=%s, RUN_MODE=%s:\n  %s\n' ...
                 'Run X09_electrode_level_statistics_cognitive.m first with RUN_MODE=''%s'' ' ...
                 'before running this companion script for that Task/RUN_MODE.'], ...
                TASK_LABEL, RUN_MODE_LABEL, matFile, RUN_MODE);
        end

        Sx = load(matFile, 'ElectrodeBandPowerT', 'ElectrodeOmnibusT', 'PARAMS5B');
        if ~all(isfield(Sx, {'ElectrodeBandPowerT','ElectrodeOmnibusT','PARAMS5B'}))
            error('Phase5B_X09B:MissingVariables', ...
                '%s does not contain ElectrodeBandPowerT/ElectrodeOmnibusT/PARAMS5B as expected from X09.', matFile);
        end
        ElectrodeBandPowerT = Sx.ElectrodeBandPowerT;
        ElectrodeOmnibusT_X09 = Sx.ElectrodeOmnibusT;
        clear Sx;

        if isempty(ElectrodeBandPowerT)
            warning('Task=%s, RUN_MODE=%s: ElectrodeBandPowerT from X09 is empty -- skipping.', TASK_LABEL, RUN_MODE_LABEL);
            continue;
        end
        ElectrodeBandPowerT.ChannelUC = upper(ElectrodeBandPowerT.Channel);

        fprintf('Loaded X09''s saved ElectrodeBandPowerT (%d rows) and ElectrodeOmnibusT (%d rows) for Task=%s, RUN_MODE=%s.\n', ...
            height(ElectrodeBandPowerT), height(ElectrodeOmnibusT_X09), TASK_LABEL, RUN_MODE_LABEL);

        %% ----------------------------------------------------------------
        %  REFIT all Electrode x Band x Metric models, cross-check against
        %  X09's own saved omnibus result, and extract EMM + contrasts.
        % ----------------------------------------------------------------
        estRows = cell(0, 15);
        crossCheckFails = {};
        nRefitOK = 0; nRefitFailed = 0; nCrossCheckOK = 0;

        for e = 1:nElec
            elecName = ELECTRODES{e};
            for m = 1:nMetrics
                metricName = METRICS{m};
                for b = 1:nBands
                    bandName = ALL_BANDS{b};
                    subT = ElectrodeBandPowerT(strcmp(ElectrodeBandPowerT.ChannelUC, upper(elecName)) & strcmp(ElectrodeBandPowerT.Band, bandName), :);
                    if isempty(subT)
                        nRefitFailed = nRefitFailed + 1;
                        continue;
                    end

                    subT.Subject = categorical(local_subject_numeric(subT.Subject));
                    subT.Condition = categorical(cellstr(subT.Condition), CONDITION_ORDER);
                    formula = sprintf('%s ~ Condition + (1|Subject)', metricName);

                    try
                        lme = fitlme(subT, formula, 'FitMethod', 'REML');
                    catch ME %#ok<NASGU>
                        nRefitFailed = nRefitFailed + 1;
                        continue;
                    end

                    % --- cross-check against X09's own saved omnibus row ---
                    x09Row = ElectrodeOmnibusT_X09(strcmp(ElectrodeOmnibusT_X09.Electrode, elecName) & ...
                        strcmp(ElectrodeOmnibusT_X09.Metric, metricName) & strcmp(ElectrodeOmnibusT_X09.Band, bandName), :);
                    if height(x09Row) ~= 1
                        crossCheckFails{end+1} = sprintf('%s/%s/%s: expected exactly 1 matching row in X09''s ElectrodeOmnibusT, found %d.', ...
                            elecName, metricName, bandName, height(x09Row)); %#ok<AGROW>
                        continue;
                    end
                    if ~strcmp(x09Row.Status{1}, 'OK')
                        % X09 itself failed to fit this model -- nothing to cross-check or export for it.
                        nRefitFailed = nRefitFailed + 1;
                        continue;
                    end
                    try
                        aov = anova(lme, 'DFMethod', 'satterthwaite');
                        rowIdx = find(strcmp(aov.Term, 'Condition'), 1);
                        fStat = aov.FStat(rowIdx); pVal = aov.pValue(rowIdx);
                    catch ME %#ok<NASGU>
                        nRefitFailed = nRefitFailed + 1;
                        continue;
                    end
                    fMismatch = abs(fStat - x09Row.FStat(1)) > FIT_TOL_F * max(1, abs(x09Row.FStat(1)));
                    pMismatch = abs(pVal - x09Row.pValue(1)) > FIT_TOL_P;
                    if fMismatch || pMismatch
                        crossCheckFails{end+1} = sprintf('%s/%s/%s: refit F=%.6g p=%.6g vs X09-saved F=%.6g p=%.6g (mismatch beyond tolerance).', ...
                            elecName, metricName, bandName, fStat, pVal, x09Row.FStat(1), x09Row.pValue(1)); %#ok<AGROW>
                        continue;
                    end
                    nCrossCheckOK = nCrossCheckOK + 1;

                    % --- EMM + contrast extraction (identical logic to
                    %     the intervention pipeline's own equivalent estRows loop) ---
                    cn = lme.CoefficientNames; beta = lme.Coefficients.Estimate; covB = lme.CoefficientCovariance;
                    L = zeros(3, numel(cn));
                    L(:, strcmp(cn, '(Intercept)')) = 1;
                    L(2, strcmp(cn, 'Condition_N')) = 1;
                    L(3, strcmp(cn, 'Condition_R')) = 1;
                    emm = L * beta;
                    ests = nan(1, 3); ses = nan(1, 3); ps = nan(1, 3);
                    for c = 1:3
                        C = local_condition_contrast(cn, 'Condition', REF_LEVEL, CONTRASTS{c, 2}, CONTRASTS{c, 1});
                        ests(c) = C * beta; ses(c) = sqrt(C * covB * C');
                        ps(c) = coefTest(lme, C, 0, 'DFMethod', 'satterthwaite');
                    end
                    estRows(end+1, :) = {elecName, metricName, bandName, emm(1), emm(2), emm(3), ...
                        ests(1), ests(2), ests(3), ses(1), ses(2), ses(3), ps(1), ps(2), ps(3)}; %#ok<AGROW>
                    nRefitOK = nRefitOK + 1;
                end
            end
            if mod(e, 10) == 0 || e == nElec
                fprintf('  ... %d/%d electrodes processed (Task=%s, RUN_MODE=%s; %d exported OK, %d cross-check fails so far)\n', ...
                    e, nElec, TASK_LABEL, RUN_MODE_LABEL, nRefitOK, numel(crossCheckFails));
            end
        end

        ElectrodeModelEstimatesT = local_build_table(estRows, [1 2 3], 4:15, ...
            {'Electrode','Metric','Band','EMM_E','EMM_N','EMM_R', ...
             'Est_EminusN','Est_EminusR','Est_NminusR','SE_EminusN','SE_EminusR','SE_NminusR', ...
             'p_EminusN_uncorrected','p_EminusR_uncorrected','p_NminusR_uncorrected'});

        fprintf('\nRefit summary (Task=%s, RUN_MODE=%s): %d models refit and exported OK, %d refit/skip failures, %d cross-check mismatches.\n', ...
            TASK_LABEL, RUN_MODE_LABEL, nRefitOK, nRefitFailed, numel(crossCheckFails));

        if ~isempty(crossCheckFails)
            grandCrossCheckFails = grandCrossCheckFails + numel(crossCheckFails);
            fprintf('*** CROSS-CHECK MISMATCHES (refit result does not match X09''s own saved result) ***\n');
            for i = 1:numel(crossCheckFails)
                fprintf('  - %s\n', crossCheckFails{i});
            end
            error('Phase5B_X09B:CrossCheckFailed', ...
                ['Task=%s, RUN_MODE=%s: %d refit(s) did not match X09''s own saved ElectrodeOmnibusT result ' ...
                 'within tolerance (F: rel %.1e, p: abs %.1e). Stopping before saving any output for this ' ...
                 'Task/RUN_MODE -- an EMM export that does not provably correspond to X09''s own fit would not ' ...
                 'be trustworthy. See mismatches above.'], TASK_LABEL, RUN_MODE_LABEL, numel(crossCheckFails), FIT_TOL_F, FIT_TOL_P);
        else
            fprintf('ALL %d refit models cross-checked EXACTLY against X09''s own saved omnibus result (within tolerance). No mismatches.\n', nCrossCheckOK);
        end

        %% ----------------------------------------------------------------
        %  SAVE -- append-only into X09's own .mat; new .csv
        % ----------------------------------------------------------------
        save(matFile, 'ElectrodeModelEstimatesT', '-append');
        writetable(ElectrodeModelEstimatesT, fullfile(ELEC_DIR, 'electrode_model_estimates.csv'));
        fprintf('\nSaved electrode_model_estimates.csv (%d rows) and appended ElectrodeModelEstimatesT into %s.\n', ...
            height(ElectrodeModelEstimatesT), matFile);

        %% ----------------------------------------------------------------
        %  REPORT
        % ----------------------------------------------------------------
        reportLines = {};
        reportLines{end+1} = repmat('=', 1, 60);
        reportLines{end+1} = sprintf('X09B -- ELECTRODE-LEVEL MODEL-BASED ESTIMATES REPORT (Task=%s, RUN_MODE=%s)', TASK_LABEL, RUN_MODE_LABEL);
        reportLines{end+1} = repmat('=', 1, 60);
        reportLines{end+1} = '';
        reportLines{end+1} = 'Companion to X09_electrode_level_statistics_cognitive.m -- does not modify X09';
        reportLines{end+1} = 'or X10, and does not change any FDR/Holm/inferential result. Re-fits the SAME';
        reportLines{end+1} = '868 models X09 already fit (from X09''s own saved, already-filtered';
        reportLines{end+1} = 'ElectrodeBandPowerT), purely to export EMM_E/N/R + Est_/SE_/p_uncorrected for';
        reportLines{end+1} = 'EVERY electrode (ungated, descriptive -- same status as Intervention''s own';
        reportLines{end+1} = 'electrode_model_estimates.csv). The uncorrected p-values here are NOT a';
        reportLines{end+1} = 'substitute for X09''s own FDR-gated, Holm-corrected ElectrodePairwiseT.';
        reportLines{end+1} = '';
        reportLines{end+1} = sprintf('Models refit and exported: %d / %d expected (%d electrodes x %d bands x %d metrics)', ...
            nRefitOK, nElec*nBands*nMetrics, nElec, nBands, nMetrics);
        reportLines{end+1} = sprintf('Refit/skip failures (model also not OK in X09, or newly failed to refit): %d', nRefitFailed);
        reportLines{end+1} = sprintf('Cross-check against X09''s own saved omnibus result: %d checked, %d mismatches.', nCrossCheckOK, numel(crossCheckFails));
        if isempty(crossCheckFails)
            reportLines{end+1} = 'Every refit model reproduced X09''s own F-statistic and p-value within tolerance.';
        end
        reportLines{end+1} = '';
        reportLines{end+1} = 'Output:';
        reportLines{end+1} = ['  ' fullfile(ELEC_DIR, 'electrode_model_estimates.csv')];
        reportLines{end+1} = ['  ' matFile ' (ElectrodeModelEstimatesT appended)'];
        reportLines{end+1} = repmat('=', 1, 60);

        reportText = strjoin(reportLines, '\n');
        fprintf('%s\n', reportText);

        fid1 = fopen(fullfile(REPORT_DIR, sprintf('phase5B_electrode_model_estimates_report_%s_%s.txt', TASK_LABEL, RUN_MODE_LABEL)), 'w');
        fprintf(fid1, '%s\n', reportText);
        fclose(fid1);
    end
end

fprintf('\n=== X09B complete for all RUN_MODEs x all three Tasks. ===\n');
if grandCrossCheckFails == 0
    fprintf('No Task/RUN_MODE combination had any cross-check mismatch.\n');
end

%% ======================= LOCAL FUNCTIONS (all after top-level code) =======================

function n = local_subject_numeric(subjVal)
% Identical implementation to X09_electrode_level_statistics_cognitive.m.
    if iscategorical(subjVal)
        subjVal = cellstr(subjVal);
    end
    if iscell(subjVal)
        n = cellfun(@(s) str2double(regexprep(char(s), '\D', '')), subjVal);
    elseif isnumeric(subjVal)
        n = double(subjVal);
    else
        n = str2double(regexprep(char(subjVal), '\D', ''));
    end
end

function cRow = local_condition_contrast(coeffNames, varName, refLevel, level1, level2)
% Identical implementation to X09_electrode_level_statistics_cognitive.m.
    cRow = zeros(1, numel(coeffNames));
    if ~strcmp(level1, refLevel)
        name1 = sprintf('%s_%s', varName, level1);
        idx1 = find(strcmp(coeffNames, name1));
        if isempty(idx1)
            error('Phase5B_X09B:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name1);
        end
        cRow(idx1) = cRow(idx1) - 1;
    end
    if ~strcmp(level2, refLevel)
        name2 = sprintf('%s_%s', varName, level2);
        idx2 = find(strcmp(coeffNames, name2));
        if isempty(idx2)
            error('Phase5B_X09B:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name2);
        end
        cRow(idx2) = cRow(idx2) + 1;
    end
end

function T = local_build_table(rows, textColIdx, numColIdx, varNames)
% Identical implementation to X09_electrode_level_statistics_cognitive.m.
    nCols = numel(varNames);
    if size(rows, 2) ~= nCols
        error('Phase5B_X09B:BuildTableColumnMismatch', ...
            'local_build_table: rows has %d columns but %d variable names were given.', size(rows, 2), nCols);
    end
    cols = cell(1, nCols);
    for j = 1:nCols
        if ismember(j, textColIdx)
            cols{j} = rows(:, j);
        elseif ismember(j, numColIdx)
            if isempty(rows)
                cols{j} = zeros(0, 1);
            else
                cols{j} = cell2mat(rows(:, j));
            end
        else
            error('Phase5B_X09B:BuildTableColumnMismatch', 'Column %d not classified as text or numeric.', j);
        end
    end
    T = table(cols{:}, 'VariableNames', varNames);
end
