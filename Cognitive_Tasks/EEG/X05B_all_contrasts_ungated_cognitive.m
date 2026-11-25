%% ========================================================================
% X05B_all_contrasts_ungated_cognitive.m
%
% Purpose: companion to X05_primary_cognitive_statistics.m. X05 only
% reports pairwise contrasts for the Band x Metric combinations whose
% tiered-family FDR already survived; this script reports all 3 pairwise
% contrasts for every one of the 14 Band x Metric models, regardless of
% whether that model's omnibus test survived FDR, matching the
% intervention pipeline's own equivalent output and supporting a
% forest-plot figure that needs a point + CI for every band. It does not
% change X05's own models, FDR, or Holm logic, and does not touch
% X06/X08/X09/X12/X24/X25; it only adds a new, read-only-derived file per
% Task per RUN_MODE.
%
% Method: for each RUN_MODE x Task (Full/Nback/SART; 6 independent
% passes, never pooled): (1) re-derives the same Task+RUN_MODE-filtered
% GlobalBandPowerT X05 itself fits on, independently re-applied step for
% step from X04's original, unfiltered GlobalBandPowerT (filter to this
% TaskType, exclude Baseline, then, for sensitivity_excl_noisy only,
% exclude EXCLUDE_NOISY_SUBJECTS); (2) re-fits the same 14 models
% (<Metric> ~ Condition + (1|Subject), REML, REF_LEVEL='E'); (3)
% cross-checks every successful refit's omnibus F/p-value against X05's
% own already-saved OmnibusT for that exact Band/Metric, within
% floating-point tolerance, and stops with an error on any mismatch
% beyond tolerance; (4) computes all 3 pairwise contrasts for every
% band/metric model (ungated), using this pipeline's own contrast
% convention (CONTRASTS = {N-E, R-E, R-N}, REF_LEVEL='E'), with
% Holm-Bonferroni applied within each band/metric model's own 3 contrasts
% (pHolm_withinModel, for reference only -- these rows carry no
% inferential/confirmatory claim); (5) appends AllContrastsT_cognitive
% into the same phase4_omnibus_results.mat X05 already wrote (OmnibusT
% and SensAllBandsT are never rewritten), and writes
% phase4_all_contrasts_ungated.csv/.xlsx to the same folder.
%
% Requirements: MATLAB with Statistics and Machine Learning Toolbox
% (fitlme, anova(LinearMixedModel), coefTest).
% X05_primary_cognitive_statistics.m must already have been run for every
% (Task, RUN_MODE) combination this script attempts.
%
% Inputs (read-only, not modified):
%   <ANALYSIS_ROOT>/04_subject_level/phase3_subject_level_bandpower.mat
%     (GlobalBandPowerT, PARAMS3)
%   <ANALYSIS_ROOT>/<Full|Nback|SART>/05_primary_stats/<RUN_MODE>/
%     phase4_omnibus_results.mat (OmnibusT -- produced by X05, read-only)
%
% Outputs (added, per Task per RUN_MODE; nothing pre-existing is overwritten):
%   <ANALYSIS_ROOT>/<Full|Nback|SART>/05_primary_stats/<RUN_MODE>/
%       phase4_all_contrasts_ungated.csv
%       phase4_all_contrasts_ungated.xlsx
%       phase4_omnibus_results.mat          (AllContrastsT_cognitive appended)
%   <ANALYSIS_ROOT>/09_reports/
%       phase4_all_contrasts_ungated_report_<Full|Nback|SART>_<RunModeLabel>.txt
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
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');

bandPowerFile = fullfile(SUBJ_DIR, 'phase3_subject_level_bandpower.mat');

if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

run(fullfile(fileparts(mfilename('fullpath')), 'EXCLUSION_CONFIG_cognitive.m'));

fprintf('========================================================\n');
fprintf('X05B -- ALL-CONTRASTS, UNGATED (every Band x Metric model), COGNITIVE-TASK DATA\n');
fprintf('Companion to X05. Re-derives X05''s own filtered input and re-fits the SAME 14\n');
fprintf('models X05 already fit, per Task per RUN_MODE, purely to export all 3 pairwise\n');
fprintf('contrasts for EVERY band/metric model (ungated, descriptive).\n');
fprintf('========================================================\n');
fprintf('ANALYSIS_ROOT   : %s\n', ANALYSIS_ROOT);
fprintf('--------------------------------------------------------\n\n');

if ~exist(bandPowerFile, 'file')
    error('Phase4_X05B:MissingInput', 'X04 band-power file not found: %s', bandPowerFile);
end
if isempty(ver('stats'))
    error('Phase4_X05B:MissingToolbox', ...
        'Statistics and Machine Learning Toolbox not detected. This script requires fitlme, anova(LinearMixedModel), and coefTest.');
end

%% ------------------------------------------------------------------
%  1. SHARED (Task- and RUN_MODE-invariant) REFERENCE DATA
% ------------------------------------------------------------------

S3 = load(bandPowerFile, 'GlobalBandPowerT', 'PARAMS3');
GlobalBandPowerT_ALL = S3.GlobalBandPowerT;
PARAMS3 = S3.PARAMS3;
clear S3;

requiredCols = {'Subject','TaskType','BlockNumber','Condition','Band', ...
                'AbsolutePower','LogAbsolutePower','RelativePower'};
missingCols = requiredCols(~ismember(requiredCols, GlobalBandPowerT_ALL.Properties.VariableNames));
if ~isempty(missingCols)
    error('Phase4_X05B:MissingColumns', 'GlobalBandPowerT is missing expected column(s): %s', strjoin(missingCols, ', '));
end

ALL_BANDS = fieldnames(PARAMS3.BandDefinitionsHz);
nBands = numel(ALL_BANDS);
METRICS = {'LogAbsolutePower', 'RelativePower'};
nMetrics = numel(METRICS);
CONDITION_ORDER = {'E', 'N', 'R'};
REF_LEVEL = CONDITION_ORDER{1};
% This script's own convention is reproduced VERBATIM from X05's own
% CONTRASTS declaration (N-E, R-E, R-N) -- NOT relabeled to Intervention's
% E-N/E-R/N-R convention. See disclosure above.
CONTRASTS = {'N','E','N - E'; 'R','E','R - E'; 'R','N','R - N'};

TASK_TYPES = {'full', 'nback', 'sart'};
TASK_DISPLAY = struct('full', 'Full', 'nback', 'Nback', 'sart', 'SART');

FIT_TOL_F = 1e-4;     % relative tolerance for cross-checking refit FStat vs X05's saved FStat
FIT_TOL_P = 1e-6;     % absolute tolerance for cross-checking refit pValue vs X05's saved pValue

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

        STATS_DIR = fullfile(ANALYSIS_ROOT, TASK_LABEL, '05_primary_stats', RUN_MODE);
        omnibusFile = fullfile(STATS_DIR, 'phase4_omnibus_results.mat');

        fprintf('\n########################################################\n');
        fprintf('### TASK = %s | RUN_MODE = %s\n', TASK_LABEL, RUN_MODE_LABEL);
        fprintf('########################################################\n\n');

        if ~exist(omnibusFile, 'file')
            error('Phase4_X05B:MissingX05Output', ...
                ['X05''s saved output not found for Task=%s, RUN_MODE=%s:\n  %s\n' ...
                 'Run X05_primary_cognitive_statistics.m first with RUN_MODE=''%s'' ' ...
                 'before running this companion script for that Task/RUN_MODE.'], ...
                TASK_LABEL, RUN_MODE_LABEL, omnibusFile, RUN_MODE);
        end

        Sx = load(omnibusFile, 'OmnibusT');
        if ~isfield(Sx, 'OmnibusT')
            error('Phase4_X05B:MissingVariables', '%s does not contain OmnibusT as expected from X05.', omnibusFile);
        end
        OmnibusT_X05 = Sx.OmnibusT;
        clear Sx;

        %% --- Re-derive X05's own Task+RUN_MODE-filtered GlobalBandPowerT,
        %      step for step, from X04's original, unfiltered table -----
        thisTaskT = GlobalBandPowerT_ALL(strcmp(GlobalBandPowerT_ALL.TaskType, TASK_KEY), :);
        nBaselineExcluded = sum(strcmp(thisTaskT.Condition, 'Baseline'));
        GlobalBandPowerT = thisTaskT(ismember(thisTaskT.Condition, CONDITION_ORDER), :);

        nExcludedRowsSensitivity = 0;
        if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
            subjNumForExcl = local_subject_numeric(GlobalBandPowerT.Subject);
            excludeMaskSensitivity = ismember(subjNumForExcl, EXCLUDE_NOISY_SUBJECTS);
            nExcludedRowsSensitivity = sum(excludeMaskSensitivity);
            GlobalBandPowerT = GlobalBandPowerT(~excludeMaskSensitivity, :);
        end
        fprintf('Re-derived Task=%s, RUN_MODE=%s filtered input: %d total rows, %d Baseline excluded, %d sensitivity-excluded, %d rows remain.\n', ...
            TASK_LABEL, RUN_MODE_LABEL, height(thisTaskT), nBaselineExcluded, nExcludedRowsSensitivity, height(GlobalBandPowerT));

        if isempty(GlobalBandPowerT)
            warning('Task=%s, RUN_MODE=%s: no E/N/R rows after filtering -- skipping.', TASK_LABEL, RUN_MODE_LABEL);
            continue;
        end

        %% ----------------------------------------------------------------
        %  REFIT all Band x Metric models, cross-check against X05's own
        %  saved omnibus result, and extract all 3 contrasts (ungated).
        % ----------------------------------------------------------------
        allRows = cell(0, 11);
        crossCheckFails = {};
        nRefitOK = 0; nRefitFailed = 0; nCrossCheckOK = 0;

        for b = 1:nBands
            bandName = ALL_BANDS{b};
            for m = 1:nMetrics
                metricName = METRICS{m};
                subT = GlobalBandPowerT(strcmp(GlobalBandPowerT.Band, bandName), :);
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

                % --- cross-check against X05's own saved omnibus row ---
                x05Row = OmnibusT_X05(strcmp(OmnibusT_X05.Metric, metricName) & strcmp(OmnibusT_X05.Band, bandName), :);
                if height(x05Row) ~= 1
                    crossCheckFails{end+1} = sprintf('%s/%s: expected exactly 1 matching row in X05''s OmnibusT, found %d.', ...
                        metricName, bandName, height(x05Row)); %#ok<AGROW>
                    continue;
                end
                if isfield(x05Row, 'Status') && ~strcmp(x05Row.Status{1}, 'OK')
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
                fMismatch = abs(fStat - x05Row.FStat(1)) > FIT_TOL_F * max(1, abs(x05Row.FStat(1)));
                pMismatch = abs(pVal - x05Row.pValue(1)) > FIT_TOL_P;
                if fMismatch || pMismatch
                    crossCheckFails{end+1} = sprintf('%s/%s: refit F=%.6g p=%.6g vs X05-saved F=%.6g p=%.6g (mismatch beyond tolerance).', ...
                        metricName, bandName, fStat, pVal, x05Row.FStat(1), x05Row.pValue(1)); %#ok<AGROW>
                    continue;
                end
                nCrossCheckOK = nCrossCheckOK + 1;

                % --- all 3 contrasts, ungated, Holm within this model only
                %     (identical logic to the intervention pipeline's own
                %     equivalent allRows loop) ---
                cn = lme.CoefficientNames;
                cp = nan(3, 1); cdat = cell(3, 8);
                for c = 1:3
                    C = local_condition_contrast(cn, 'Condition', REF_LEVEL, CONTRASTS{c, 2}, CONTRASTS{c, 1});
                    est = C * lme.Coefficients.Estimate;
                    se = sqrt(C * lme.CoefficientCovariance * C');
                    [p, ~, ~, df2c] = coefTest(lme, C, 0, 'DFMethod', 'satterthwaite');
                    tc = tinv(0.975, df2c);
                    cp(c) = p;
                    cdat(c, :) = {CONTRASTS{c, 3}, est, se, df2c, est / se, p, est - tc * se, est + tc * se};
                end
                adj = local_holm_bonferroni(cp, 0.05);
                for c = 1:3
                    allRows(end+1, :) = [{metricName, bandName}, cdat(c, 1:6), {adj(c)}, cdat(c, 7:8)]; %#ok<AGROW>
                end
                nRefitOK = nRefitOK + 1;
            end
        end

        AllContrastsT_cognitive = local_build_table(allRows, [1 2 3], [4 5 6 7 8 9 10 11], ...
            {'Metric','Band','Contrast','Estimate','SE','DF','tStat','pValue','pHolm_withinModel','CI_Low','CI_High'});

        fprintf('\nRefit summary (Task=%s, RUN_MODE=%s): %d models refit and exported OK, %d refit/skip failures, %d cross-check mismatches.\n', ...
            TASK_LABEL, RUN_MODE_LABEL, nRefitOK, nRefitFailed, numel(crossCheckFails));

        if ~isempty(crossCheckFails)
            grandCrossCheckFails = grandCrossCheckFails + numel(crossCheckFails);
            fprintf('*** CROSS-CHECK MISMATCHES (refit result does not match X05''s own saved result) ***\n');
            for i = 1:numel(crossCheckFails)
                fprintf('  - %s\n', crossCheckFails{i});
            end
            error('Phase4_X05B:CrossCheckFailed', ...
                ['Task=%s, RUN_MODE=%s: %d refit(s) did not match X05''s own saved OmnibusT result within ' ...
                 'tolerance (F: rel %.1e, p: abs %.1e). Stopping before saving any output for this Task/RUN_MODE.'], ...
                TASK_LABEL, RUN_MODE_LABEL, numel(crossCheckFails), FIT_TOL_F, FIT_TOL_P);
        else
            fprintf('ALL %d refit models cross-checked EXACTLY against X05''s own saved omnibus result (within tolerance). No mismatches.\n', nCrossCheckOK);
        end

        %% ----------------------------------------------------------------
        %  SAVE -- append-only into X05's own .mat; new .csv/.xlsx
        % ----------------------------------------------------------------
        save(omnibusFile, 'AllContrastsT_cognitive', '-append');
        writetable(AllContrastsT_cognitive, fullfile(STATS_DIR, 'phase4_all_contrasts_ungated.csv'));
        writetable(AllContrastsT_cognitive, fullfile(STATS_DIR, 'phase4_all_contrasts_ungated.xlsx'));
        fprintf('\nSaved phase4_all_contrasts_ungated.csv/.xlsx (%d rows) and appended AllContrastsT_cognitive into %s.\n', ...
            height(AllContrastsT_cognitive), omnibusFile);

        %% ----------------------------------------------------------------
        %  REPORT
        % ----------------------------------------------------------------
        reportLines = {};
        reportLines{end+1} = repmat('=', 1, 60);
        reportLines{end+1} = sprintf('X05B -- ALL-CONTRASTS UNGATED REPORT (Task=%s, RUN_MODE=%s)', TASK_LABEL, RUN_MODE_LABEL);
        reportLines{end+1} = repmat('=', 1, 60);
        reportLines{end+1} = '';
        reportLines{end+1} = 'Companion to X05_primary_cognitive_statistics.m -- does not modify X05, X06,';
        reportLines{end+1} = 'X08, X09, X12, or X24/X25, and does not change any FDR/Holm/inferential result.';
        reportLines{end+1} = 'Re-fits the SAME 14 models X05 already fit (from an independently re-derived,';
        reportLines{end+1} = 'step-for-step-identical copy of X05''s own filtered input), purely to export all';
        reportLines{end+1} = '3 pairwise contrasts (this pipeline''s own N-E/R-E/R-N convention) for EVERY';
        reportLines{end+1} = 'band/metric model (ungated, descriptive -- same status as Intervention''s own';
        reportLines{end+1} = 'phase4_all_contrasts_ungated.xlsx). pHolm_withinModel is NOT a substitute for';
        reportLines{end+1} = 'X05''s own tiered-family BH-FDR-gated phase4_pairwise_contrasts.';
        reportLines{end+1} = '';
        reportLines{end+1} = sprintf('Models refit and exported: %d / %d expected (%d bands x %d metrics)', ...
            nRefitOK, nBands*nMetrics, nBands, nMetrics);
        reportLines{end+1} = sprintf('Refit/skip failures: %d', nRefitFailed);
        reportLines{end+1} = sprintf('Cross-check against X05''s own saved omnibus result: %d checked, %d mismatches.', nCrossCheckOK, numel(crossCheckFails));
        if isempty(crossCheckFails)
            reportLines{end+1} = 'Every refit model reproduced X05''s own F-statistic and p-value within tolerance.';
        end
        reportLines{end+1} = '';
        reportLines{end+1} = 'Output:';
        reportLines{end+1} = ['  ' fullfile(STATS_DIR, 'phase4_all_contrasts_ungated.csv')];
        reportLines{end+1} = ['  ' fullfile(STATS_DIR, 'phase4_all_contrasts_ungated.xlsx')];
        reportLines{end+1} = ['  ' omnibusFile ' (AllContrastsT_cognitive appended)'];
        reportLines{end+1} = repmat('=', 1, 60);

        reportText = strjoin(reportLines, '\n');
        fprintf('%s\n', reportText);

        fid1 = fopen(fullfile(REPORT_DIR, sprintf('phase4_all_contrasts_ungated_report_%s_%s.txt', TASK_LABEL, RUN_MODE_LABEL)), 'w');
        fprintf(fid1, '%s\n', reportText);
        fclose(fid1);
    end
end

fprintf('\n=== X05B complete for all RUN_MODEs x all three Tasks. ===\n');
if grandCrossCheckFails == 0
    fprintf('No Task/RUN_MODE combination had any cross-check mismatch.\n');
end

%% ======================= LOCAL FUNCTIONS (all after top-level code) =======================

function n = local_subject_numeric(subjVal)
% Identical implementation to X05_primary_cognitive_statistics.m.
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

function [adjP, reject] = local_holm_bonferroni(pvals, alphaFDR)
% Identical implementation to X05_primary_cognitive_statistics.m.
    pvals = pvals(:);
    m = numel(pvals);
    [sortedP, sortIdx] = sort(pvals);
    adjSorted = nan(m, 1);
    runningMax = 0;
    for i = 1:m
        val = sortedP(i) * (m - i + 1);
        runningMax = max(runningMax, val);
        adjSorted(i) = runningMax;
    end
    adjP = nan(m, 1);
    adjP(sortIdx) = adjSorted;
    adjP = min(adjP, 1);
    reject = adjP <= alphaFDR;
end

function cRow = local_condition_contrast(coeffNames, varName, refLevel, level1, level2)
% Identical implementation to X05_primary_cognitive_statistics.m.
    cRow = zeros(1, numel(coeffNames));
    if ~strcmp(level1, refLevel)
        name1 = sprintf('%s_%s', varName, level1);
        idx1 = find(strcmp(coeffNames, name1));
        if isempty(idx1)
            error('Phase4_X05B:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name1);
        end
        cRow(idx1) = cRow(idx1) - 1;
    end
    if ~strcmp(level2, refLevel)
        name2 = sprintf('%s_%s', varName, level2);
        idx2 = find(strcmp(coeffNames, name2));
        if isempty(idx2)
            error('Phase4_X05B:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name2);
        end
        cRow(idx2) = cRow(idx2) + 1;
    end
end

function T = local_build_table(rows, textColIdx, numColIdx, varNames)
% Identical implementation to X05_primary_cognitive_statistics.m.
    nCols = numel(varNames);
    if size(rows, 2) ~= nCols
        error('Phase4_X05B:BuildTableColumnMismatch', ...
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
            error('Phase4_X05B:BuildTableColumnMismatch', 'Column %d not classified as text or numeric.', j);
        end
    end
    T = table(cols{:}, 'VariableNames', varNames);
end
