%% ========================================================================
% X11B_electrode_model_estimates_cognitive.m
%
% Purpose: companion to X11_median_robust_family_cognitive.m, following
% the same pattern as X09B and X05B. X11's own pairwise-contrasts output
% only contains rows for (Electrode x Band x Metric) combinations whose
% own FDR family already survived BH-FDR -- enough for a forest plot, but
% not enough for a continuous full-scalp topomap, which needs an
% estimate/t value at every electrode. This script reads X11's already-
% saved per-electrode source table and refits, without re-running X11 or
% any earlier script. It changes nothing about X11's own models, FDR, or
% Holm logic, and does not touch X05/X05B/X09/X09B/X12/X17/X18 or their
% companions; it only adds a new, read-only-derived file per Task per
% RUN_MODE, electrode level only (ROI/Global are already small enough --
% 5 ROIs / 1 Global -- that the existing gated output is effectively
% "full scalp" at those levels).
%
% Method: for each RUN_MODE x Task (Full/Nback/SART; 6 independent
% passes, never pooled): (1) loads that Task+RUN_MODE's already-saved
% median_robust_data.mat (MedianRobustElectrodeT) and
% phase11_omnibus_band_results.mat (OmnibusBandT, filtered to
% Level=='Electrode'); errors out if X11 has not already been run for
% that combination; (2) re-fits the same Electrode x Band x Metric models
% X11 already fit (<Metric> ~ Condition + (1|Subject), REML,
% REF_LEVEL='E'), for all 62 electrodes x 9 bands x 4 band-dependent
% metrics (2232 models; broadband TotalPower metrics have no Band
% dimension and are not included here); (3) cross-checks every successful
% refit's omnibus F/DF1/DF2/p against X11's own saved OmnibusBandT row
% for that exact Unit/Metric/Band, within floating-point tolerance, and
% stops with an error on any mismatch beyond tolerance; (4) computes all
% 3 pairwise contrasts (N-E, R-E, R-N, X11's own convention) for every
% electrode, ungated, with Holm-Bonferroni applied within each
% electrode/band/metric model's own 3 contrasts (pHolm_withinModel, for
% reference only -- carries no inferential/confirmatory claim); (5)
% appends ElectrodeAllContrastsT_X11 into the same median_robust_data.mat
% X11 already wrote, and writes electrode_all_contrasts_ungated.csv to
% the same folder.
%
% Requirements: MATLAB with Statistics and Machine Learning Toolbox.
% X11_median_robust_family_cognitive.m must already have been run for
% every (Task, RUN_MODE) combination this script attempts.
%
% Inputs (read-only, not modified):
%   <ANALYSIS_ROOT>/04_subject_level/phase3_subject_level_bandpower.mat
%     (PARAMS3 only -- band list, shared, Task-invariant)
%   <ANALYSIS_ROOT>/<Full|Nback|SART>/11_median_robust_family/<RUN_MODE>/
%     median_robust_data.mat (MedianRobustElectrodeT -- produced by X11,
%     read-only)
%     phase11_omnibus_band_results.mat (OmnibusBandT -- produced by X11,
%     read-only)
%
% Outputs (added, per Task per RUN_MODE; nothing pre-existing is overwritten):
%   <ANALYSIS_ROOT>/<Full|Nback|SART>/11_median_robust_family/<RUN_MODE>/
%       electrode_all_contrasts_ungated.csv
%       median_robust_data.mat                (ElectrodeAllContrastsT_X11 appended)
%   <ANALYSIS_ROOT>/09_reports/
%       phase11B_electrode_all_contrasts_report_<Full|Nback|SART>_<RunModeLabel>.txt
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
fprintf('X11B -- ELECTRODE-LEVEL ALL-CONTRASTS, UNGATED (every electrode), COGNITIVE-TASK DATA\n');
fprintf('Companion to X11. Re-fits the SAME electrode-level models X11 already fit (from X11''s\n');
fprintf('own saved MedianRobustElectrodeT), per Task per RUN_MODE, purely to export all 3\n');
fprintf('pairwise contrasts for EVERY electrode (ungated, descriptive) -- needed for a full-scalp\n');
fprintf('topomap figure.\n');
fprintf('========================================================\n');
fprintf('ANALYSIS_ROOT   : %s\n', ANALYSIS_ROOT);
fprintf('--------------------------------------------------------\n\n');

if ~exist(bandPowerFile, 'file')
    error('Phase11B:MissingInput', 'X04 band-power file not found: %s', bandPowerFile);
end
if isempty(ver('stats'))
    error('Phase11B:MissingToolbox', ...
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

METRICS_BAND = {'AbsMedian', 'RelMedian', 'CV_Absolute', 'CV_Relative'};
nMetrics = numel(METRICS_BAND);
CONDITION_ORDER = {'E', 'N', 'R'};
REF_LEVEL = CONDITION_ORDER{1};
% Reproduced VERBATIM from X11's own CONTRASTS declaration (N-E, R-E,
% R-N) -- NOT relabeled. See disclosure above.
CONTRASTS = {'N','E','N - E'; 'R','E','R - E'; 'R','N','R - N'};

TASK_TYPES = {'full', 'nback', 'sart'};
TASK_DISPLAY = struct('full', 'Full', 'nback', 'Nback', 'sart', 'SART');

FIT_TOL_F = 1e-4;     % relative tolerance for cross-checking refit FStat vs X11's saved FStat
FIT_TOL_P = 1e-6;     % absolute tolerance for cross-checking refit pValue vs X11's saved pValue

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

        STATS_DIR = fullfile(ANALYSIS_ROOT, TASK_LABEL, '11_median_robust_family', RUN_MODE);
        dataFile = fullfile(STATS_DIR, 'median_robust_data.mat');
        omnibusFile = fullfile(STATS_DIR, 'phase11_omnibus_band_results.mat');

        fprintf('\n########################################################\n');
        fprintf('### TASK = %s | RUN_MODE = %s\n', TASK_LABEL, RUN_MODE_LABEL);
        fprintf('########################################################\n\n');

        if ~exist(dataFile, 'file') || ~exist(omnibusFile, 'file')
            error('Phase11B:MissingX11Output', ...
                ['X11''s saved output not found for Task=%s, RUN_MODE=%s:\n  %s\n  %s\n' ...
                 'Run X11_median_robust_family_cognitive.m first with RUN_MODE=''%s'' ' ...
                 'before running this companion script for that Task/RUN_MODE.'], ...
                TASK_LABEL, RUN_MODE_LABEL, dataFile, omnibusFile, RUN_MODE);
        end

        Sd = load(dataFile, 'MedianRobustElectrodeT');
        if ~isfield(Sd, 'MedianRobustElectrodeT')
            error('Phase11B:MissingVariables', '%s does not contain MedianRobustElectrodeT as expected from X11.', dataFile);
        end
        MedianRobustElectrodeT = Sd.MedianRobustElectrodeT;
        clear Sd;

        So = load(omnibusFile, 'OmnibusBandT');
        if ~isfield(So, 'OmnibusBandT')
            error('Phase11B:MissingVariables', '%s does not contain OmnibusBandT as expected from X11.', omnibusFile);
        end
        OmnibusBandT_X11 = So.OmnibusBandT(strcmp(So.OmnibusBandT.Level, 'Electrode'), :);
        clear So;

        if isempty(MedianRobustElectrodeT)
            warning('Task=%s, RUN_MODE=%s: MedianRobustElectrodeT from X11 is empty -- skipping.', TASK_LABEL, RUN_MODE_LABEL);
            continue;
        end
        ELECTRODES = unique(MedianRobustElectrodeT.Channel);
        nElec = numel(ELECTRODES);

        fprintf('Loaded X11''s saved MedianRobustElectrodeT (%d rows, %d electrodes) and electrode-level OmnibusBandT (%d rows) for Task=%s, RUN_MODE=%s.\n', ...
            height(MedianRobustElectrodeT), nElec, height(OmnibusBandT_X11), TASK_LABEL, RUN_MODE_LABEL);

        %% ----------------------------------------------------------------
        %  REFIT all Electrode x Band x Metric models, cross-check against
        %  X11's own saved omnibus result, and extract all 3 contrasts
        %  (ungated).
        % ----------------------------------------------------------------
        allRows = cell(0, 12);
        crossCheckFails = {};
        nRefitOK = 0; nRefitFailed = 0; nCrossCheckOK = 0;

        for e = 1:nElec
            elecName = ELECTRODES{e};
            for m = 1:nMetrics
                metricName = METRICS_BAND{m};
                for b = 1:nBands
                    bandName = ALL_BANDS{b};
                    subT = MedianRobustElectrodeT(strcmp(MedianRobustElectrodeT.Channel, elecName) & strcmp(MedianRobustElectrodeT.Band, bandName), :);
                    if isempty(subT) || ~ismember(metricName, subT.Properties.VariableNames)
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

                    % --- cross-check against X11's own saved electrode omnibus row ---
                    x11Row = OmnibusBandT_X11(strcmp(OmnibusBandT_X11.Unit, elecName) & ...
                        strcmp(OmnibusBandT_X11.Metric, metricName) & strcmp(OmnibusBandT_X11.Band, bandName), :);
                    if height(x11Row) ~= 1
                        % Not an error by itself -- X11 only writes a row here if that
                        % electrode's own fit succeeded AND was part of its FDR family
                        % pass (it always is, for every successfully fit electrode, since
                        % OmnibusBandT is ungated at the omnibus level already -- see
                        % X11's own Section 10). A genuinely missing row means X11's own
                        % fit failed for this electrode/band/metric; skip silently (not a
                        % cross-check failure), matching X11's own integrityIssues pattern.
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
                    fMismatch = abs(fStat - x11Row.FStat(1)) > FIT_TOL_F * max(1, abs(x11Row.FStat(1)));
                    pMismatch = abs(pVal - x11Row.pValue(1)) > FIT_TOL_P;
                    if fMismatch || pMismatch
                        crossCheckFails{end+1} = sprintf('%s/%s/%s: refit F=%.6g p=%.6g vs X11-saved F=%.6g p=%.6g (mismatch beyond tolerance).', ...
                            elecName, metricName, bandName, fStat, pVal, x11Row.FStat(1), x11Row.pValue(1)); %#ok<AGROW>
                        continue;
                    end
                    nCrossCheckOK = nCrossCheckOK + 1;

                    % --- all 3 contrasts, ungated, Holm within this model only ---
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
                        allRows(end+1, :) = [{elecName, metricName, bandName}, cdat(c, 1:6), {adj(c)}, cdat(c, 7:8)]; %#ok<AGROW>
                    end
                    nRefitOK = nRefitOK + 1;
                end
            end
            if mod(e, 10) == 0 || e == nElec
                fprintf('  ... %d/%d electrodes processed (Task=%s, RUN_MODE=%s; %d models exported OK, %d cross-check fails so far)\n', ...
                    e, nElec, TASK_LABEL, RUN_MODE_LABEL, nRefitOK, numel(crossCheckFails));
            end
        end

        ElectrodeAllContrastsT_X11 = local_build_table(allRows, [1 2 3 4], [5 6 7 8 9 10 11 12], ...
            {'Electrode','Metric','Band','Contrast','Estimate','SE','DF','tStat','pValue','pHolm_withinModel','CI_Low','CI_High'});

        fprintf('\nRefit summary (Task=%s, RUN_MODE=%s): %d models refit and exported OK, %d refit/skip failures, %d cross-check mismatches.\n', ...
            TASK_LABEL, RUN_MODE_LABEL, nRefitOK, nRefitFailed, numel(crossCheckFails));

        if ~isempty(crossCheckFails)
            grandCrossCheckFails = grandCrossCheckFails + numel(crossCheckFails);
            fprintf('*** CROSS-CHECK MISMATCHES (refit result does not match X11''s own saved result) ***\n');
            for i = 1:numel(crossCheckFails)
                fprintf('  - %s\n', crossCheckFails{i});
            end
            error('Phase11B:CrossCheckFailed', ...
                ['Task=%s, RUN_MODE=%s: %d refit(s) did not match X11''s own saved OmnibusBandT result ' ...
                 'within tolerance (F: rel %.1e, p: abs %.1e). Stopping before saving any output for this ' ...
                 'Task/RUN_MODE.'], TASK_LABEL, RUN_MODE_LABEL, numel(crossCheckFails), FIT_TOL_F, FIT_TOL_P);
        else
            fprintf('ALL %d refit models cross-checked EXACTLY against X11''s own saved omnibus result (within tolerance). No mismatches.\n', nCrossCheckOK);
        end

        %% ----------------------------------------------------------------
        %  SAVE -- append-only into X11's own .mat; new .csv
        % ----------------------------------------------------------------
        save(dataFile, 'ElectrodeAllContrastsT_X11', '-append');
        writetable(ElectrodeAllContrastsT_X11, fullfile(STATS_DIR, 'electrode_all_contrasts_ungated.csv'));
        fprintf('\nSaved electrode_all_contrasts_ungated.csv (%d rows) and appended ElectrodeAllContrastsT_X11 into %s.\n', ...
            height(ElectrodeAllContrastsT_X11), dataFile);

        %% ----------------------------------------------------------------
        %  REPORT
        % ----------------------------------------------------------------
        reportLines = {};
        reportLines{end+1} = repmat('=', 1, 60);
        reportLines{end+1} = sprintf('X11B -- ELECTRODE-LEVEL ALL-CONTRASTS UNGATED REPORT (Task=%s, RUN_MODE=%s)', TASK_LABEL, RUN_MODE_LABEL);
        reportLines{end+1} = repmat('=', 1, 60);
        reportLines{end+1} = '';
        reportLines{end+1} = 'Companion to X11_median_robust_family_cognitive.m -- does not modify X11 and';
        reportLines{end+1} = 'does not change any FDR/Holm/inferential result. Re-fits the SAME electrode-level';
        reportLines{end+1} = 'models X11 already fit (from X11''s own saved, already-filtered MedianRobustElectrodeT),';
        reportLines{end+1} = 'purely to export all 3 pairwise contrasts for EVERY electrode (ungated, descriptive --';
        reportLines{end+1} = 'needed for a full-scalp topomap). pHolm_withinModel is NOT a substitute for X11''s';
        reportLines{end+1} = 'own (Band x Metric)-family BH-FDR-gated phase11_pairwise_band_contrasts.';
        reportLines{end+1} = '';
        reportLines{end+1} = sprintf('Models refit and exported: %d / %d expected (%d electrodes x %d bands x %d metrics)', ...
            nRefitOK, nElec*nBands*nMetrics, nElec, nBands, nMetrics);
        reportLines{end+1} = sprintf('Refit/skip failures (model also not OK in X11, or newly failed to refit): %d', nRefitFailed);
        reportLines{end+1} = sprintf('Cross-check against X11''s own saved omnibus result: %d checked, %d mismatches.', nCrossCheckOK, numel(crossCheckFails));
        if isempty(crossCheckFails)
            reportLines{end+1} = 'Every refit model reproduced X11''s own F-statistic and p-value within tolerance.';
        end
        reportLines{end+1} = '';
        reportLines{end+1} = 'Output:';
        reportLines{end+1} = ['  ' fullfile(STATS_DIR, 'electrode_all_contrasts_ungated.csv')];
        reportLines{end+1} = ['  ' dataFile ' (ElectrodeAllContrastsT_X11 appended)'];
        reportLines{end+1} = repmat('=', 1, 60);

        reportText = strjoin(reportLines, '\n');
        fprintf('%s\n', reportText);

        fid1 = fopen(fullfile(REPORT_DIR, sprintf('phase11B_electrode_all_contrasts_report_%s_%s.txt', TASK_LABEL, RUN_MODE_LABEL)), 'w');
        fprintf(fid1, '%s\n', reportText);
        fclose(fid1);
    end
end

fprintf('\n=== X11B complete for all RUN_MODEs x all three Tasks. ===\n');
if grandCrossCheckFails == 0
    fprintf('No Task/RUN_MODE combination had any cross-check mismatch.\n');
end

%% ======================= LOCAL FUNCTIONS (all after top-level code) =======================

function n = local_subject_numeric(subjVal)
% Identical implementation to X11_median_robust_family_cognitive.m.
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
% Identical implementation to X05/X11's own local_holm_bonferroni.
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
% Identical implementation to X05/X11's own local_condition_contrast.
    cRow = zeros(1, numel(coeffNames));
    if ~strcmp(level1, refLevel)
        name1 = sprintf('%s_%s', varName, level1);
        idx1 = find(strcmp(coeffNames, name1));
        if isempty(idx1)
            error('Phase11B:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name1);
        end
        cRow(idx1) = cRow(idx1) - 1;
    end
    if ~strcmp(level2, refLevel)
        name2 = sprintf('%s_%s', varName, level2);
        idx2 = find(strcmp(coeffNames, name2));
        if isempty(idx2)
            error('Phase11B:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name2);
        end
        cRow(idx2) = cRow(idx2) + 1;
    end
end

function T = local_build_table(rows, textColIdx, numColIdx, varNames)
% Identical implementation to X05/X11's own local_build_table.
    nCols = numel(varNames);
    if size(rows, 2) ~= nCols
        error('Phase11B:BuildTableColumnMismatch', ...
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
            error('Phase11B:BuildTableColumnMismatch', 'Column %d not classified as text or numeric.', j);
        end
    end
    T = table(cols{:}, 'VariableNames', varNames);
end
