%% ========================================================================
% X12B_electrode_ratio_estimates_cognitive.m
%
% Purpose: companion to X12_band_ratio_indices_cognitive.m, same pattern
% as X11B. X12's own pairwise-contrasts output only contains rows for
% (Electrode x RatioType x Variant) combinations whose own FDR family
% already survived BH-FDR -- not enough for a full-scalp topomap. This
% script re-derives X12's electrode-level inputs and refits, without
% re-running X11 or X12. It changes nothing about X12's own models, FDR,
% or Holm logic, and does not touch X05/X05B/X09/X09B/X11/X11B/X17/X18 or
% their companions.
%
% Method: for each RUN_MODE x Task (Full/Nback/SART; 6 independent
% passes, never pooled): (1) re-derives X12's own Task+RUN_MODE-filtered
% electrode-level BandPowerT from X04's original, unfiltered BandPowerT
% (filter to TaskType, exclude Baseline, then, for
% sensitivity_excl_noisy only, exclude EXCLUDE_NOISY_SUBJECTS); (2) loads
% that Task+RUN_MODE's already-saved median_robust_data.mat
% (MedianRobustElectrodeT, produced by X11; errors out if X11 has not
% already been run); (3) computes ElecRatio_mean (from re-derived
% BandPowerT.AbsolutePower) and ElecRatio_median (from X11's
% MedianRobustElectrodeT.AbsMedian), using X12's own band-ratio formulas
% (ThetaBeta=Theta/Beta, AlphaTheta=Alpha/Theta, AlphaBeta=Alpha/Beta,
% EngagementIndex=Beta/(Alpha+Theta)); (4) re-fits the same Electrode x
% RatioType x Variant models X12 already fit (RatioValue ~ Condition +
% (1|Subject), REML, REF_LEVEL='E'), for all 62 electrodes x 4 RatioTypes
% x 2 Variants (496 models); (5) cross-checks every successful refit's
% omnibus F/p against X12's own saved electrode-omnibus results for that
% exact Channel/RatioType/Variant, within floating-point tolerance, and
% stops with an error on any mismatch beyond tolerance; (6) computes all
% 3 pairwise contrasts (N-E, R-E, R-N, X12's own convention) for every
% electrode, ungated, with Holm-Bonferroni applied within each
% electrode/RatioType/Variant model's own 3 contrasts
% (pHolm_withinModel, for reference only); (7) appends
% ElectrodeAllContrastsT_X12 into the same
% phase12_band_ratio_electrode_omnibus.mat X12 already wrote, and writes
% electrode_all_contrasts_ungated.csv to the same folder.
%
% Requirements: MATLAB with Statistics and Machine Learning Toolbox.
% X11_median_robust_family_cognitive.m and
% X12_band_ratio_indices_cognitive.m must already have been run for every
% (Task, RUN_MODE) combination this script attempts.
%
% Inputs (read-only, not modified):
%   <ANALYSIS_ROOT>/04_subject_level/phase3_subject_level_bandpower.mat
%     (BandPowerT, PARAMS3 -- shared, Task-invariant source)
%   <ANALYSIS_ROOT>/<Full|Nback|SART>/11_median_robust_family/<RUN_MODE>/
%     median_robust_data.mat (MedianRobustElectrodeT -- produced by X11)
%   <ANALYSIS_ROOT>/<Full|Nback|SART>/12_band_ratio_indices/<RUN_MODE>/
%     phase12_band_ratio_electrode_omnibus.mat (ElecOmnibusT -- produced
%     by X12, read-only)
%
% Outputs (added, per Task per RUN_MODE; nothing pre-existing is overwritten):
%   <ANALYSIS_ROOT>/<Full|Nback|SART>/12_band_ratio_indices/<RUN_MODE>/
%       electrode_all_contrasts_ungated.csv
%       phase12_band_ratio_electrode_omnibus.mat        (ElectrodeAllContrastsT_X12 appended)
%   <ANALYSIS_ROOT>/09_reports/
%       phase12B_electrode_all_contrasts_report_<Full|Nback|SART>_<RunModeLabel>.txt
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
fprintf('X12B -- ELECTRODE-LEVEL RATIO ALL-CONTRASTS, UNGATED (every electrode), COGNITIVE-TASK DATA\n');
fprintf('Companion to X12. Re-derives X12''s own filtered electrode input, re-fits the SAME\n');
fprintf('Electrode x RatioType x Variant models X12 already fit, per Task per RUN_MODE, purely\n');
fprintf('to export all 3 pairwise contrasts for EVERY electrode (ungated, descriptive).\n');
fprintf('========================================================\n');
fprintf('ANALYSIS_ROOT   : %s\n', ANALYSIS_ROOT);
fprintf('--------------------------------------------------------\n\n');

if ~exist(bandPowerFile, 'file')
    error('Phase12B:MissingInput', 'X04 band-power file not found: %s', bandPowerFile);
end
if isempty(ver('stats'))
    error('Phase12B:MissingToolbox', ...
        'Statistics and Machine Learning Toolbox not detected. This script requires fitlme, anova(LinearMixedModel), and coefTest.');
end

%% ------------------------------------------------------------------
%  1. SHARED (Task- and RUN_MODE-invariant) REFERENCE DATA
% ------------------------------------------------------------------

S3 = load(bandPowerFile, 'BandPowerT', 'PARAMS3');
BandPowerT_ALL = S3.BandPowerT;
PARAMS3 = S3.PARAMS3; %#ok<NASGU>
clear S3;

requiredCols = {'Subject','TaskType','BlockNumber','Condition','Channel','Band','AbsolutePower'};
missingCols = requiredCols(~ismember(requiredCols, BandPowerT_ALL.Properties.VariableNames));
if ~isempty(missingCols)
    error('Phase12B:MissingColumns', 'BandPowerT is missing expected column(s): %s', strjoin(missingCols, ', '));
end

RATIO_TYPES = {'ThetaBeta', 'AlphaTheta', 'AlphaBeta', 'EngagementIndex'};
VARIANTS    = {'MeanBased', 'MedianBased'};
CONDITION_ORDER = {'E', 'N', 'R'};
REF_LEVEL = CONDITION_ORDER{1};
% Reproduced VERBATIM from X12's own CONTRASTS declaration (N-E, R-E,
% R-N) -- NOT relabeled. See disclosure above.
CONTRASTS = {'N','E','N - E'; 'R','E','R - E'; 'R','N','R - N'};

TASK_TYPES = {'full', 'nback', 'sart'};
TASK_DISPLAY = struct('full', 'Full', 'nback', 'Nback', 'sart', 'SART');

FIT_TOL_F = 1e-4;     % relative tolerance for cross-checking refit FStat vs X12's saved FStat
FIT_TOL_P = 1e-6;     % absolute tolerance for cross-checking refit pValue vs X12's saved pValue

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

        MEDIAN_DIR = fullfile(ANALYSIS_ROOT, TASK_LABEL, '11_median_robust_family', RUN_MODE);
        STATS_DIR  = fullfile(ANALYSIS_ROOT, TASK_LABEL, '12_band_ratio_indices', RUN_MODE);
        medianFile  = fullfile(MEDIAN_DIR, 'median_robust_data.mat');
        omnibusFile = fullfile(STATS_DIR, 'phase12_band_ratio_electrode_omnibus.mat');

        fprintf('\n########################################################\n');
        fprintf('### TASK = %s | RUN_MODE = %s\n', TASK_LABEL, RUN_MODE_LABEL);
        fprintf('########################################################\n\n');

        if ~exist(medianFile, 'file')
            error('Phase12B:MissingX11Output', ...
                ['X11''s saved output not found for Task=%s, RUN_MODE=%s:\n  %s\n' ...
                 'Run X11_median_robust_family_cognitive.m first with RUN_MODE=''%s'' ' ...
                 'before running this companion script for that Task/RUN_MODE.'], ...
                TASK_LABEL, RUN_MODE_LABEL, medianFile, RUN_MODE);
        end
        if ~exist(omnibusFile, 'file')
            error('Phase12B:MissingX12Output', ...
                ['X12''s saved output not found for Task=%s, RUN_MODE=%s:\n  %s\n' ...
                 'Run X12_band_ratio_indices_cognitive.m first with RUN_MODE=''%s'' ' ...
                 'before running this companion script for that Task/RUN_MODE.'], ...
                TASK_LABEL, RUN_MODE_LABEL, omnibusFile, RUN_MODE);
        end

        Sm = load(medianFile, 'MedianRobustElectrodeT');
        if ~isfield(Sm, 'MedianRobustElectrodeT')
            error('Phase12B:MissingVariables', '%s does not contain MedianRobustElectrodeT as expected from X11.', medianFile);
        end
        MedianRobustElectrodeT = Sm.MedianRobustElectrodeT;
        clear Sm;

        So = load(omnibusFile, 'ElecOmnibusT');
        if ~isfield(So, 'ElecOmnibusT')
            error('Phase12B:MissingVariables', '%s does not contain ElecOmnibusT as expected from X12.', omnibusFile);
        end
        ElecOmnibusT_X12 = So.ElecOmnibusT;
        clear So;

        %% --- Re-derive X12's own Task+RUN_MODE-filtered electrode BandPowerT,
        %      step for step, from X04's original, unfiltered table -----
        thisTaskElecT = BandPowerT_ALL(strcmp(BandPowerT_ALL.TaskType, TASK_KEY), :);
        nBaselineExcludedElec = sum(strcmp(thisTaskElecT.Condition, 'Baseline'));
        BandPowerT = thisTaskElecT(ismember(thisTaskElecT.Condition, CONDITION_ORDER), :);

        nExcludedRowsSensitivityElec = 0;
        if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
            subjNumForExclElec = local_subject_numeric(BandPowerT.Subject);
            excludeMaskElec = ismember(subjNumForExclElec, EXCLUDE_NOISY_SUBJECTS);
            nExcludedRowsSensitivityElec = sum(excludeMaskElec);
            BandPowerT = BandPowerT(~excludeMaskElec, :);
        end
        fprintf('Re-derived Task=%s, RUN_MODE=%s filtered electrode input: %d total rows, %d Baseline excluded, %d sensitivity-excluded, %d rows remain.\n', ...
            TASK_LABEL, RUN_MODE_LABEL, height(thisTaskElecT), nBaselineExcludedElec, nExcludedRowsSensitivityElec, height(BandPowerT));

        if isempty(BandPowerT) || isempty(MedianRobustElectrodeT)
            warning('Task=%s, RUN_MODE=%s: no E/N/R rows after filtering -- skipping.', TASK_LABEL, RUN_MODE_LABEL);
            continue;
        end

        %% --- Compute the four ratios, both variants, electrode level
        %      (VERBATIM logic from X12's own local_compute_band_ratios) ---
        [ElecRatio_mean, nSkipMean, nNonFiniteMean] = local_compute_band_ratios(BandPowerT, {'Subject','Condition','Channel'}, 'Band', 'AbsolutePower');
        [ElecRatio_median, nSkipMed, nNonFiniteMed] = local_compute_band_ratios(MedianRobustElectrodeT, {'Subject','Condition','Channel'}, 'Band', 'AbsMedian');
        ElecRatio_mean.Variant   = repmat({'MeanBased'},   height(ElecRatio_mean),   1);
        ElecRatio_median.Variant = repmat({'MedianBased'}, height(ElecRatio_median), 1);
        ElecRatioT = [ElecRatio_mean; ElecRatio_median];
        fprintf('Computed ElecRatioT (Task=%s, RUN_MODE=%s): %d rows (MeanBased skipped=%d/nonfinite=%d; MedianBased skipped=%d/nonfinite=%d).\n', ...
            TASK_LABEL, RUN_MODE_LABEL, height(ElecRatioT), nSkipMean, nNonFiniteMean, nSkipMed, nNonFiniteMed);

        ELECTRODES = unique(ElecRatioT.Channel);
        nElec = numel(ELECTRODES);

        %% ----------------------------------------------------------------
        %  REFIT all Electrode x RatioType x Variant models, cross-check
        %  against X12's own saved omnibus result, and extract all 3
        %  contrasts (ungated).
        % ----------------------------------------------------------------
        allRows = cell(0, 12);
        crossCheckFails = {};
        nRefitOK = 0; nRefitFailed = 0; nCrossCheckOK = 0;

        for e = 1:nElec
            elecName = ELECTRODES{e};
            for vi = 1:numel(VARIANTS)
                variantName = VARIANTS{vi};
                for ri = 1:numel(RATIO_TYPES)
                    ratioName = RATIO_TYPES{ri};
                    subT = ElecRatioT(strcmp(ElecRatioT.Channel, elecName) & strcmp(ElecRatioT.Variant, variantName) & strcmp(ElecRatioT.RatioType, ratioName), :);
                    if isempty(subT)
                        nRefitFailed = nRefitFailed + 1;
                        continue;
                    end

                    subT.Subject = categorical(local_subject_numeric(subT.Subject));
                    subT.Condition = categorical(cellstr(subT.Condition), CONDITION_ORDER);

                    try
                        lme = fitlme(subT, 'RatioValue ~ Condition + (1|Subject)', 'FitMethod', 'REML');
                    catch ME %#ok<NASGU>
                        nRefitFailed = nRefitFailed + 1;
                        continue;
                    end

                    % --- cross-check against X12's own saved electrode omnibus row ---
                    x12Row = ElecOmnibusT_X12(strcmp(ElecOmnibusT_X12.Channel, elecName) & ...
                        strcmp(ElecOmnibusT_X12.Variant, variantName) & strcmp(ElecOmnibusT_X12.RatioType, ratioName), :);
                    if height(x12Row) ~= 1
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
                    fMismatch = abs(fStat - x12Row.FStat(1)) > FIT_TOL_F * max(1, abs(x12Row.FStat(1)));
                    pMismatch = abs(pVal - x12Row.pValue(1)) > FIT_TOL_P;
                    if fMismatch || pMismatch
                        crossCheckFails{end+1} = sprintf('%s/%s/%s: refit F=%.6g p=%.6g vs X12-saved F=%.6g p=%.6g (mismatch beyond tolerance).', ...
                            elecName, variantName, ratioName, fStat, pVal, x12Row.FStat(1), x12Row.pValue(1)); %#ok<AGROW>
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
                        allRows(end+1, :) = [{elecName, variantName, ratioName}, cdat(c, 1:6), {adj(c)}, cdat(c, 7:8)]; %#ok<AGROW>
                    end
                    nRefitOK = nRefitOK + 1;
                end
            end
            if mod(e, 10) == 0 || e == nElec
                fprintf('  ... %d/%d electrodes processed (Task=%s, RUN_MODE=%s; %d models exported OK, %d cross-check fails so far)\n', ...
                    e, nElec, TASK_LABEL, RUN_MODE_LABEL, nRefitOK, numel(crossCheckFails));
            end
        end

        ElectrodeAllContrastsT_X12 = local_build_table(allRows, [1 2 3 4], [5 6 7 8 9 10 11 12], ...
            {'Channel','Variant','RatioType','Contrast','Estimate','SE','DF','tStat','pValue','pHolm_withinModel','CI_Low','CI_High'});

        fprintf('\nRefit summary (Task=%s, RUN_MODE=%s): %d models refit and exported OK, %d refit/skip failures, %d cross-check mismatches.\n', ...
            TASK_LABEL, RUN_MODE_LABEL, nRefitOK, nRefitFailed, numel(crossCheckFails));

        if ~isempty(crossCheckFails)
            grandCrossCheckFails = grandCrossCheckFails + numel(crossCheckFails);
            fprintf('*** CROSS-CHECK MISMATCHES (refit result does not match X12''s own saved result) ***\n');
            for i = 1:numel(crossCheckFails)
                fprintf('  - %s\n', crossCheckFails{i});
            end
            error('Phase12B:CrossCheckFailed', ...
                ['Task=%s, RUN_MODE=%s: %d refit(s) did not match X12''s own saved ElecOmnibusT result ' ...
                 'within tolerance (F: rel %.1e, p: abs %.1e). Stopping before saving any output for this ' ...
                 'Task/RUN_MODE.'], TASK_LABEL, RUN_MODE_LABEL, numel(crossCheckFails), FIT_TOL_F, FIT_TOL_P);
        else
            fprintf('ALL %d refit models cross-checked EXACTLY against X12''s own saved omnibus result (within tolerance). No mismatches.\n', nCrossCheckOK);
        end

        %% ----------------------------------------------------------------
        %  SAVE -- append-only into X12's own .mat; new .csv
        % ----------------------------------------------------------------
        save(omnibusFile, 'ElectrodeAllContrastsT_X12', '-append');
        writetable(ElectrodeAllContrastsT_X12, fullfile(STATS_DIR, 'electrode_all_contrasts_ungated.csv'));
        fprintf('\nSaved electrode_all_contrasts_ungated.csv (%d rows) and appended ElectrodeAllContrastsT_X12 into %s.\n', ...
            height(ElectrodeAllContrastsT_X12), omnibusFile);

        %% ----------------------------------------------------------------
        %  REPORT
        % ----------------------------------------------------------------
        reportLines = {};
        reportLines{end+1} = repmat('=', 1, 60);
        reportLines{end+1} = sprintf('X12B -- ELECTRODE-LEVEL RATIO ALL-CONTRASTS UNGATED REPORT (Task=%s, RUN_MODE=%s)', TASK_LABEL, RUN_MODE_LABEL);
        reportLines{end+1} = repmat('=', 1, 60);
        reportLines{end+1} = '';
        reportLines{end+1} = 'Companion to X12_band_ratio_indices_cognitive.m -- does not modify X11 or X12 and';
        reportLines{end+1} = 'does not change any FDR/Holm/inferential result. Re-derives X12''s own filtered';
        reportLines{end+1} = 'electrode input and re-fits the SAME electrode-level models X12 already fit, purely';
        reportLines{end+1} = 'to export all 3 pairwise contrasts for EVERY electrode (ungated, descriptive --';
        reportLines{end+1} = 'needed for a full-scalp topomap). pHolm_withinModel is NOT a substitute for X12''s';
        reportLines{end+1} = 'own (RatioType x Variant)-family BH-FDR-gated phase12_band_ratio_electrode_pairwise.';
        reportLines{end+1} = '';
        reportLines{end+1} = sprintf('Models refit and exported: %d / %d expected (%d electrodes x %d ratio types x %d variants)', ...
            nRefitOK, nElec*numel(RATIO_TYPES)*numel(VARIANTS), nElec, numel(RATIO_TYPES), numel(VARIANTS));
        reportLines{end+1} = sprintf('Refit/skip failures: %d', nRefitFailed);
        reportLines{end+1} = sprintf('Cross-check against X12''s own saved omnibus result: %d checked, %d mismatches.', nCrossCheckOK, numel(crossCheckFails));
        if isempty(crossCheckFails)
            reportLines{end+1} = 'Every refit model reproduced X12''s own F-statistic and p-value within tolerance.';
        end
        reportLines{end+1} = '';
        reportLines{end+1} = 'Output:';
        reportLines{end+1} = ['  ' fullfile(STATS_DIR, 'electrode_all_contrasts_ungated.csv')];
        reportLines{end+1} = ['  ' omnibusFile ' (ElectrodeAllContrastsT_X12 appended)'];
        reportLines{end+1} = repmat('=', 1, 60);

        reportText = strjoin(reportLines, '\n');
        fprintf('%s\n', reportText);

        fid1 = fopen(fullfile(REPORT_DIR, sprintf('phase12B_electrode_all_contrasts_report_%s_%s.txt', TASK_LABEL, RUN_MODE_LABEL)), 'w');
        fprintf(fid1, '%s\n', reportText);
        fclose(fid1);
    end
end

fprintf('\n=== X12B complete for all RUN_MODEs x all three Tasks. ===\n');
if grandCrossCheckFails == 0
    fprintf('No Task/RUN_MODE combination had any cross-check mismatch.\n');
end

%% ======================= LOCAL FUNCTIONS (all after top-level code) =======================

function n = local_subject_numeric(subjVal)
% Identical implementation to X12_band_ratio_indices_cognitive.m.
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

function key = local_make_key(T, keyVars)
% Identical implementation to X12's own local_make_key: builds a unique
% string key per row from the keyVars columns, for containers.Map lookups.
    n = height(T);
    key = cell(n, 1);
    parts = cell(n, numel(keyVars));
    for j = 1:numel(keyVars)
        col = T.(keyVars{j});
        if iscell(col)
            parts(:, j) = col;
        elseif isnumeric(col)
            parts(:, j) = cellstr(num2str(col));
        else
            parts(:, j) = cellstr(col);
        end
    end
    for i = 1:n
        key{i} = strjoin(parts(i, :), '||');
    end
end

function [RatioLongT, nSkippedMissingBand, nNonFiniteRatios] = local_compute_band_ratios(sourceT, keyVars, bandCol, valueCol)
% Identical implementation to X12_band_ratio_indices_cognitive.m's own
% local_compute_band_ratios. Computes the four band-ratio indices
% (ThetaBeta, AlphaTheta, AlphaBeta, EngagementIndex) from a long table
% sourceT. One Theta/Alpha/Beta triplet is required per unique combination
% of keyVars; a combination missing any of the three is skipped and
% counted (never silently dropped without a trace).
    keyAll = local_make_key(sourceT, keyVars);
    bandVals = cellstr(sourceT.(bandCol));
    valCol = sourceT.(valueCol);

    thetaMask = strcmp(bandVals, 'Theta');
    alphaMask = strcmp(bandVals, 'Alpha');
    betaMask  = strcmp(bandVals, 'Beta');

    thetaMap = containers.Map(keyAll(thetaMask), num2cell(valCol(thetaMask)));
    alphaMap = containers.Map(keyAll(alphaMask), num2cell(valCol(alphaMask)));
    betaMap  = containers.Map(keyAll(betaMask),  num2cell(valCol(betaMask)));

    baseT = unique(sourceT(:, keyVars), 'rows');
    baseKeys = local_make_key(baseT, keyVars);
    nBase = height(baseT);

    ratioNames = {'ThetaBeta', 'AlphaTheta', 'AlphaBeta', 'EngagementIndex'};
    nRatios = numel(ratioNames);

    rows = cell(nBase * nRatios, numel(keyVars) + 2);
    rowCount = 0;
    nSkippedMissingBand = 0;
    nNonFiniteRatios = 0;

    for i = 1:nBase
        k = baseKeys{i};
        if ~(isKey(thetaMap, k) && isKey(alphaMap, k) && isKey(betaMap, k))
            nSkippedMissingBand = nSkippedMissingBand + 1;
            continue;
        end
        th = thetaMap(k); al = alphaMap(k); be = betaMap(k);
        ratioVals = [th / be, al / th, al / be, be / (al + th)];
        baseRowCells = table2cell(baseT(i, :));
        for r = 1:nRatios
            if ~isfinite(ratioVals(r))
                nNonFiniteRatios = nNonFiniteRatios + 1;
            end
            rowCount = rowCount + 1;
            rows(rowCount, :) = [baseRowCells, {ratioNames{r}}, {ratioVals(r)}];
        end
    end
    rows = rows(1:rowCount, :);

    varNames = [keyVars, {'RatioType', 'RatioValue'}];
    textColIdx = 1:(numel(keyVars) + 1);
    numColIdx = numel(keyVars) + 2;
    RatioLongT = local_build_table(rows, textColIdx, numColIdx, varNames);
end

function [adjP, reject] = local_holm_bonferroni(pvals, alphaFDR)
% Identical implementation to X05/X11/X12's own local_holm_bonferroni.
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
% Identical implementation to X05/X11/X12's own local_condition_contrast.
    cRow = zeros(1, numel(coeffNames));
    if ~strcmp(level1, refLevel)
        name1 = sprintf('%s_%s', varName, level1);
        idx1 = find(strcmp(coeffNames, name1));
        if isempty(idx1)
            error('Phase12B:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name1);
        end
        cRow(idx1) = cRow(idx1) - 1;
    end
    if ~strcmp(level2, refLevel)
        name2 = sprintf('%s_%s', varName, level2);
        idx2 = find(strcmp(coeffNames, name2));
        if isempty(idx2)
            error('Phase12B:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name2);
        end
        cRow(idx2) = cRow(idx2) + 1;
    end
end

function T = local_build_table(rows, textColIdx, numColIdx, varNames)
% Identical implementation to X05/X11/X12's own local_build_table.
    nCols = numel(varNames);
    if size(rows, 2) ~= nCols
        error('Phase12B:BuildTableColumnMismatch', ...
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
            error('Phase12B:BuildTableColumnMismatch', 'Column %d not classified as text or numeric.', j);
        end
    end
    T = table(cols{:}, 'VariableNames', varNames);
end
