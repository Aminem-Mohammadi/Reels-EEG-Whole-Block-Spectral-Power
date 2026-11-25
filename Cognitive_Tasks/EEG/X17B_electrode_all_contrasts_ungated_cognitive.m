%% ========================================================================
% X17B_electrode_all_contrasts_ungated_cognitive.m
%
% COMPANION SCRIPT TO X17_half_split_condition_time_cognitive.m -- NOT a
% modification of X17, NOT a new numbered phase. A purely additive
% companion, same pattern and reason as X11B/X12B: closing a
% schema gap needed to build a full-scalp (ALL 62 electrodes) topographic
% figure for X17's own half-split Condition x TimeHalf family.
%
% WHY THIS SCRIPT EXISTS
%   X17's own phase17_pairwise_contrasts.csv only contains rows for
%   (Electrode x Band x Metric x Term) combinations whose OWN per-term FDR
%   family already survived BH-FDR. Not enough for a full-scalp topomap.
%   X17 DOES already save its own per-electrode, long-format source table
%   (half_split_data.mat -> ElectrodeT, with TimeHalf already derived/
%   attached, and already Task+RUN_MODE-filtered+sensitivity-excluded) to
%   disk -- verified by reading X17's own Section 9 save() call directly.
%   This script reads ONLY that already-saved table and refits; NO rerun
%   of X17 or any earlier script is needed, same situation as X11B's.
%   This script changes NOTHING about X17's own models, FDR, or Holm
%   logic, and does not touch X05/X05B/X09/X09B/X11/X11B/X12/X12B/X18/X18B
%   or their companions. Electrode level ONLY (Global level is a single
%   unit -- already effectively "full scalp" via the existing gated CSV).
%
% WHAT THIS SCRIPT DOES
%   For each RUN_MODE in RUN_MODES_cognitive x each Task in {Nback, SART}
%   (4 independent passes, never pooled -- Full is categorically out of
%   scope for this half-split analysis, identical restriction to X17
%   itself: TASK_TYPES = {'nback','sart'}):
%     1. Loads that Task+RUN_MODE's ALREADY-SAVED half_split_data.mat
%        (produced by X17) -- specifically ElectrodeT (Subject/TaskType/
%        BlockNumber/Condition/TimeHalf/Channel/Band/LogAbsolutePower/
%        LogAbsolutePower/RelativePower) -- and that Task+RUN_MODE's
%        ALREADY-SAVED phase17_omnibus_results.mat (OmnibusT), filtered to
%        Level=='Electrode', for cross-checking. REQUIRES X17 to already
%        have been run for that Task+RUN_MODE; errors out if missing.
%     2. Re-fits the SAME Electrode x Band x Metric two-factor models X17
%        already fit (formula "<Metric> ~ Condition*TimeHalf + (1|Subject)",
%        REML, Condition ref='E', TimeHalf ref='1', identical to X17's own
%        local_fit_two_factor_omnibus, reproduced verbatim below), for ALL
%        62 electrodes x 9 bands x 2 metrics (1116 models).
%     3. CROSS-CHECKS every successful refit's THREE omnibus F/DF1/DF2/p
%        (Condition, TimeHalf, Condition:TimeHalf) against X17's own
%        already-saved OmnibusT rows (Level='Electrode') for that exact
%        Unit/Metric/Band/Term, within floating-point tolerance. Any
%        mismatch beyond tolerance is a hard STOP -- same safeguard as
%        X09B/X05B/X11B/X12B.
%     4. From each successful refit, computes ALL THREE pairwise-contrast
%        families, ungated, for EVERY electrode, using X17's OWN
%        contrast-construction functions reproduced VERBATIM (NOT
%        re-derived differently):
%          - ConditionMarginal: N-E/R-E/R-N, averaged over TimeHalf,
%            Holm-3 (local_run_pairwise_condition_marginal)
%          - TimeHalfMarginal: TimeHalf=2 - TimeHalf=1, averaged over
%            Condition, single comparison, no correction
%            (local_run_pairwise_timehalf_marginal)
%          - InteractionSimple: N-E/R-E/R-N evaluated separately AT
%            TimeHalf=1 and TimeHalf=2 (6 contrasts), Holm-6
%            (local_run_pairwise_condition_simple)
%        -- all ungated (computed for every electrode regardless of that
%        term's own FDR outcome), tagged with a SurvivingTerm-style label
%        for direct comparability with X17's own column, but renamed
%        "ContrastFamily" since nothing here actually "survived" anything.
%     5. Saves ElectrodeAllContrastsT_X17:
%          - appended into the SAME half_split_data.mat X17 already wrote
%            (via save(...,'-append'));
%          - electrode_all_contrasts_ungated.csv, written to the SAME
%            17_half_split_condition_time/<RUN_MODE> folder X17 already
%            writes to.
%
% WHAT THIS SCRIPT DOES NOT DO
%   - Does not change X17's FDR, Holm, or any inferential conclusion.
%   - Does not touch Global-level output.
%   - Does not run X17 itself, and does not attempt Full (same categorical
%     exclusion as X17 itself -- logged, not silently skipped).
%
% REQUIREMENTS
%   - MATLAB with Statistics and Machine Learning Toolbox.
%   - X17_half_split_condition_time_cognitive.m must already have been run
%     for EVERY (Task, RUN_MODE) combination this script attempts.
%
% INPUT (read-only, not modified):
%   <ANALYSIS_ROOT>/<Nback|SART>/17_half_split_condition_time/<RUN_MODE>/
%     half_split_data.mat (ElectrodeT -- produced by X17, read-only)
%     phase17_omnibus_results.mat (OmnibusT -- produced by X17, read-only)
%
% OUTPUT (added, per Task per RUN_MODE; nothing pre-existing is overwritten):
%   <ANALYSIS_ROOT>/<Nback|SART>/17_half_split_condition_time/<RUN_MODE>/
%       electrode_all_contrasts_ungated.csv   (new)
%       half_split_data.mat                   (ElectrodeAllContrastsT_X17 appended)
%   <ANALYSIS_ROOT>/09_reports/
%       phase17B_electrode_all_contrasts_report_<Nback|SART>_<RunModeLabel>.txt
%
% The report, per Task per RUN_MODE, includes the cross-check section
% confirming the refit omnibus results match X17's own saved OmnibusT
% exactly.
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
fprintf('X17B -- ELECTRODE-LEVEL ALL-CONTRASTS, UNGATED (every electrode), COGNITIVE-TASK DATA\n');
fprintf('Companion to X17. Re-fits the SAME electrode-level two-factor models X17 already fit\n');
fprintf('(from X17''s own saved ElectrodeT), per Task per RUN_MODE, purely to export all three\n');
fprintf('contrast families (ConditionMarginal/TimeHalfMarginal/InteractionSimple) for EVERY\n');
fprintf('electrode (ungated, descriptive) -- needed for a full-scalp topomap figure.\n');
fprintf('========================================================\n');
fprintf('ANALYSIS_ROOT   : %s\n', ANALYSIS_ROOT);
fprintf('--------------------------------------------------------\n\n');

if ~exist(bandPowerFile, 'file')
    error('Phase17B:MissingInput', 'X04 band-power file not found: %s', bandPowerFile);
end
if isempty(ver('stats'))
    error('Phase17B:MissingToolbox', ...
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

METRICS = {'LogAbsolutePower', 'RelativePower'};
nMetrics = numel(METRICS);
CONDITION_ORDER = {'E', 'N', 'R'};
REF_LEVEL = CONDITION_ORDER{1};
TIMEHALF_ORDER = {'1', '2'};
REF_TIMEHALF = TIMEHALF_ORDER{1};
NONREF_TIMEHALF = TIMEHALF_ORDER{2};
% Reproduced VERBATIM from X17's own CONTRASTS declaration (N-E, R-E,
% R-N) -- NOT relabeled. See disclosure above.
CONTRASTS = {'N','E','N - E'; 'R','E','R - E'; 'R','N','R - N'};

TASK_TYPES = {'nback', 'sart'};   % Full is OUT OF SCOPE -- identical restriction to X17 itself
TASK_DISPLAY = struct('nback', 'Nback', 'sart', 'SART');

FIT_TOL_F = 1e-4;     % relative tolerance for cross-checking refit FStat vs X17's saved FStat
FIT_TOL_P = 1e-6;     % absolute tolerance for cross-checking refit pValue vs X17's saved pValue

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

        STATS_DIR = fullfile(ANALYSIS_ROOT, TASK_LABEL, '17_half_split_condition_time', RUN_MODE);
        dataFile = fullfile(STATS_DIR, 'half_split_data.mat');
        omnibusFile = fullfile(STATS_DIR, 'phase17_omnibus_results.mat');

        fprintf('\n########################################################\n');
        fprintf('### TASK = %s | RUN_MODE = %s\n', TASK_LABEL, RUN_MODE_LABEL);
        fprintf('########################################################\n\n');

        if ~exist(dataFile, 'file') || ~exist(omnibusFile, 'file')
            error('Phase17B:MissingX17Output', ...
                ['X17''s saved output not found for Task=%s, RUN_MODE=%s:\n  %s\n  %s\n' ...
                 'Run X17_half_split_condition_time_cognitive.m first with RUN_MODE=''%s'' ' ...
                 'before running this companion script for that Task/RUN_MODE.'], ...
                TASK_LABEL, RUN_MODE_LABEL, dataFile, omnibusFile, RUN_MODE);
        end

        Sd = load(dataFile, 'ElectrodeT');
        if ~isfield(Sd, 'ElectrodeT')
            error('Phase17B:MissingVariables', '%s does not contain ElectrodeT as expected from X17.', dataFile);
        end
        ElectrodeT = Sd.ElectrodeT;
        clear Sd;

        So = load(omnibusFile, 'OmnibusT');
        if ~isfield(So, 'OmnibusT')
            error('Phase17B:MissingVariables', '%s does not contain OmnibusT as expected from X17.', omnibusFile);
        end
        OmnibusT_X17 = So.OmnibusT(strcmp(So.OmnibusT.Level, 'Electrode'), :);
        clear So;

        if isempty(ElectrodeT)
            warning('Task=%s, RUN_MODE=%s: ElectrodeT from X17 is empty -- skipping.', TASK_LABEL, RUN_MODE_LABEL);
            continue;
        end
        ELECTRODES = unique(ElectrodeT.Channel);
        nElec = numel(ELECTRODES);

        fprintf('Loaded X17''s saved ElectrodeT (%d rows, %d electrodes) and electrode-level OmnibusT (%d rows) for Task=%s, RUN_MODE=%s.\n', ...
            height(ElectrodeT), nElec, height(OmnibusT_X17), TASK_LABEL, RUN_MODE_LABEL);

        %% ----------------------------------------------------------------
        %  REFIT all Electrode x Band x Metric two-factor models,
        %  cross-check against X17's own saved omnibus result (all 3
        %  terms), and extract all 3 contrast families (ungated).
        % ----------------------------------------------------------------
        allRows = cell(0, 13);
        crossCheckFails = {};
        nRefitOK = 0; nRefitFailed = 0; nCrossCheckOK = 0;

        TERM_KEYS = {'ConditionMain', 'TimeHalfMain', 'Interaction'};

        for e = 1:nElec
            elecName = ELECTRODES{e};
            for m = 1:nMetrics
                metricName = METRICS{m};
                for b = 1:nBands
                    bandName = ALL_BANDS{b};
                    subT = ElectrodeT(strcmp(ElectrodeT.Channel, elecName) & strcmp(ElectrodeT.Band, bandName), :);
                    fr = local_fit_two_factor_omnibus(subT, metricName, CONDITION_ORDER, TIMEHALF_ORDER);
                    if ~fr.ok
                        nRefitFailed = nRefitFailed + 1;
                        continue;
                    end

                    % --- cross-check all 3 terms against X17's own saved electrode omnibus rows ---
                    mismatchHere = false;
                    for tI = 1:3
                        x17Row = OmnibusT_X17(strcmp(OmnibusT_X17.Unit, elecName) & strcmp(OmnibusT_X17.Metric, metricName) & ...
                            strcmp(OmnibusT_X17.Band, bandName) & strcmp(OmnibusT_X17.Term, TERM_KEYS{tI}), :);
                        if height(x17Row) ~= 1
                            mismatchHere = true;
                            crossCheckFails{end+1} = sprintf('%s/%s/%s/%s: expected exactly 1 matching row in X17''s OmnibusT, found %d.', ...
                                elecName, metricName, bandName, TERM_KEYS{tI}, height(x17Row)); %#ok<AGROW>
                            continue;
                        end
                        fMismatch = abs(fr.F(tI) - x17Row.FStat(1)) > FIT_TOL_F * max(1, abs(x17Row.FStat(1)));
                        pMismatch = abs(fr.p(tI) - x17Row.pValue(1)) > FIT_TOL_P;
                        if fMismatch || pMismatch
                            mismatchHere = true;
                            crossCheckFails{end+1} = sprintf('%s/%s/%s/%s: refit F=%.6g p=%.6g vs X17-saved F=%.6g p=%.6g (mismatch beyond tolerance).', ...
                                elecName, metricName, bandName, TERM_KEYS{tI}, fr.F(tI), fr.p(tI), x17Row.FStat(1), x17Row.pValue(1)); %#ok<AGROW>
                        end
                    end
                    if mismatchHere
                        continue;
                    end
                    nCrossCheckOK = nCrossCheckOK + 1;

                    % --- all 3 contrast families, ungated (VERBATIM from X17's own functions) ---
                    rows3 = local_run_pairwise_condition_marginal(fr.lme, REF_LEVEL, NONREF_TIMEHALF, CONTRASTS, 0.05);
                    for c = 1:3
                        allRows(end+1, :) = [{elecName, metricName, bandName, 'ConditionMarginal'}, rows3(c, :)]; %#ok<AGROW>
                    end
                    row1 = local_run_pairwise_timehalf_marginal(fr.lme, REF_LEVEL, CONDITION_ORDER, REF_TIMEHALF, NONREF_TIMEHALF);
                    allRows(end+1, :) = [{elecName, metricName, bandName, 'TimeHalfMarginal'}, row1]; %#ok<AGROW>
                    rows6 = local_run_pairwise_condition_simple(fr.lme, REF_LEVEL, REF_TIMEHALF, NONREF_TIMEHALF, TIMEHALF_ORDER, CONTRASTS, 0.05);
                    for c = 1:6
                        allRows(end+1, :) = [{elecName, metricName, bandName, 'InteractionSimple'}, rows6(c, :)]; %#ok<AGROW>
                    end
                    nRefitOK = nRefitOK + 1;
                end
            end
            if mod(e, 10) == 0 || e == nElec
                fprintf('  ... %d/%d electrodes processed (Task=%s, RUN_MODE=%s; %d models exported OK, %d cross-check fails so far)\n', ...
                    e, nElec, TASK_LABEL, RUN_MODE_LABEL, nRefitOK, numel(crossCheckFails));
            end
        end

        ElectrodeAllContrastsT_X17 = local_build_table(allRows, [1 2 3 4 5], [6 7 8 9 10 11 12 13], ...
            {'Electrode','Metric','Band','ContrastFamily','Contrast','Estimate','SE','DF','tStat','pValue','pHolm_withinFamily','CI_Low','CI_High'});

        fprintf('\nRefit summary (Task=%s, RUN_MODE=%s): %d models refit and exported OK (x10 rows each), %d refit/skip failures, %d cross-check mismatches.\n', ...
            TASK_LABEL, RUN_MODE_LABEL, nRefitOK, nRefitFailed, numel(crossCheckFails));

        if ~isempty(crossCheckFails)
            grandCrossCheckFails = grandCrossCheckFails + numel(crossCheckFails);
            fprintf('*** CROSS-CHECK MISMATCHES (refit result does not match X17''s own saved result) ***\n');
            for i = 1:numel(crossCheckFails)
                fprintf('  - %s\n', crossCheckFails{i});
            end
            error('Phase17B:CrossCheckFailed', ...
                ['Task=%s, RUN_MODE=%s: %d refit(s) did not match X17''s own saved OmnibusT result ' ...
                 'within tolerance (F: rel %.1e, p: abs %.1e). Stopping before saving any output for this ' ...
                 'Task/RUN_MODE.'], TASK_LABEL, RUN_MODE_LABEL, numel(crossCheckFails), FIT_TOL_F, FIT_TOL_P);
        else
            fprintf('ALL %d refit models cross-checked EXACTLY against X17''s own saved omnibus result (within tolerance, all 3 terms). No mismatches.\n', nCrossCheckOK);
        end

        %% ----------------------------------------------------------------
        %  SAVE -- append-only into X17's own .mat; new .csv
        % ----------------------------------------------------------------
        save(dataFile, 'ElectrodeAllContrastsT_X17', '-append');
        writetable(ElectrodeAllContrastsT_X17, fullfile(STATS_DIR, 'electrode_all_contrasts_ungated.csv'));
        fprintf('\nSaved electrode_all_contrasts_ungated.csv (%d rows) and appended ElectrodeAllContrastsT_X17 into %s.\n', ...
            height(ElectrodeAllContrastsT_X17), dataFile);

        %% ----------------------------------------------------------------
        %  REPORT
        % ----------------------------------------------------------------
        reportLines = {};
        reportLines{end+1} = repmat('=', 1, 60);
        reportLines{end+1} = sprintf('X17B -- ELECTRODE-LEVEL ALL-CONTRASTS UNGATED REPORT (Task=%s, RUN_MODE=%s)', TASK_LABEL, RUN_MODE_LABEL);
        reportLines{end+1} = repmat('=', 1, 60);
        reportLines{end+1} = '';
        reportLines{end+1} = 'Companion to X17_half_split_condition_time_cognitive.m -- does not modify X17 and';
        reportLines{end+1} = 'does not change any FDR/Holm/inferential result. Re-fits the SAME electrode-level';
        reportLines{end+1} = 'two-factor models X17 already fit (from X17''s own saved, already-filtered ElectrodeT),';
        reportLines{end+1} = 'purely to export all 3 contrast families (ConditionMarginal/TimeHalfMarginal/';
        reportLines{end+1} = 'InteractionSimple) for EVERY electrode (ungated, descriptive -- needed for a';
        reportLines{end+1} = 'full-scalp topomap). pHolm_withinFamily is NOT a substitute for X17''s own';
        reportLines{end+1} = 'per-term-type BH-FDR-gated phase17_pairwise_contrasts.';
        reportLines{end+1} = '';
        reportLines{end+1} = sprintf('Models refit and exported: %d / %d expected (%d electrodes x %d bands x %d metrics; 10 rows per model)', ...
            nRefitOK, nElec*nBands*nMetrics, nElec, nBands, nMetrics);
        reportLines{end+1} = sprintf('Refit/skip failures: %d', nRefitFailed);
        reportLines{end+1} = sprintf('Cross-check against X17''s own saved omnibus result (all 3 terms): %d checked, %d mismatches.', nCrossCheckOK, numel(crossCheckFails));
        if isempty(crossCheckFails)
            reportLines{end+1} = 'Every refit model reproduced X17''s own F-statistics and p-values (all 3 terms) within tolerance.';
        end
        reportLines{end+1} = '';
        reportLines{end+1} = 'Output:';
        reportLines{end+1} = ['  ' fullfile(STATS_DIR, 'electrode_all_contrasts_ungated.csv')];
        reportLines{end+1} = ['  ' dataFile ' (ElectrodeAllContrastsT_X17 appended)'];
        reportLines{end+1} = repmat('=', 1, 60);

        reportText = strjoin(reportLines, '\n');
        fprintf('%s\n', reportText);

        fid1 = fopen(fullfile(REPORT_DIR, sprintf('phase17B_electrode_all_contrasts_report_%s_%s.txt', TASK_LABEL, RUN_MODE_LABEL)), 'w');
        fprintf(fid1, '%s\n', reportText);
        fclose(fid1);
    end
end

fprintf('\n=== X17B complete for all RUN_MODEs x both Tasks (Nback, SART -- Full out of scope). ===\n');
if grandCrossCheckFails == 0
    fprintf('No Task/RUN_MODE combination had any cross-check mismatch.\n');
end

%% ======================= LOCAL FUNCTIONS (all after top-level code) =======================
%  local_fit_two_factor_omnibus, local_find_interaction_coeff_idx,
%  local_marginal_condition_contrast, local_simple_condition_contrast,
%  local_marginal_timehalf_contrast, local_run_pairwise_condition_marginal,
%  local_run_pairwise_condition_simple, local_run_pairwise_timehalf_marginal,
%  local_condition_contrast, local_holm_bonferroni, local_subject_numeric,
%  local_build_table are ALL reproduced VERBATIM from
%  X17_half_split_condition_time_cognitive.m -- this companion script
%  deliberately does not re-derive this contrast arithmetic independently.

function n = local_subject_numeric(subjVal)
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
    cRow = zeros(1, numel(coeffNames));
    if ~strcmp(level1, refLevel)
        name1 = sprintf('%s_%s', varName, level1);
        idx1 = find(strcmp(coeffNames, name1));
        if isempty(idx1)
            error('Phase17B:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name1);
        end
        cRow(idx1) = cRow(idx1) - 1;
    end
    if ~strcmp(level2, refLevel)
        name2 = sprintf('%s_%s', varName, level2);
        idx2 = find(strcmp(coeffNames, name2));
        if isempty(idx2)
            error('Phase17B:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name2);
        end
        cRow(idx2) = cRow(idx2) + 1;
    end
end

function T = local_build_table(rows, textColIdx, numColIdx, varNames)
    nCols = numel(varNames);
    if size(rows, 2) ~= nCols
        error('Phase17B:BuildTableColumnMismatch', ...
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
            error('Phase17B:BuildTableColumnMismatch', 'Column %d not classified as text or numeric.', j);
        end
    end
    T = table(cols{:}, 'VariableNames', varNames);
end

function fr = local_fit_two_factor_omnibus(subT, metricName, conditionOrder, timeHalfOrder)
% Reproduced VERBATIM from X17_half_split_condition_time_cognitive.m.
    fr = struct('ok', false, 'lme', [], 'F', nan(1,3), 'DF1', nan(1,3), 'DF2', nan(1,3), 'p', nan(1,3), ...
        'warnMsg', '', 'errMsg', '');
    if isempty(subT) || ~ismember(metricName, subT.Properties.VariableNames)
        fr.errMsg = 'Empty subT or missing metric column.';
        return;
    end
    subT.Subject = categorical(local_subject_numeric(subT.Subject));
    subT.Condition = categorical(cellstr(subT.Condition), conditionOrder);
    subT.TimeHalf = categorical(cellstr(subT.TimeHalf), timeHalfOrder);
    formula = sprintf('%s ~ Condition*TimeHalf + (1|Subject)', metricName);
    try
        lastwarn('');
        lme = fitlme(subT, formula, 'FitMethod', 'REML');
        [warnMsg, ~] = lastwarn();
        fr.warnMsg = warnMsg;
    catch ME
        fr.errMsg = ME.message;
        return;
    end
    try
        aov = anova(lme, 'DFMethod', 'satterthwaite');
    catch ME
        fr.errMsg = sprintf('anova(lme) failed: %s', ME.message);
        return;
    end
    termNames = {'Condition', 'TimeHalf', 'Condition:TimeHalf'};
    okAll = true;
    Fv = nan(1,3); DF1v = nan(1,3); DF2v = nan(1,3); pv = nan(1,3);
    for tI = 1:3
        rowIdx = find(strcmp(aov.Term, termNames{tI}), 1);
        if isempty(rowIdx)
            okAll = false;
            fr.errMsg = sprintf('anova(lme) has no "%s" row.', termNames{tI});
            break;
        end
        Fv(tI) = aov.FStat(rowIdx);
        DF1v(tI) = aov.DF1(rowIdx);
        DF2v(tI) = aov.DF2(rowIdx);
        pv(tI) = aov.pValue(rowIdx);
    end
    if ~okAll
        return;
    end
    fr.ok = true;
    fr.lme = lme;
    fr.F = Fv; fr.DF1 = DF1v; fr.DF2 = DF2v; fr.p = pv;
end

function idx = local_find_interaction_coeff_idx(coeffNames, condDummyName, timeHalfDummyName)
% Reproduced VERBATIM from X17_half_split_condition_time_cognitive.m.
    name1 = sprintf('%s:%s', condDummyName, timeHalfDummyName);
    name2 = sprintf('%s:%s', timeHalfDummyName, condDummyName);
    idx = find(strcmp(coeffNames, name1), 1);
    if isempty(idx)
        idx = find(strcmp(coeffNames, name2), 1);
    end
    if isempty(idx)
        error('Phase17B:InteractionCoeffNotFound', ...
            'Could not find interaction coefficient "%s" or "%s" among the model''s CoefficientNames: %s.', ...
            name1, name2, strjoin(coeffNames, ', '));
    end
end

function cRow = local_marginal_condition_contrast(coeffNames, refCond, nonRefTimeHalf, level1, level2)
% Reproduced VERBATIM from X17_half_split_condition_time_cognitive.m.
    cRow = local_condition_contrast(coeffNames, 'Condition', refCond, level1, level2);
    timeHalfDummyName = sprintf('TimeHalf_%s', nonRefTimeHalf);
    if ~strcmp(level1, refCond)
        condDummyName1 = sprintf('Condition_%s', level1);
        idx1 = local_find_interaction_coeff_idx(coeffNames, condDummyName1, timeHalfDummyName);
        cRow(idx1) = cRow(idx1) - 0.5;
    end
    if ~strcmp(level2, refCond)
        condDummyName2 = sprintf('Condition_%s', level2);
        idx2 = local_find_interaction_coeff_idx(coeffNames, condDummyName2, timeHalfDummyName);
        cRow(idx2) = cRow(idx2) + 0.5;
    end
end

function cRow = local_simple_condition_contrast(coeffNames, refCond, refTimeHalf, nonRefTimeHalf, timeHalfLevel, level1, level2)
% Reproduced VERBATIM from X17_half_split_condition_time_cognitive.m.
    cRow = local_condition_contrast(coeffNames, 'Condition', refCond, level1, level2);
    if strcmp(timeHalfLevel, refTimeHalf)
        return;
    end
    timeHalfDummyName = sprintf('TimeHalf_%s', nonRefTimeHalf);
    if ~strcmp(level1, refCond)
        condDummyName1 = sprintf('Condition_%s', level1);
        idx1 = local_find_interaction_coeff_idx(coeffNames, condDummyName1, timeHalfDummyName);
        cRow(idx1) = cRow(idx1) - 1;
    end
    if ~strcmp(level2, refCond)
        condDummyName2 = sprintf('Condition_%s', level2);
        idx2 = local_find_interaction_coeff_idx(coeffNames, condDummyName2, timeHalfDummyName);
        cRow(idx2) = cRow(idx2) + 1;
    end
end

function cRow = local_marginal_timehalf_contrast(coeffNames, refCond, condLevels, refTimeHalf, nonRefTimeHalf)
% Reproduced VERBATIM from X17_half_split_condition_time_cognitive.m.
    cRow = zeros(1, numel(coeffNames));
    timeHalfDummyName = sprintf('TimeHalf_%s', nonRefTimeHalf);
    idxTH = find(strcmp(coeffNames, timeHalfDummyName), 1);
    if isempty(idxTH)
        error('Phase17B:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', timeHalfDummyName);
    end
    cRow(idxTH) = cRow(idxTH) + 1;
    nCondLevels = numel(condLevels);
    for k = 1:nCondLevels
        lvl = condLevels{k};
        if strcmp(lvl, refCond); continue; end
        condDummyName = sprintf('Condition_%s', lvl);
        idxInt = local_find_interaction_coeff_idx(coeffNames, condDummyName, timeHalfDummyName);
        cRow(idxInt) = cRow(idxInt) + (1 / nCondLevels);
    end
end

function rows3 = local_run_pairwise_condition_marginal(lme, refCond, nonRefTimeHalf, CONTRASTS, alphaFDR)
% Reproduced VERBATIM from X17_half_split_condition_time_cognitive.m.
    coeffNames = lme.CoefficientNames;
    beta = lme.Coefficients.Estimate;
    covB = lme.CoefficientCovariance;
    contrastP = nan(3, 1);
    contrastData = cell(3, 8);
    for c = 1:3
        lvl2 = CONTRASTS{c, 1}; lvl1 = CONTRASTS{c, 2}; label = CONTRASTS{c, 3};
        C = local_marginal_condition_contrast(coeffNames, refCond, nonRefTimeHalf, lvl1, lvl2);
        est = C * beta;
        se = sqrt(C * covB * C');
        [p, ~, ~, df2c] = coefTest(lme, C, 0, 'DFMethod', 'satterthwaite');
        tstat = est / se;
        tcrit = tinv(0.975, df2c);
        ciLow = est - tcrit * se;
        ciHigh = est + tcrit * se;
        contrastP(c) = p;
        contrastData(c, :) = {label, est, se, df2c, tstat, p, ciLow, ciHigh};
    end
    [adjP, ~] = local_holm_bonferroni(contrastP, alphaFDR);
    rows3 = cell(3, 9);
    for c = 1:3
        rows3(c, :) = [contrastData(c, 1:6), {adjP(c)}, contrastData(c, 7:8)];
    end
end

function rows6 = local_run_pairwise_condition_simple(lme, refCond, refTimeHalf, nonRefTimeHalf, timeHalfOrder, CONTRASTS, alphaFDR)
% Reproduced VERBATIM from X17_half_split_condition_time_cognitive.m.
    coeffNames = lme.CoefficientNames;
    beta = lme.Coefficients.Estimate;
    covB = lme.CoefficientCovariance;
    contrastP = nan(6, 1);
    contrastData = cell(6, 8);
    rowIdx = 0;
    for h = 1:numel(timeHalfOrder)
        timeHalfLevel = timeHalfOrder{h};
        for c = 1:3
            rowIdx = rowIdx + 1;
            lvl2 = CONTRASTS{c, 1}; lvl1 = CONTRASTS{c, 2};
            label = sprintf('%s @ TimeHalf=%s', CONTRASTS{c, 3}, timeHalfLevel);
            C = local_simple_condition_contrast(coeffNames, refCond, refTimeHalf, nonRefTimeHalf, timeHalfLevel, lvl1, lvl2);
            est = C * beta;
            se = sqrt(C * covB * C');
            [p, ~, ~, df2c] = coefTest(lme, C, 0, 'DFMethod', 'satterthwaite');
            tstat = est / se;
            tcrit = tinv(0.975, df2c);
            ciLow = est - tcrit * se;
            ciHigh = est + tcrit * se;
            contrastP(rowIdx) = p;
            contrastData(rowIdx, :) = {label, est, se, df2c, tstat, p, ciLow, ciHigh};
        end
    end
    [adjP, ~] = local_holm_bonferroni(contrastP, alphaFDR);
    rows6 = cell(6, 9);
    for r = 1:6
        rows6(r, :) = [contrastData(r, 1:6), {adjP(r)}, contrastData(r, 7:8)];
    end
end

function row1 = local_run_pairwise_timehalf_marginal(lme, refCond, condLevels, refTimeHalf, nonRefTimeHalf)
% Reproduced VERBATIM from X17_half_split_condition_time_cognitive.m.
    coeffNames = lme.CoefficientNames;
    beta = lme.Coefficients.Estimate;
    covB = lme.CoefficientCovariance;
    C = local_marginal_timehalf_contrast(coeffNames, refCond, condLevels, refTimeHalf, nonRefTimeHalf);
    est = C * beta;
    se = sqrt(C * covB * C');
    [p, ~, ~, df2c] = coefTest(lme, C, 0, 'DFMethod', 'satterthwaite');
    tstat = est / se;
    tcrit = tinv(0.975, df2c);
    ciLow = est - tcrit * se;
    ciHigh = est + tcrit * se;
    label = sprintf('TimeHalf=%s - TimeHalf=%s', nonRefTimeHalf, refTimeHalf);
    row1 = {label, est, se, df2c, tstat, p, p, ciLow, ciHigh};
end
