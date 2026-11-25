%% ========================================================================
% X04_aggregate_cognitive_subject_level.m
%
% Purpose: convert X03's window-level spectral results into one
% analysis-ready observation per (Subject x TaskType x BlockNumber)
% recording -- the inferential unit for every later statistical phase,
% computed separately per Task -- while preserving channel-level and
% frequency-resolved information, without ever treating individual
% windows as independent observations.
%
% Method: arithmetic-mean window aggregation; LogAbsolutePower =
% ln(mean(AbsolutePower)) (never mean(ln(...))); arithmetic-mean
% global-scalp aggregation with a median-across-channels secondary
% cross-check; channels are never averaged away at the per-channel
% product level. TaskType and BlockNumber are carried as key columns
% through every structure (A: frequency-resolved, B: per-channel band
% power, C: whole-scalp band power, D: recording QC); the
% aggregation/grouping key is Subject x TaskType x BlockNumber x
% Condition. Diagnostic 5 (extreme-value flagging, median +/- 3 robust
% MADs) is computed separately per task type, pooling all conditions
% (including Baseline) within that task type -- Full/N-back/SART are
% structurally different recording types (different durations, task
% demands, and plausibly different baseline power scales), so a single
% MAD computed across all three would be dominated by whichever task has
% the most spread; the rule itself (median +/- 3*1.4826*MAD, flag only,
% never exclude) is unchanged. Integrity checks (unique-key integrity, no
% duplicate Subject x Task x Block x Condition rows, no unexpected
% NaN/Inf) compute their own expected counts from the data rather than
% against a hard-coded literal, and the resulting recording counts are
% reported explicitly in the report/console for review.
%
% This script does not compute PSD (reads X03's already-computed output;
% no .set files touched, no EEGLAB/Signal Processing Toolbox needed),
% perform inferential statistics of any kind, compare E vs N vs R or
% Full vs N-back vs SART or model a Task x Condition interaction,
% redefine frequency bands or the relative-power denominator (both
% loaded from X03's own saved parameters), impute or fabricate data for
% any missing Subject x Task x Block combination, or exclude any
% recording or subject based on QC/variability metrics (all such metrics
% are reported, never acted on automatically).
%
% Requirements: base MATLAB only (no EEGLAB, no Signal Processing
% Toolbox, no Statistics and Machine Learning Toolbox); MATLAB R2016b or
% later (local functions in a script file). X03 outputs must already
% exist: cognitive_psd.mat and phase2_processing_log.mat.
%
% Inputs: X03's PSD output.
%
% Outputs (nothing under X03's output is ever written to):
%   <ANALYSIS_ROOT>/04_subject_level/
%       phase3_subject_level.mat            (Structure A: frequency-resolved)
%       phase3_subject_level_bandpower.mat  (Structures B & C: band power)
%       phase3_recording_summary.xlsx/.mat  (Structure D: recording QC)
%       phase3_parameters.mat
%       diagnostic1_subject_bandpower_distributions.png
%       diagnostic2_mean_vs_median_window_aggregation.png
%       diagnostic3_valid_windows_by_condition.png
%       diagnostic4_alpha_power_by_condition.png
%       diagnostic5_extreme_value_flags.csv   (includes a TaskType column)
%   <ANALYSIS_ROOT>/09_reports/
%       phase3_aggregation_report_cognitive.txt
%
% RUN_MODE: not applicable; this script runs once.
% ========================================================================

clear; clc;

%% ------------------------------------------------------------------
%  0. USER-EDITABLE SETTINGS
% ------------------------------------------------------------------

config;
DATA_ROOT_FOR_SIBLING = DATA_DIR;
[DATA_ROOT_PARENT, ~] = fileparts(DATA_ROOT_FOR_SIBLING);
ANALYSIS_ROOT = OUT_DIR;

PSD_DIR     = fullfile(ANALYSIS_ROOT, '03_psd');
SUBJ_DIR    = fullfile(ANALYSIS_ROOT, '04_subject_level');
REPORT_DIR  = fullfile(ANALYSIS_ROOT, '09_reports');
if ~exist(SUBJ_DIR, 'dir');   mkdir(SUBJ_DIR);   end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

PSD_MAT_FILE = fullfile(PSD_DIR, 'cognitive_psd.mat');
LOG_MAT_FILE = fullfile(PSD_DIR, 'phase2_processing_log.mat');

fprintf('========================================================\n');
fprintf('PHASE 3 -- COGNITIVE-TASK SUBJECT-LEVEL AGGREGATION\n');
fprintf('========================================================\n');
fprintf('ANALYSIS_ROOT : %s\n', ANALYSIS_ROOT);
fprintf('OUTPUT        : %s\n', SUBJ_DIR);
fprintf('--------------------------------------------------------\n\n');

if ~exist(PSD_MAT_FILE, 'file')
    error(['X03 PSD file not found: %s\nRun X03_compute_cognitive_psd.m first.'], PSD_MAT_FILE);
end
if ~exist(LOG_MAT_FILE, 'file')
    error(['X03 processing log not found: %s\nRun X03_compute_cognitive_psd.m first.'], LOG_MAT_FILE);
end

%% ------------------------------------------------------------------
%  1. LOAD PHASE 2 OUTPUT (source of truth; PSD is never recomputed)
% ------------------------------------------------------------------

psdFile = load(PSD_MAT_FILE, 'psd_data', 'FrequencyVectorHz', 'ChannelLabelsReference', 'PARAMS');
psd_data   = psdFile.psd_data;
freqVecRef = psdFile.FrequencyVectorHz(:);
chanLabelsRef = psdFile.ChannelLabelsReference(:);
PARAMS2 = psdFile.PARAMS;

logFile = load(LOG_MAT_FILE, 'ProcessingLogT');
ProcessingLogT = logFile.ProcessingLogT;

nRec = numel(psd_data);
bandNames = fieldnames(PARAMS2.BandDefinitionsHz);
nBands = numel(bandNames);
nChanRef = numel(chanLabelsRef);
nFreqRef = numel(freqVecRef);

fprintf('Loaded Phase 2 output: %d recordings, %d bands, %d channels, %d frequency bins.\n\n', ...
    nRec, nBands, nChanRef, nFreqRef);

%% ------------------------------------------------------------------
%  2. PHASE 3 AGGREGATION PARAMETERS (methodological choices, saved
%     verbatim for reproducibility -- IDENTICAL values/rationale to
%     04_aggregate_subject_level.m)
% ------------------------------------------------------------------

PARAMS3 = struct();
PARAMS3.ScriptName  = 'X04_aggregate_cognitive_subject_level.m';
PARAMS3.GeneratedOn = datestr(now);
PARAMS3.SourcePhase2File = PSD_MAT_FILE;

PARAMS3.WindowAggregationStatistic = 'arithmetic_mean';
PARAMS3.AbsolutePowerLogTransform = true;
PARAMS3.LogTransformConvention = 'natural log (ln), applied to the recording-level AbsolutePower AFTER window-averaging';
PARAMS3.RelativePowerLogTransform = false;
PARAMS3.RelativePowerLogTransformRationale = ['Relative power is a bounded proportion, typically far less skewed ' ...
    'than absolute power. No log-transform is applied here (unchanged from the intervention pipeline).'];
PARAMS3.GlobalScalpAggregation = 'arithmetic_mean_across_channels';
PARAMS3.GlobalScalpSecondary   = 'median_across_channels (saved alongside, not primary)';
PARAMS3.ChannelHandling = 'Channels are NEVER averaged away at the per-channel product level (Structures A and B). A separate whole-scalp summary (Structure C) is derived FROM the per-channel values, not a replacement for them.';
PARAMS3.BandDefinitionsHz = PARAMS2.BandDefinitionsHz;
PARAMS3.PrimaryBands = PARAMS2.PrimaryBands;
PARAMS3.ExploratoryBands = PARAMS2.ExploratoryBands;
PARAMS3.RelativePowerDenominatorRangeHz = PARAMS2.RelativePowerDenominatorRangeHz;
PARAMS3.NotchGuardBandExcludedHz = PARAMS2.NotchGuardBandExcludedHz;
PARAMS3.BandEdgeConvention = PARAMS2.BandEdgeConvention;
PARAMS3.ValidWindowHandling = ['Only windows with WindowValid=true (per Phase 2 QC) contribute to any ' ...
    'aggregation. Each recording contributes exactly ONE row per channel per band (and one whole-scalp ' ...
    'row per band) regardless of how many valid windows it had.'];
PARAMS3.ExtremeValueRule = 'median +/- 3*1.4826*MAD, computed per band, SEPARATELY PER TASK TYPE, pooling all conditions within that task type (never per-condition); flag only, never exclude. (Adapted from the intervention pipeline''s single-task pooling rule -- see script header item 2.)';
PARAMS3.TaskDimension = {'full', 'nback', 'sart'};
PARAMS3.ConditionsIncluded = {'Baseline', 'E', 'N', 'R'};
PARAMS3.BaselinePolicy = 'Baseline recordings are aggregated here exactly like E/N/R (no special handling in this script). Exclusion from the inferential model happens only in X05/X06/X08/X09.';

save(fullfile(SUBJ_DIR, 'phase3_parameters.mat'), 'PARAMS3');

fprintf('Window-aggregation statistic : %s\n', PARAMS3.WindowAggregationStatistic);
fprintf('Absolute-power log transform : %d (%s)\n', PARAMS3.AbsolutePowerLogTransform, PARAMS3.LogTransformConvention);
fprintf('Global scalp aggregation     : %s\n\n', PARAMS3.GlobalScalpAggregation);

%% ------------------------------------------------------------------
%  3. MAIN LOOP: AGGREGATE EACH RECORDING
% ------------------------------------------------------------------

subject_psd = struct([]);

bandPowerCols = {'Subject', 'TaskType', 'BlockNumber', 'Condition', 'Channel', 'Band', ...
    'AbsolutePower', 'LogAbsolutePower', 'RelativePower', ...
    'WindowMedian_Absolute', 'WindowSD_Absolute', 'WindowIQR_Absolute', 'WindowSkewness_Absolute', ...
    'WindowMedian_Relative', 'WindowSD_Relative', 'WindowIQR_Relative', 'WindowSkewness_Relative'};
bandPowerRows = cell(0, numel(bandPowerCols));

globalCols = {'Subject', 'TaskType', 'BlockNumber', 'Condition', 'Band', ...
    'AbsolutePower', 'LogAbsolutePower', 'RelativePower', ...
    'AbsolutePower_MedianAcrossChannels', 'RelativePower_MedianAcrossChannels'};
globalRows = cell(0, numel(globalCols));

qcCols = {'Subject', 'TaskType', 'BlockNumber', 'Condition', 'ValidWindows', 'RejectedWindows', ...
    'CandidateWindows', 'PercentValid', 'DurationSec', ...
    'RejectedNaNInf', 'RejectedAmplitude', 'RejectedFlat'};
qcRows = cell(0, numel(qcCols));

integrityIssues = {};

for i = 1:nRec
    rec = psd_data(i);
    subjStr = rec.Subject;
    taskStr = rec.TaskType;
    blockNum = rec.BlockNumber;
    condStr = rec.Condition;

    fprintf('[%d/%d] Subject %s, Task %s, Block %g, Condition %s ... ', i, nRec, subjStr, taskStr, blockNum, condStr);

    if numel(rec.ChannelLabels) ~= nChanRef || ~isequal(rec.ChannelLabels, chanLabelsRef)
        integrityIssues{end+1} = sprintf('Sub%s Task%s Block%g: channel labels differ from the reference set.', subjStr, taskStr, blockNum); %#ok<AGROW>
    end
    if numel(rec.FrequenciesHz) ~= nFreqRef || max(abs(rec.FrequenciesHz - freqVecRef)) > 1e-9
        integrityIssues{end+1} = sprintf('Sub%s Task%s Block%g: frequency vector differs from the reference vector.', subjStr, taskStr, blockNum); %#ok<AGROW>
    end

    valid = rec.WindowValid;
    nValidThis = sum(valid);

    logIdx = find(strcmp(ProcessingLogT.Subject, subjStr) & strcmp(ProcessingLogT.TaskType, taskStr) & ProcessingLogT.BlockNumber == blockNum);
    if numel(logIdx) ~= 1
        integrityIssues{end+1} = sprintf('Sub%s Task%s Block%g: expected exactly 1 matching Phase 2 processing-log row, found %d.', ...
            subjStr, taskStr, blockNum, numel(logIdx)); %#ok<AGROW>
        candW = rec.NumCandidateWindows; validW = rec.NumValidWindows; pctV = NaN; durS = NaN;
        rNaN = NaN; rAmp = NaN; rFlat = NaN;
    else
        lr = ProcessingLogT(logIdx, :);
        candW = lr.NumCandidateWindows; validW = lr.NumValidWindows; pctV = lr.PercentValid; durS = lr.DurationSec;
        rNaN = lr.RejectedNaNInf; rAmp = lr.RejectedAmplitude; rFlat = lr.RejectedFlat;
        if iscell(candW); candW = candW{1}; end
        if iscell(validW); validW = validW{1}; end
        if iscell(pctV); pctV = pctV{1}; end
        if iscell(durS); durS = durS{1}; end
        if iscell(rNaN); rNaN = rNaN{1}; end
        if iscell(rAmp); rAmp = rAmp{1}; end
        if iscell(rFlat); rFlat = rFlat{1}; end
    end
    qcRows(end+1, :) = {subjStr, taskStr, blockNum, condStr, validW, candW - validW, candW, pctV, durS, rNaN, rAmp, rFlat}; %#ok<AGROW>

    if nValidThis == 0
        integrityIssues{end+1} = sprintf('Sub%s Task%s Block%g: 0 valid windows -- cannot aggregate into Structures A/B/C; recorded in the QC table only.', subjStr, taskStr, blockNum); %#ok<AGROW>
        fprintf('SKIPPED for aggregation (0 valid windows; QC row still recorded)\n');
        continue;
    end

    aggAbsPSD = rec.MeanPSD;
    meanTotalPower = mean(rec.TotalPower_1to45Hz(valid, :), 1);
    aggRelPSD = aggAbsPSD ./ meanTotalPower(:);

    entryA = struct();
    entryA.Subject = subjStr;
    entryA.TaskType = taskStr;
    entryA.BlockNumber = blockNum;
    entryA.Condition = condStr;
    entryA.ChannelLabels = chanLabelsRef;
    entryA.FrequenciesHz = freqVecRef;
    entryA.AggregatedAbsolutePSD = aggAbsPSD;
    entryA.AggregatedRelativePSD = aggRelPSD;
    entryA.NumValidWindows = nValidThis;
    subject_psd = concat_struct_array(subject_psd, entryA);

    for b = 1:nBands
        band = bandNames{b};
        absVals = rec.AbsoluteBandPower.(band)(valid, :);
        relVals = rec.RelativeBandPower.(band)(valid, :);

        absMean = mean(absVals, 1);
        relMean = mean(relVals, 1);
        absLog  = log(absMean);
        if any(absMean <= 0)
            integrityIssues{end+1} = sprintf('Sub%s Task%s Block%g Band %s: non-positive AbsolutePower found in %d channel(s) -- LogAbsolutePower set to NaN there, not silently patched.', ...
                subjStr, taskStr, blockNum, band, sum(absMean <= 0)); %#ok<AGROW>
            absLog(absMean <= 0) = NaN;
        end

        absMedian = median(absVals, 1);
        absSD     = std(absVals, 0, 1);
        absIQR    = local_iqr(absVals);
        absSkew   = local_skewness(absVals);

        relMedian = median(relVals, 1);
        relSD     = std(relVals, 0, 1);
        relIQR    = local_iqr(relVals);
        relSkew   = local_skewness(relVals);

        for ch = 1:nChanRef
            bandPowerRows(end+1, :) = {subjStr, taskStr, blockNum, condStr, chanLabelsRef{ch}, band, ...
                absMean(ch), absLog(ch), relMean(ch), ...
                absMedian(ch), absSD(ch), absIQR(ch), absSkew(ch), ...
                relMedian(ch), relSD(ch), relIQR(ch), relSkew(ch)}; %#ok<AGROW>
        end

        globalAbs = mean(absMean);
        globalRel = mean(relMean);
        globalLog = log(globalAbs);
        if globalAbs <= 0
            integrityIssues{end+1} = sprintf('Sub%s Task%s Block%g Band %s: non-positive whole-scalp AbsolutePower -- LogAbsolutePower set to NaN.', subjStr, taskStr, blockNum, band); %#ok<AGROW>
            globalLog = NaN;
        end
        globalAbsMedian = median(absMean);
        globalRelMedian = median(relMean);

        globalRows(end+1, :) = {subjStr, taskStr, blockNum, condStr, band, ...
            globalAbs, globalLog, globalRel, globalAbsMedian, globalRelMedian}; %#ok<AGROW>
    end

    fprintf('OK\n');
end

fprintf('\nMain loop complete: %d/%d recordings aggregated.\n\n', numel(subject_psd), nRec);

%% ------------------------------------------------------------------
%  4. ASSEMBLE TABLES
% ------------------------------------------------------------------

BandPowerT = cell2table(bandPowerRows, 'VariableNames', bandPowerCols);
BandPowerT = normalize_text_columns(BandPowerT);

GlobalBandPowerT = cell2table(globalRows, 'VariableNames', globalCols);
GlobalBandPowerT = normalize_text_columns(GlobalBandPowerT);

RecordingQCT = cell2table(qcRows, 'VariableNames', qcCols);
RecordingQCT = normalize_text_columns(RecordingQCT);

%% ------------------------------------------------------------------
%  5. DATA INTEGRITY CHECKS (reported, never silently corrected -- see
%     header item 3 for why exact recording-count literals are NOT
%     hard-coded here)
% ------------------------------------------------------------------

fprintf('--- Data integrity checks ---\n');

% (a) Subject x TaskType x BlockNumber unique
subjKey = strcat(RecordingQCT.Subject, '_', RecordingQCT.TaskType, '_blk', string(RecordingQCT.BlockNumber));
if numel(unique(subjKey)) ~= height(RecordingQCT)
    integrityIssues{end+1} = 'DUPLICATE Subject x TaskType x BlockNumber combinations found in the subject-level output.';
end

% (b) exactly one Condition per recording
for i = 1:numel(subject_psd)
    if isempty(subject_psd(i).Condition) || ~ischar(subject_psd(i).Condition)
        integrityIssues{end+1} = sprintf('Sub%s Task%s Block%g: Condition is empty or non-char.', subject_psd(i).Subject, subject_psd(i).TaskType, subject_psd(i).BlockNumber); %#ok<AGROW>
    end
end

% (c) no subject x task receives two different conditions for the SAME block
[uKey, ~, uIdx] = unique(subjKey);
for k = 1:numel(uKey)
    condsHere = unique(RecordingQCT.Condition(uIdx == k));
    if numel(condsHere) > 1
        integrityIssues{end+1} = sprintf('%s: multiple different Condition values found (%s) for the same Subject x TaskType x Block.', ...
            uKey{k}, strjoin(condsHere, ', ')); %#ok<AGROW>
    end
end

% (d) recording count is REPORTED (not hard-gated against a literal --
%     see header item 3) and cross-checked for internal consistency: every
%     TaskType should have the SAME number of subject-level recordings as
%     every other TaskType (all three task types are recorded from the
%     same blocks/subjects, just a different processing stream).
nSubjLevel = numel(subject_psd);
taskCountsHere = struct();
taskTypesObserved = unique(RecordingQCT.TaskType);
for t = 1:numel(taskTypesObserved)
    taskCountsHere.(taskTypesObserved{t}) = sum(strcmp(RecordingQCT.TaskType, taskTypesObserved{t}));
end
taskCountVals = struct2array(taskCountsHere);
if numel(unique(taskCountVals)) > 1
    integrityIssues{end+1} = sprintf('Task types do not have equal subject-level recording counts (%s) -- review for missing/extra recordings.', ...
        strjoin(cellfun(@(t) sprintf('%s=%d', t, taskCountsHere.(t)), taskTypesObserved, 'UniformOutput', false), ', '));
end

% (e) condition counts per TaskType are REPORTED (not gated against a
%     literal -- see header item 3), but checked for a basic structural
%     invariant: Baseline should occur exactly once per Subject x
%     TaskType (Block 1 only), and each of E/N/R should occur at most
%     once per Subject x TaskType (a subject sees each intervention
%     condition at most once per task type, matching the incomplete-block
%     design already established for the intervention data).
for t = 1:numel(taskTypesObserved)
    tt = taskTypesObserved{t};
    subT = RecordingQCT(strcmp(RecordingQCT.TaskType, tt), :);
    subjHere = unique(subT.Subject);
    for si = 1:numel(subjHere)
        condsThisSubjTask = subT.Condition(strcmp(subT.Subject, subjHere{si}));
        nBaseline = sum(strcmp(condsThisSubjTask, 'Baseline'));
        if nBaseline > 1
            integrityIssues{end+1} = sprintf('Sub%s Task%s: %d Baseline recordings found (expected at most 1).', subjHere{si}, tt, nBaseline); %#ok<AGROW>
        end
        for cc = {'E','N','R'}
            nThis = sum(strcmp(condsThisSubjTask, cc{1}));
            if nThis > 1
                integrityIssues{end+1} = sprintf('Sub%s Task%s: %d recordings of condition %s found (expected at most 1).', subjHere{si}, tt, nThis, cc{1}); %#ok<AGROW>
            end
        end
    end
end

% (f) no unexpected NaN/Inf in the core aggregated fields
nNaNInf_Abs = sum(~isfinite(BandPowerT.AbsolutePower));
nNaNInf_Rel = sum(~isfinite(BandPowerT.RelativePower));
nNaNInf_Log = sum(~isfinite(BandPowerT.LogAbsolutePower));
if nNaNInf_Abs > 0
    integrityIssues{end+1} = sprintf('%d non-finite AbsolutePower values found in BandPowerT.', nNaNInf_Abs);
end
if nNaNInf_Rel > 0
    integrityIssues{end+1} = sprintf('%d non-finite RelativePower values found in BandPowerT.', nNaNInf_Rel);
end

% (g) valid-window counts preserved
if abs(sum(RecordingQCT.ValidWindows) - sum(ProcessingLogT.NumValidWindows)) > 0
    integrityIssues{end+1} = 'Sum of ValidWindows in Phase 3 output does not match the sum in the Phase 2 processing log.';
end

for k = 1:numel(integrityIssues)
    fprintf('  ISSUE: %s\n', integrityIssues{k});
end
if isempty(integrityIssues)
    fprintf('  All integrity checks passed.\n');
end
fprintf('\n');

%% ------------------------------------------------------------------
%  6. SAVE OUTPUTS
% ------------------------------------------------------------------

save(fullfile(SUBJ_DIR, 'phase3_subject_level.mat'), ...
    'subject_psd', 'freqVecRef', 'chanLabelsRef', 'PARAMS3', '-v7.3');

save(fullfile(SUBJ_DIR, 'phase3_subject_level_bandpower.mat'), ...
    'BandPowerT', 'GlobalBandPowerT', 'PARAMS3');

writetable(RecordingQCT, fullfile(SUBJ_DIR, 'phase3_recording_summary.xlsx'));
save(fullfile(SUBJ_DIR, 'phase3_recording_summary.mat'), 'RecordingQCT');

%% ------------------------------------------------------------------
%  7. DIAGNOSTICS
% ------------------------------------------------------------------

try
    fig1 = figure('Visible', 'off', 'Position', [100 100 1400 700]);
    nPrimary = numel(PARAMS3.PrimaryBands);
    for b = 1:nPrimary
        subplot(2, ceil(nPrimary/2), b);
        band = PARAMS3.PrimaryBands{b};
        vals = GlobalBandPowerT.AbsolutePower(strcmp(GlobalBandPowerT.Band, band));
        histogram(vals);
        title(sprintf('%s (whole-scalp AbsolutePower, all tasks pooled)', band), 'Interpreter', 'none');
        xlabel('uV^2'); ylabel('Recordings');
    end
    sgtitle('Diagnostic 1: Distribution of subject-level band power (primary bands, all tasks pooled)');
    saveas(fig1, fullfile(SUBJ_DIR, 'diagnostic1_subject_bandpower_distributions.png'));
    close(fig1);
catch ME1
    warning('Could not generate Diagnostic 1 figure: %s', ME1.message);
end

try
    fig2 = figure('Visible', 'off');
    scatter(GlobalBandPowerT.AbsolutePower_MedianAcrossChannels, GlobalBandPowerT.AbsolutePower, 15, 'filled', 'MarkerFaceAlpha', 0.4);
    hold on;
    lims = [0, max([GlobalBandPowerT.AbsolutePower; GlobalBandPowerT.AbsolutePower_MedianAcrossChannels])*1.05];
    plot(lims, lims, 'k--');
    hold off;
    xlabel('Median-across-channels AbsolutePower (uV^2)');
    ylabel('Mean-across-channels AbsolutePower (uV^2)');
    title('Diagnostic 2: Mean vs median aggregation (whole-scalp, all recordings x bands, all tasks pooled)');
    axis(repmat(lims,1,2));
    saveas(fig2, fullfile(SUBJ_DIR, 'diagnostic2_mean_vs_median_window_aggregation.png'));
    close(fig2);
catch ME2
    warning('Could not generate Diagnostic 2 figure: %s', ME2.message);
end

try
    fig3 = figure('Visible', 'off', 'Position', [100 100 1200 400]);
    conds = {'Baseline', 'E', 'N', 'R'};
    taskTypesToPlot = {'full', 'nback', 'sart'};
    for tI = 1:numel(taskTypesToPlot)
        subplot(1, numel(taskTypesToPlot), tI);
        hold on;
        for c = 1:numel(conds)
            vals = RecordingQCT.ValidWindows(strcmp(RecordingQCT.TaskType, taskTypesToPlot{tI}) & strcmp(RecordingQCT.Condition, conds{c}));
            if isempty(vals); continue; end
            jitterX = c + (rand(size(vals))-0.5)*0.3;
            scatter(jitterX, vals, 25, 'filled');
            errorbar(c, mean(vals), std(vals), 'k', 'LineWidth', 1.5, 'CapSize', 12);
        end
        hold off;
        xticks(1:numel(conds)); xticklabels(conds); xlim([0.5 numel(conds)+0.5]);
        ylabel('Valid windows per recording');
        title(sprintf('%s', taskTypesToPlot{tI}), 'Interpreter', 'none');
    end
    sgtitle('Diagnostic 3: Valid-window counts by condition, per task type (descriptive, not a test)');
    saveas(fig3, fullfile(SUBJ_DIR, 'diagnostic3_valid_windows_by_condition.png'));
    close(fig3);
catch ME3
    warning('Could not generate Diagnostic 3 figure: %s', ME3.message);
end

try
    fig4 = figure('Visible', 'off', 'Position', [100 100 1200 400]);
    conds = {'Baseline', 'E', 'N', 'R'};
    taskTypesToPlot = {'full', 'nback', 'sart'};
    for tI = 1:numel(taskTypesToPlot)
        subplot(1, numel(taskTypesToPlot), tI);
        hold on;
        for c = 1:numel(conds)
            vals = GlobalBandPowerT.AbsolutePower(strcmp(GlobalBandPowerT.Band, 'Alpha') & strcmp(GlobalBandPowerT.TaskType, taskTypesToPlot{tI}) & strcmp(GlobalBandPowerT.Condition, conds{c}));
            if isempty(vals); continue; end
            jitterX = c + (rand(size(vals))-0.5)*0.3;
            scatter(jitterX, vals, 25, 'filled');
            errorbar(c, mean(vals), std(vals)/sqrt(numel(vals)), 'k', 'LineWidth', 1.5, 'CapSize', 12);
        end
        hold off;
        xticks(1:numel(conds)); xticklabels(conds); xlim([0.5 numel(conds)+0.5]);
        ylabel('Whole-scalp Alpha AbsolutePower (uV^2)');
        title(sprintf('%s', taskTypesToPlot{tI}), 'Interpreter', 'none');
    end
    sgtitle({'Diagnostic 4: Alpha power by condition, per task type -- DESCRIPTIVE ONLY', 'error bars = mean +/- SEM; no statistical test performed'});
    saveas(fig4, fullfile(SUBJ_DIR, 'diagnostic4_alpha_power_by_condition.png'));
    close(fig4);
catch ME4
    warning('Could not generate Diagnostic 4 figure: %s', ME4.message);
end

% --- Diagnostic 5: extreme participant-level values -- SEPARATELY PER
%     TASK TYPE (see header item 2), pooling conditions within task type,
%     flag only, never exclude ---
extremeRows = {};
for t = 1:numel(taskTypesObserved)
    tt = taskTypesObserved{t};
    for b = 1:nBands
        band = bandNames{b};
        maskThis = strcmp(GlobalBandPowerT.Band, band) & strcmp(GlobalBandPowerT.TaskType, tt);
        vals = GlobalBandPowerT.AbsolutePower(maskThis);
        subjs = GlobalBandPowerT.Subject(maskThis);
        blocks = GlobalBandPowerT.BlockNumber(maskThis);
        conds = GlobalBandPowerT.Condition(maskThis);
        if isempty(vals); continue; end
        med = median(vals);
        madVal = median(abs(vals - med));
        robustSD = 1.4826 * madVal;
        if robustSD == 0
            continue;
        end
        z = (vals - med) / robustSD;
        flagged = find(abs(z) > 3);
        for f = flagged(:)'
            extremeRows(end+1, :) = {subjs{f}, tt, blocks(f), conds{f}, band, vals(f), z(f)}; %#ok<AGROW>
        end
    end
end
if isempty(extremeRows)
    ExtremeT = cell2table(cell(0,7), 'VariableNames', {'Subject','TaskType','BlockNumber','Condition','Band','AbsolutePower','RobustZ'});
else
    ExtremeT = cell2table(extremeRows, 'VariableNames', {'Subject','TaskType','BlockNumber','Condition','Band','AbsolutePower','RobustZ'});
end
ExtremeT = normalize_text_columns(ExtremeT);
writetable(ExtremeT, fullfile(SUBJ_DIR, 'diagnostic5_extreme_value_flags.csv'));

%% ------------------------------------------------------------------
%  8. FINAL PHASE 3 REPORT
% ------------------------------------------------------------------

report_lines = {};
report_lines{end+1} = '========================================================';
report_lines{end+1} = 'PHASE 3 -- COGNITIVE-TASK SUBJECT-LEVEL AGGREGATION REPORT';
report_lines{end+1} = '========================================================';
report_lines{end+1} = '';
report_lines{end+1} = sprintf('Input recordings          : %d', nRec);
report_lines{end+1} = sprintf('Subject-level recordings  : %d', nSubjLevel);
report_lines{end+1} = '';
report_lines{end+1} = 'Recording counts by TaskType (reported, not gated against a hard-coded literal -- see header item 3):';
for t = 1:numel(taskTypesObserved)
    report_lines{end+1} = sprintf('  %s : %d', taskTypesObserved{t}, taskCountsHere.(taskTypesObserved{t})); %#ok<AGROW>
end
report_lines{end+1} = '';
report_lines{end+1} = 'Condition counts by TaskType (includes Baseline):';
for t = 1:numel(taskTypesObserved)
    tt = taskTypesObserved{t};
    subT = RecordingQCT(strcmp(RecordingQCT.TaskType, tt), :);
    report_lines{end+1} = sprintf('  %-6s: Baseline=%d E=%d N=%d R=%d', tt, ...
        sum(strcmp(subT.Condition,'Baseline')), sum(strcmp(subT.Condition,'E')), ...
        sum(strcmp(subT.Condition,'N')), sum(strcmp(subT.Condition,'R'))); %#ok<AGROW>
end
report_lines{end+1} = '';
report_lines{end+1} = sprintf('Aggregation method: %s (window-level); LogAbsolutePower = ln(recording-level AbsolutePower)', ...
    PARAMS3.WindowAggregationStatistic);
report_lines{end+1} = '';
report_lines{end+1} = sprintf('Absolute power transformation: %s -- %s', ...
    iif(PARAMS3.AbsolutePowerLogTransform, 'ln (natural log) applied', 'none applied'), PARAMS3.LogTransformConvention);
report_lines{end+1} = '';
report_lines{end+1} = 'Valid-window summary (per recording, all tasks pooled):';
report_lines{end+1} = sprintf('  Median  : %g', median(RecordingQCT.ValidWindows));
report_lines{end+1} = sprintf('  Minimum : %g', min(RecordingQCT.ValidWindows));
report_lines{end+1} = sprintf('  Maximum : %g', max(RecordingQCT.ValidWindows));
report_lines{end+1} = '';
report_lines{end+1} = sprintf('Participant-level QC problems (integrity issues found): %d', numel(integrityIssues));
for k = 1:numel(integrityIssues)
    report_lines{end+1} = sprintf('  - %s', integrityIssues{k}); %#ok<AGROW>
end
report_lines{end+1} = '';
report_lines{end+1} = sprintf('NaN/Inf: AbsolutePower=%d, RelativePower=%d, LogAbsolutePower=%d (Log NaNs occur only when AbsolutePower<=0, individually flagged above)', ...
    nNaNInf_Abs, nNaNInf_Rel, nNaNInf_Log);
report_lines{end+1} = '';
report_lines{end+1} = sprintf('Extreme value flags (Diagnostic 5, per task type, flagged not removed): %d', height(ExtremeT));
report_lines{end+1} = '';
report_lines{end+1} = 'Output files:';
report_lines{end+1} = sprintf('  %s', fullfile(SUBJ_DIR, 'phase3_subject_level.mat'));
report_lines{end+1} = sprintf('  %s', fullfile(SUBJ_DIR, 'phase3_subject_level_bandpower.mat'));
report_lines{end+1} = sprintf('  %s', fullfile(SUBJ_DIR, 'phase3_recording_summary.xlsx'));
report_lines{end+1} = sprintf('  %s', fullfile(SUBJ_DIR, 'phase3_parameters.mat'));
report_lines{end+1} = sprintf('  %s', fullfile(SUBJ_DIR, 'diagnostic1_subject_bandpower_distributions.png'));
report_lines{end+1} = sprintf('  %s', fullfile(SUBJ_DIR, 'diagnostic2_mean_vs_median_window_aggregation.png'));
report_lines{end+1} = sprintf('  %s', fullfile(SUBJ_DIR, 'diagnostic3_valid_windows_by_condition.png'));
report_lines{end+1} = sprintf('  %s', fullfile(SUBJ_DIR, 'diagnostic4_alpha_power_by_condition.png'));
report_lines{end+1} = sprintf('  %s', fullfile(SUBJ_DIR, 'diagnostic5_extreme_value_flags.csv'));
report_lines{end+1} = sprintf('  %s', fullfile(REPORT_DIR, 'phase3_aggregation_report_cognitive.txt'));
report_lines{end+1} = '========================================================';

fullReport = strjoin(report_lines, '\n');
fprintf('%s\n', fullReport);

fid = fopen(fullfile(REPORT_DIR, 'phase3_aggregation_report_cognitive.txt'), 'w');
fprintf(fid, '%s\n', fullReport);
fclose(fid);

fprintf('\nPhase 3 complete. Do not proceed to X05 until this output has been reviewed.\n');


%% ========================================================================
%  LOCAL FUNCTIONS
%  (must follow ALL top-level executable statements per MATLAB's rule for
%  local functions in a script file)
% ========================================================================

function outArr = concat_struct_array(inArr, entry)
    if isempty(inArr)
        outArr = entry;
    else
        outArr = inArr;
        outArr(end+1) = entry; %#ok<AGROW>
    end
end

function T = normalize_text_columns(T)
    varNames = T.Properties.VariableNames;
    for v = 1:numel(varNames)
        col = T.(varNames{v});
        if ischar(col)
            T.(varNames{v}) = cellstr(col);
        end
    end
end

function v = local_iqr(X)
    v = local_prctile(X, 75) - local_prctile(X, 25);
end

function p = local_prctile(X, pct)
    [n, m] = size(X);
    p = nan(1, m);
    if n == 0
        return;
    end
    Xs = sort(X, 1);
    if n == 1
        p = Xs(1, :);
        return;
    end
    posPct = (pct/100) * n + 0.5;
    posPct = min(max(posPct, 1), n);
    lo = floor(posPct);
    hi = ceil(posPct);
    frac = posPct - lo;
    p = Xs(lo, :) + frac .* (Xs(hi, :) - Xs(lo, :));
end

function s = local_skewness(X)
    n = size(X, 1);
    mu = mean(X, 1);
    sigma = std(X, 1, 1);
    sigma(sigma == 0) = NaN;
    s = mean((X - mu).^3, 1) ./ (sigma.^3);
end

function out = iif(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end
