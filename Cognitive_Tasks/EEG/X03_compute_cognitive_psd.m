%% ========================================================================
% X03_compute_cognitive_psd.m
%
% Purpose: convert each already-preprocessed, continuous cognitive-task
% EEG recording (Baseline and E/N/R blocks, for each of Full/N-back/SART)
% into window-level and recording-level spectral measures (frequency-
% resolved PSD, absolute band power, relative band power), fully
% preserving the Subject x Task x Condition x BlockNumber x Window x
% Channel x Frequency structure, so X04 can aggregate correctly,
% separately per Task. Uses the same PSD methodology as the intervention
% pipeline's own PSD script, applied here to the cognitive-task metadata.
%
% Method:
%   - Analysis window length : 4 s  (2048 samples @ 512 Hz)
%   - Window overlap         : 50 % (2 s step)
%   - Taper                  : Hamming, applied to the full window
%   - Spectral estimator     : MATLAB pwelch(), one segment per window
%   - Retained frequency ranges : 1-45 Hz AND 55-100 Hz; 45-55 Hz notch
%                               guard band excluded entirely
%   - Frequency resolution   : Fs / NFFT = 512 / 2048 = 0.25 Hz
%   - Bands (Hz, both edges inclusive): Delta[1,4] Theta[4,8] Alpha[8,13]
%     Beta[13,30] Gamma_low_exploratory[30,45] Gamma_high1_exploratory
%     [55,80] Gamma_high2_exploratory[80,100]
%   - Relative-power denominator : total power integrated over 1-45 Hz only
%   - Window QC   : objective, condition-blind, fixed a priori criteria
%                   only (NaN/Inf; |amplitude| > 150 uV; flat channel
%                   std < 0.01 uV)
%   - Duration handling : uses ALL available valid windows per recording;
%                   never truncates to a common/shortest duration, and
%                   never rejects a recording merely because its duration
%                   differs from another recording or task type.
%
% This script does not filter/re-reference/run ICA/remove ICs/
% interpolate/resample/remove channels/modify any original .set file,
% perform inferential statistics, perform E vs N vs R (or Task vs Task)
% comparisons, or decide ROI definitions (deferred to X07). Individual
% windows are never treated as independent observations for any
% comparison (Subject, within Task, is the inferential unit); Full,
% N-back and SART are never compared to one another or modeled as a
% Task x Condition interaction.
%
% Requirements: EEGLAB on the MATLAB path, initialized enough for
% pop_loadset; Signal Processing Toolbox (pwelch, hamming); base MATLAB
% table/statistics functions; MATLAB R2016b or later (local functions in
% a script file). The X02 output file cognitive_metadata.mat must already
% exist.
%
% Inputs: DATA_DIR (preprocessed .set/.fdt files), X02's metadata output.
%
% Outputs:
%   <ANALYSIS_ROOT>/03_psd/
%       phase2_parameters.mat
%       cognitive_psd.mat
%       phase2_processing_log.csv / .xlsx / .mat
%       diagnostic1_example_psd_curves.png
%       diagnostic2_valid_window_count_hist.png
%       diagnostic4_band_power_summary.csv
%   <ANALYSIS_ROOT>/09_reports/
%       phase2_psd_report_cognitive.txt
%
% RUN_MODE: not applicable; this script runs once.
% ========================================================================

clear; clc;

%% ------------------------------------------------------------------
%  0. USER-EDITABLE SETTINGS
% ------------------------------------------------------------------

% >>> Must match the ANALYSIS_ROOT produced by X02_build_cognitive_metadata.m <<<
config;
DATA_ROOT_FOR_SIBLING = DATA_DIR;
[DATA_ROOT_PARENT, ~] = fileparts(DATA_ROOT_FOR_SIBLING);
ANALYSIS_ROOT = OUT_DIR;

META_DIR   = fullfile(ANALYSIS_ROOT, '02_metadata');
PSD_DIR    = fullfile(ANALYSIS_ROOT, '03_psd');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');
if ~exist(PSD_DIR, 'dir');    mkdir(PSD_DIR);    end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

META_MAT_FILE = fullfile(META_DIR, 'cognitive_metadata.mat');

fprintf('========================================================\n');
fprintf('PHASE 2 -- COGNITIVE-TASK CONTINUOUS EEG -> PSD\n');
fprintf('========================================================\n');
fprintf('ANALYSIS_ROOT : %s\n', ANALYSIS_ROOT);
fprintf('OUTPUT        : %s\n', PSD_DIR);
fprintf('--------------------------------------------------------\n\n');

if exist('pop_loadset', 'file') ~= 2
    error(['EEGLAB (pop_loadset) not found on the MATLAB path. ' ...
           'Start EEGLAB first, e.g.:  eeglab(''nogui'');  or  eeglab;']);
end
if exist('pwelch', 'file') ~= 2
    error(['pwelch() not found. Phase 2 requires the Signal Processing ' ...
           'Toolbox for Welch PSD estimation.']);
end
if ~exist(META_MAT_FILE, 'file')
    error(['Phase 1 metadata file not found: %s\nRun ' ...
           'X02_build_cognitive_metadata.m first.'], META_MAT_FILE);
end

%% ------------------------------------------------------------------
%  1. SPECTRAL-ANALYSIS PARAMETERS (saved verbatim for reproducibility;
%     IDENTICAL values to 03_compute_intervention_psd.m -- see that
%     script for full rationale of every choice)
% ------------------------------------------------------------------

PARAMS = struct();
PARAMS.ScriptName          = 'X03_compute_cognitive_psd.m';
PARAMS.GeneratedOn         = datestr(now);

PARAMS.ExpectedSamplingRateHz = 512;
PARAMS.ExpectedNumChannels    = 62;

PARAMS.WindowLengthSec        = 4;
PARAMS.WindowOverlapFraction  = 0.5;
PARAMS.WindowLengthSamples    = round(PARAMS.WindowLengthSec * PARAMS.ExpectedSamplingRateHz);
PARAMS.WindowStepSamples      = round(PARAMS.WindowLengthSamples * (1 - PARAMS.WindowOverlapFraction));

PARAMS.TaperType              = 'hamming';
PARAMS.NFFT                   = PARAMS.WindowLengthSamples;
PARAMS.FrequencyResolutionHz  = PARAMS.ExpectedSamplingRateHz / PARAMS.NFFT;

PARAMS.PrimaryAnalysisFreqRangeHz = [1 45];
PARAMS.HighGammaFreqRangeHz       = [55 100];
PARAMS.NotchGuardBandExcludedHz   = [45 55];
PARAMS.RetainedFreqRangesHz = {PARAMS.PrimaryAnalysisFreqRangeHz, PARAMS.HighGammaFreqRangeHz};

PARAMS.BandDefinitionsHz = struct( ...
    'Delta',                  [ 1   4], ...
    'Theta',                  [ 4   8], ...
    'Alpha',                  [ 8  13], ...
    'Beta',                   [13  30], ...
    'Beta_low_exploratory',   [13  20], ...
    'Beta_high_exploratory',  [20  30], ...
    'Gamma_low_exploratory',  [30  45], ...
    'Gamma_high1_exploratory',[55  80], ...
    'Gamma_high2_exploratory',[80 100]);
% NOTE: Beta_low_exploratory and Beta_high_exploratory are additional
% sub-divisions of the existing, unchanged locked Beta[13,30] band, added
% exactly like the pre-existing exploratory
% Gamma bands -- own separate BH-FDR family each (see X05/X08/X09/X11),
% never pooled with the locked Primary/Secondary tier bands, and Beta
% [13,30] itself is left completely untouched.
PARAMS.PrimaryBands    = {'Delta', 'Theta', 'Alpha', 'Beta'};
PARAMS.ExploratoryBands = {'Beta_low_exploratory', 'Beta_high_exploratory', ...
    'Gamma_low_exploratory', 'Gamma_high1_exploratory', 'Gamma_high2_exploratory'};
PARAMS.BandEdgeConvention = 'both edges inclusive at every boundary; adjacent bands share their boundary frequency sample (exact under trapz, not a double-count)';

PARAMS.RelativePowerDenominatorRangeHz = PARAMS.PrimaryAnalysisFreqRangeHz;

PARAMS.WindowQC = struct( ...
    'RejectIfNaNInf',            true, ...
    'AmplitudeThresholdUV',      150, ...
    'FlatChannelStdThresholdUV', 0.01);

PARAMS.AssumedDataUnits = 'microvolts (uV), per standard EEGLAB convention; not independently re-verified against amplifier calibration records.';

PARAMS.DurationHandlingPolicy = [ ...
    'Use ALL available analysis windows from each recording''s full ' ...
    'duration; do NOT truncate recordings to a common/shortest duration, ' ...
    'and do NOT reject a recording merely because its duration differs ' ...
    'from the intervention recordings or from another task type''s ' ...
    'recordings. Unequal per-recording window counts are NOT allowed to ' ...
    'bias later comparisons because X04 aggregates every subject x task ' ...
    'x condition to exactly one value before any group comparison.'];

PARAMS.SaveFullWindowLevelPSD = false;

fprintf('Window length      : %g s (%d samples)\n', PARAMS.WindowLengthSec, PARAMS.WindowLengthSamples);
fprintf('Window overlap     : %g %% (step = %d samples)\n', PARAMS.WindowOverlapFraction*100, PARAMS.WindowStepSamples);
fprintf('NFFT               : %d\n', PARAMS.NFFT);
fprintf('Frequency resolution: %.4f Hz\n', PARAMS.FrequencyResolutionHz);
fprintf('Frequency ranges retained: %g-%g Hz and %g-%g Hz (excluding %g-%g Hz notch guard band)\n\n', ...
    PARAMS.PrimaryAnalysisFreqRangeHz(1), PARAMS.PrimaryAnalysisFreqRangeHz(2), ...
    PARAMS.HighGammaFreqRangeHz(1), PARAMS.HighGammaFreqRangeHz(2), ...
    PARAMS.NotchGuardBandExcludedHz(1), PARAMS.NotchGuardBandExcludedHz(2));

%% ------------------------------------------------------------------
%  2. LOAD PHASE 1 METADATA (source of truth; not reparsed independently)
% ------------------------------------------------------------------

loadedMeta = load(META_MAT_FILE, 'T', 'ref_labels');
T = loadedMeta.T;
T = T(strcmp(T.LoadStatus, 'OK'), :);   % only successfully-loaded recordings proceed to PSD
ref_labels = loadedMeta.ref_labels(:);
nRefChan = numel(ref_labels);

nRec = height(T);
fprintf('Loaded Phase 1 metadata: %d recordings with LoadStatus=OK, %d reference channels.\n\n', nRec, nRefChan);

%% ------------------------------------------------------------------
%  3. SHARED, PRECOMPUTED QUANTITIES (identical for every recording,
%     since Fs is validated == 512 Hz for all recordings in Phase 1)
% ------------------------------------------------------------------

winSamples  = PARAMS.WindowLengthSamples;
stepSamples = PARAMS.WindowStepSamples;
NFFT        = PARAMS.NFFT;
taperWin    = hamming(winSamples);

fullFreqVec = (0:floor(NFFT/2))' * (PARAMS.ExpectedSamplingRateHz / NFFT);

freqMask = false(size(fullFreqVec));
for rr = 1:numel(PARAMS.RetainedFreqRangesHz)
    rng = PARAMS.RetainedFreqRangesHz{rr};
    freqMask = freqMask | (fullFreqVec >= rng(1) & fullFreqVec <= rng(2));
end
freqVecSel  = fullFreqVec(freqMask);
nFreqSel    = numel(freqVecSel);

bandNames = fieldnames(PARAMS.BandDefinitionsHz);
nBands = numel(bandNames);

fprintf('Selected frequency bins (union of retained ranges, excluding notch guard band): %d\n\n', nFreqSel);

%% ------------------------------------------------------------------
%  4. MAIN LOOP: LOAD EACH RECORDING, WINDOW, COMPUTE PSD / BAND POWER
% ------------------------------------------------------------------

psd_data = struct([]);
kSuccess = 0;

logCols = {'Subject', 'TaskType', 'BlockNumber', 'Condition', 'DurationSec', ...
           'NumCandidateWindows', 'NumValidWindows', 'PercentValid', ...
           'AnalysisDurationUsedSec', 'RejectedNaNInf', 'RejectedAmplitude', ...
           'RejectedFlat', 'ProcessingStatus', 'ErrorMessage'};
log_rows = cell(0, numel(logCols));

for i = 1:nRec

    subjStr   = getstr(T.Subject{i});
    taskStr   = getstr(T.TaskType{i});
    blockNum  = T.BlockNumber(i);
    condStr   = getstr(T.Condition{i});
    fpath     = getstr(T.FullPath{i});
    durPhase1 = T.DurationSec(i);

    fprintf('[%d/%d] Subject %s, Task %s, Block %g, Condition %s ... ', i, nRec, subjStr, taskStr, blockNum, condStr);

    status = 'OK';
    errMsg = '';
    nCandidate = 0; nValid = 0; pctValid = NaN; analysisDurUsed = NaN;
    nRejNaNInf = 0; nRejAmp = 0; nRejFlat = 0;

    try
        [fdir, fname, fext] = fileparts(fpath);
        EEG = pop_loadset('filename', [fname fext], 'filepath', fdir);

        Fs = EEG.srate;
        if Fs ~= PARAMS.ExpectedSamplingRateHz
            error('Sampling rate mismatch: expected %d Hz, found %g Hz.', ...
                PARAMS.ExpectedSamplingRateHz, Fs);
        end

        chanLabels = {EEG.chanlocs.labels};
        if numel(chanLabels) ~= nRefChan
            error('Channel count mismatch: expected %d, found %d.', nRefChan, numel(chanLabels));
        end

        reorderIdx = match_channel_order(chanLabels, ref_labels);

        data = double(EEG.data(reorderIdx, :));
        nSamplesRec = size(data, 2);
        nChan = size(data, 1);

        if nSamplesRec < winSamples
            error('Recording shorter than one analysis window (%d samples < %d).', ...
                nSamplesRec, winSamples);
        end

        starts = 0:stepSamples:(nSamplesRec - winSamples);
        nCandidate = numel(starts);

        windowValid       = false(nCandidate, 1);
        windowRejectReason = cell(nCandidate, 1);
        windowStartSec = (starts(:)) / Fs;
        windowEndSec   = (starts(:) + winSamples) / Fs;

        absBandPower = struct();
        relBandPower = struct();
        for b = 1:nBands
            absBandPower.(bandNames{b}) = nan(nCandidate, nChan);
            relBandPower.(bandNames{b}) = nan(nCandidate, nChan);
        end
        totalPower = nan(nCandidate, nChan);

        sumPSD = zeros(nFreqSel, nChan);

        for w = 1:nCandidate
            idx = (starts(w) + 1):(starts(w) + winSamples);
            seg = data(:, idx);
            segT = seg.';

            [isBad, reasonStr, isNaNInf, isAmp, isFlat] = check_window_qc( ...
                segT, PARAMS.WindowQC.AmplitudeThresholdUV, PARAMS.WindowQC.FlatChannelStdThresholdUV);

            if isNaNInf; nRejNaNInf = nRejNaNInf + 1; end
            if isAmp;    nRejAmp    = nRejAmp    + 1; end
            if isFlat;   nRejFlat   = nRejFlat   + 1; end

            windowRejectReason{w} = reasonStr;

            if isBad
                windowValid(w) = false;
                continue;
            end
            windowValid(w) = true;

            [Pxx, ~] = pwelch(segT, taperWin, 0, NFFT, Fs);
            PxxSel = Pxx(freqMask, :);

            sumPSD = sumPSD + PxxSel;

            totalPower(w, :) = integrate_band_power(freqVecSel, PxxSel, PARAMS.RelativePowerDenominatorRangeHz);

            for b = 1:nBands
                bandRange = PARAMS.BandDefinitionsHz.(bandNames{b});
                bp = integrate_band_power(freqVecSel, PxxSel, bandRange);
                absBandPower.(bandNames{b})(w, :) = bp;
                relBandPower.(bandNames{b})(w, :) = bp ./ totalPower(w, :);
            end
        end

        nValid = sum(windowValid);
        pctValid = 100 * nValid / max(nCandidate, 1);
        analysisDurUsed = nValid * PARAMS.WindowLengthSec;

        if nValid > 0
            meanPSD = (sumPSD / nValid).';
        else
            meanPSD = nan(nChan, nFreqSel);
            warning('Recording sub%s task=%s block=%g (%s): 0 valid windows out of %d candidates.', ...
                subjStr, taskStr, blockNum, fname, nCandidate);
        end

        entry = struct();
        entry.Subject = subjStr;
        entry.TaskType = taskStr;
        entry.BlockNumber = blockNum;
        entry.Condition = condStr;
        entry.ConditionLetterRaw = getstr(T.ConditionLetterRaw{i});
        entry.FullPath = fpath;
        entry.DurationSec = durPhase1;
        entry.ChannelLabels = ref_labels;
        entry.FrequenciesHz = freqVecSel;
        entry.WindowStartSec = windowStartSec;
        entry.WindowEndSec = windowEndSec;
        entry.WindowValid = windowValid;
        entry.WindowRejectReason = windowRejectReason;
        entry.AbsoluteBandPower = absBandPower;
        entry.RelativeBandPower = relBandPower;
        entry.TotalPower_1to45Hz = totalPower;
        entry.MeanPSD = meanPSD;
        entry.NumCandidateWindows = nCandidate;
        entry.NumValidWindows = nValid;

        if PARAMS.SaveFullWindowLevelPSD
            fullPSDArray = nan(nCandidate, nChan, nFreqSel);
            for w = 1:nCandidate
                if ~windowValid(w); continue; end
                idx = (starts(w) + 1):(starts(w) + winSamples);
                segT = data(:, idx).';
                [Pxx, ~] = pwelch(segT, taperWin, 0, NFFT, Fs);
                fullPSDArray(w, :, :) = Pxx(freqMask, :).';
            end
            entry.PSD_WindowLevel = fullPSDArray;
        end

        psd_data = concat_struct_array(psd_data, entry);
        kSuccess = kSuccess + 1;

        fprintf('OK (%d/%d valid windows, %.1f%%)\n', nValid, nCandidate, pctValid);

    catch ME
        status = 'FAILED';
        errMsg = ME.message;
        fprintf('FAILED (%s)\n', errMsg);
    end

    log_rows(end+1, :) = {subjStr, taskStr, blockNum, condStr, durPhase1, nCandidate, nValid, ...
        pctValid, analysisDurUsed, nRejNaNInf, nRejAmp, nRejFlat, status, errMsg}; %#ok<AGROW>
end

kSuccess = numel(psd_data);
fprintf('\nMain loop complete: %d/%d recordings processed successfully.\n\n', kSuccess, nRec);

%% ------------------------------------------------------------------
%  5. SAVE OUTPUTS
% ------------------------------------------------------------------

save(fullfile(PSD_DIR, 'phase2_parameters.mat'), 'PARAMS');

FrequencyVectorHz = freqVecSel;
ChannelLabelsReference = ref_labels;
save(fullfile(PSD_DIR, 'cognitive_psd.mat'), ...
    'psd_data', 'FrequencyVectorHz', 'ChannelLabelsReference', 'PARAMS', '-v7.3');

ProcessingLogT = cell2table(log_rows, 'VariableNames', logCols);
ProcessingLogT = normalize_text_columns(ProcessingLogT);
writetable(ProcessingLogT, fullfile(PSD_DIR, 'phase2_processing_log.csv'));
writetable(ProcessingLogT, fullfile(PSD_DIR, 'phase2_processing_log.xlsx'));
save(fullfile(PSD_DIR, 'phase2_processing_log.mat'), 'ProcessingLogT');

%% ------------------------------------------------------------------
%  6. DIAGNOSTICS
% ------------------------------------------------------------------

nConsistencyFailures = sum(strcmp(ProcessingLogT.ProcessingStatus, 'FAILED') & ...
    contains(ProcessingLogT.ErrorMessage, {'Sampling rate mismatch', 'Channel count mismatch', 'Channel labels'}));

validCounts = arrayfun(@(x) x.NumValidWindows, psd_data);
try
    fig2 = figure('Visible', 'off');
    histogram(validCounts);
    xlabel('Number of valid windows per recording');
    ylabel('Number of recordings');
    title('Diagnostic 2: Distribution of valid-window counts (all tasks pooled)');
    saveas(fig2, fullfile(PSD_DIR, 'diagnostic2_valid_window_count_hist.png'));
    close(fig2);
catch ME2
    warning('Could not generate Diagnostic 2 figure: %s', ME2.message);
end

% --- Diagnostic 1: example PSD curves, one representative recording per
%     Task x Condition (now includes Baseline; a nTaskTypes x
%     nConditions subplot grid, extended from the intervention pipeline's
%     single row of 3 -- see header item 3 for why) ---
try
    exampleChans = {'Fz', 'Cz', 'Oz'};
    exampleChanIdx = find(ismember(ref_labels, exampleChans));
    if isempty(exampleChanIdx)
        exampleChanIdx = 1:min(3, numel(ref_labels));
    end
    taskTypesToPlot = {'full', 'nback', 'sart'};
    condsToPlot = {'Baseline', 'E', 'N', 'R'};
    fig1 = figure('Visible', 'off', 'Position', [50 50 1600 900]);
    plotIdx = 0;
    for tI = 1:numel(taskTypesToPlot)
        for c = 1:numel(condsToPlot)
            plotIdx = plotIdx + 1;
            subplot(numel(taskTypesToPlot), numel(condsToPlot), plotIdx);
            taskMask = arrayfun(@(x) strcmp(x.TaskType, taskTypesToPlot{tI}) && strcmp(x.Condition, condsToPlot{c}), psd_data);
            idxList = find(taskMask);
            if isempty(idxList)
                title(sprintf('%s / %s (none)', taskTypesToPlot{tI}, condsToPlot{c}), 'Interpreter', 'none');
                continue;
            end
            exRec = psd_data(idxList(1));
            hold on;
            plotFreqBase = exRec.FrequenciesHz(:)';
            gapIdx = find(diff(plotFreqBase) > 2 * PARAMS.FrequencyResolutionHz);
            for ch = exampleChanIdx(:)'
                plotFreq_g = plotFreqBase;
                plotPSD = 10*log10(exRec.MeanPSD(ch, :));
                for g = sort(gapIdx(:)', 'descend')
                    plotFreq_g = [plotFreq_g(1:g), NaN, plotFreq_g(g+1:end)]; %#ok<AGROW>
                    plotPSD    = [plotPSD(1:g),    NaN, plotPSD(g+1:end)]; %#ok<AGROW>
                end
                plot(plotFreq_g, plotPSD, 'LineWidth', 1.0);
            end
            hold off;
            xlabel('Hz'); ylabel('dB');
            title(sprintf('%s / %s: Sub%s', taskTypesToPlot{tI}, condsToPlot{c}, exRec.Subject), 'Interpreter', 'none');
            if plotIdx == 1
                legend(ref_labels(exampleChanIdx), 'Location', 'northeast');
            end
            grid on;
        end
    end
    sgtitle('Diagnostic 1: Example mean PSD curves (representative recordings/channels, per Task x Condition)');
    saveas(fig1, fullfile(PSD_DIR, 'diagnostic1_example_psd_curves.png'));
    close(fig1);
catch ME1
    warning('Could not generate Diagnostic 1 figure: %s', ME1.message);
end

% --- Diagnostic 4: descriptive summary of absolute/relative band power,
%     broken out by TaskType x ConditionGroup (POOLED-across-tasks 'ALL'
%     rows are also included for a quick overall sanity check). Purely
%     descriptive, NOT an inferential comparison. ---
diag4_rows = {};
taskGroupsForSummary = {'ALL', 'full', 'nback', 'sart'};
condsForSummary = {'ALL', 'Baseline', 'E', 'N', 'R'};
for tg = 1:numel(taskGroupsForSummary)
    tgLabel = taskGroupsForSummary{tg};
    for c = 1:numel(condsForSummary)
        cLabel = condsForSummary{c};
        if strcmp(tgLabel, 'ALL')
            taskMaskIdx = true(1, numel(psd_data));
        else
            taskMaskIdx = arrayfun(@(x) strcmp(x.TaskType, tgLabel), psd_data);
        end
        if strcmp(cLabel, 'ALL')
            condMaskIdx = true(1, numel(psd_data));
        else
            condMaskIdx = arrayfun(@(x) strcmp(x.Condition, cLabel), psd_data);
        end
        recIdx = find(taskMaskIdx & condMaskIdx);
        for b = 1:nBands
            allAbs = []; allRel = [];
            for r = recIdx
                v = psd_data(r).WindowValid;
                absVals = psd_data(r).AbsoluteBandPower.(bandNames{b})(v, :);
                relVals = psd_data(r).RelativeBandPower.(bandNames{b})(v, :);
                allAbs = [allAbs; absVals(:)]; %#ok<AGROW>
                allRel = [allRel; relVals(:)]; %#ok<AGROW>
            end
            diag4_rows(end+1, :) = {tgLabel, cLabel, bandNames{b}, mean(allAbs, 'omitnan'), std(allAbs, 'omitnan'), ...
                mean(allRel, 'omitnan'), std(allRel, 'omitnan'), numel(allAbs)}; %#ok<AGROW>
        end
    end
end
Diagnostic4T = cell2table(diag4_rows, 'VariableNames', ...
    {'TaskTypeGroup', 'ConditionGroup', 'Band', 'MeanAbsolutePower_uV2', 'SDAbsolutePower_uV2', ...
     'MeanRelativePower', 'SDRelativePower', 'NumChannelWindowSamples'});
Diagnostic4T = normalize_text_columns(Diagnostic4T);
writetable(Diagnostic4T, fullfile(PSD_DIR, 'diagnostic4_band_power_summary.csv'));

%% ------------------------------------------------------------------
%  7. FINAL PHASE 2 REPORT
% ------------------------------------------------------------------

nFailed = sum(strcmp(ProcessingLogT.ProcessingStatus, 'FAILED'));
totalCandidate = sum(ProcessingLogT.NumCandidateWindows);
totalValid = sum(ProcessingLogT.NumValidWindows);
overallPctValid = 100 * totalValid / max(totalCandidate, 1);

condCounts = struct('Baseline', 0, 'E', 0, 'N', 0, 'R', 0, 'Unknown', 0);
taskCounts = struct('full', 0, 'nback', 0, 'sart', 0, 'Unknown', 0);
for r = 1:numel(psd_data)
    c = psd_data(r).Condition;
    if ~isempty(c) && isfield(condCounts, c)
        condCounts.(c) = condCounts.(c) + 1;
    else
        condCounts.Unknown = condCounts.Unknown + 1;
    end
    t = psd_data(r).TaskType;
    if ~isempty(t) && isfield(taskCounts, t)
        taskCounts.(t) = taskCounts.(t) + 1;
    else
        taskCounts.Unknown = taskCounts.Unknown + 1;
    end
end

report_lines = {};
report_lines{end+1} = '========================================================';
report_lines{end+1} = 'PHASE 2 -- COGNITIVE-TASK CONTINUOUS EEG -> PSD REPORT';
report_lines{end+1} = '========================================================';
report_lines{end+1} = '';
report_lines{end+1} = sprintf('Recordings processed : %d', nRec);
report_lines{end+1} = sprintf('Recordings successful: %d', kSuccess);
report_lines{end+1} = sprintf('Recordings failed     : %d', nFailed);
report_lines{end+1} = '';
report_lines{end+1} = sprintf('Window length      : %g s (%d samples)', PARAMS.WindowLengthSec, PARAMS.WindowLengthSamples);
report_lines{end+1} = sprintf('Overlap            : %g %% (step = %d samples = %.2f s)', ...
    PARAMS.WindowOverlapFraction*100, PARAMS.WindowStepSamples, PARAMS.WindowStepSamples/PARAMS.ExpectedSamplingRateHz);
report_lines{end+1} = sprintf('FFT length         : %d', PARAMS.NFFT);
report_lines{end+1} = sprintf('Frequency resolution: %.4f Hz', PARAMS.FrequencyResolutionHz);
report_lines{end+1} = '';
report_lines{end+1} = sprintf('Frequency ranges retained: %g-%g Hz and %g-%g Hz', ...
    PARAMS.PrimaryAnalysisFreqRangeHz(1), PARAMS.PrimaryAnalysisFreqRangeHz(2), ...
    PARAMS.HighGammaFreqRangeHz(1), PARAMS.HighGammaFreqRangeHz(2));
report_lines{end+1} = sprintf('Notch guard band excluded: %g-%g Hz (never analyzed/reported)', ...
    PARAMS.NotchGuardBandExcludedHz(1), PARAMS.NotchGuardBandExcludedHz(2));
report_lines{end+1} = sprintf('Relative-power denominator: %g-%g Hz only', ...
    PARAMS.RelativePowerDenominatorRangeHz(1), PARAMS.RelativePowerDenominatorRangeHz(2));
report_lines{end+1} = '';
report_lines{end+1} = 'Frequency bands:';
for b = 1:nBands
    br = PARAMS.BandDefinitionsHz.(bandNames{b});
    report_lines{end+1} = sprintf('  %-18s: %g-%g Hz', bandNames{b}, br(1), br(2)); %#ok<AGROW>
end
report_lines{end+1} = '';
report_lines{end+1} = sprintf('Total candidate windows : %d', totalCandidate);
report_lines{end+1} = sprintf('Total valid windows     : %d', totalValid);
report_lines{end+1} = sprintf('Overall valid percentage: %.2f %%', overallPctValid);
report_lines{end+1} = '';
report_lines{end+1} = 'Duration-handling strategy:';
report_lines{end+1} = ['  ' PARAMS.DurationHandlingPolicy];
report_lines{end+1} = '';
report_lines{end+1} = 'Window-QC strategy:';
report_lines{end+1} = sprintf('  Reject window if: NaN/Inf present; OR any |amplitude| > %g uV; OR any channel std < %g uV.', ...
    PARAMS.WindowQC.AmplitudeThresholdUV, PARAMS.WindowQC.FlatChannelStdThresholdUV);
report_lines{end+1} = sprintf('  Total windows rejected for NaN/Inf   : %d', sum(ProcessingLogT.RejectedNaNInf));
report_lines{end+1} = sprintf('  Total windows rejected for amplitude : %d', sum(ProcessingLogT.RejectedAmplitude));
report_lines{end+1} = sprintf('  Total windows rejected for flatness  : %d', sum(ProcessingLogT.RejectedFlat));
report_lines{end+1} = '';
report_lines{end+1} = 'Task type counts processed:';
report_lines{end+1} = sprintf('  full  : %d', taskCounts.full);
report_lines{end+1} = sprintf('  nback : %d', taskCounts.nback);
report_lines{end+1} = sprintf('  sart  : %d', taskCounts.sart);
if taskCounts.Unknown > 0
    report_lines{end+1} = sprintf('  Unknown/unmatched : %d  (unexpected)', taskCounts.Unknown); %#ok<AGROW>
end
report_lines{end+1} = '';
report_lines{end+1} = 'Condition counts processed (includes Baseline; Baseline is excluded from the';
report_lines{end+1} = 'inferential model only in X05/X06/X08/X09, not here):';
report_lines{end+1} = sprintf('  Baseline : %d', condCounts.Baseline);
report_lines{end+1} = sprintf('  E : %d', condCounts.E);
report_lines{end+1} = sprintf('  N : %d', condCounts.N);
report_lines{end+1} = sprintf('  R : %d', condCounts.R);
if condCounts.Unknown > 0
    report_lines{end+1} = sprintf('  Unknown/unmatched : %d  (unexpected -- Phase 1 should have excluded this)', condCounts.Unknown); %#ok<AGROW>
end
report_lines{end+1} = '';
report_lines{end+1} = sprintf('Diagnostic 5 (channel/Fs consistency failures): %d', nConsistencyFailures);
report_lines{end+1} = '';
report_lines{end+1} = 'Output files:';
report_lines{end+1} = sprintf('  %s', fullfile(PSD_DIR, 'phase2_parameters.mat'));
report_lines{end+1} = sprintf('  %s', fullfile(PSD_DIR, 'cognitive_psd.mat'));
report_lines{end+1} = sprintf('  %s', fullfile(PSD_DIR, 'phase2_processing_log.csv'));
report_lines{end+1} = sprintf('  %s', fullfile(PSD_DIR, 'phase2_processing_log.xlsx'));
report_lines{end+1} = sprintf('  %s', fullfile(PSD_DIR, 'diagnostic1_example_psd_curves.png'));
report_lines{end+1} = sprintf('  %s', fullfile(PSD_DIR, 'diagnostic2_valid_window_count_hist.png'));
report_lines{end+1} = sprintf('  %s', fullfile(PSD_DIR, 'diagnostic4_band_power_summary.csv'));
report_lines{end+1} = sprintf('  %s', fullfile(REPORT_DIR, 'phase2_psd_report_cognitive.txt'));
report_lines{end+1} = '========================================================';

fullReport = strjoin(report_lines, '\n');
fprintf('%s\n', fullReport);

fid = fopen(fullfile(REPORT_DIR, 'phase2_psd_report_cognitive.txt'), 'w');
fprintf(fid, '%s\n', fullReport);
fclose(fid);

fprintf('\nPhase 2 complete. Do not proceed to X04 until this output has been reviewed.\n');


%% ========================================================================
%  LOCAL FUNCTIONS
%  (must follow ALL top-level executable statements per MATLAB's rule for
%  local functions in a script file)
% ========================================================================

function s = getstr(x)
    if iscell(x)
        if isempty(x)
            s = '';
            return;
        end
        x = x{1};
    end
    if ischar(x)
        s = x;
    elseif isstring(x)
        s = char(x);
    elseif isnumeric(x) || islogical(x)
        s = num2str(x);
    else
        s = char(x);
    end
end

function idx = match_channel_order(chanLabels, refLabels)
    idx = nan(numel(refLabels), 1);
    for k = 1:numel(refLabels)
        hit = find(strcmpi(chanLabels, refLabels{k}), 1);
        if isempty(hit)
            error('Channel "%s" (from the Phase 1 reference channel order) not found in this recording.', refLabels{k});
        end
        idx(k) = hit;
    end
    if numel(unique(idx)) ~= numel(idx)
        error('Channel matching produced duplicate indices; channel labels are not uniquely resolvable.');
    end
    idx = round(idx);
end

function [isBad, reasonStr, isNaNInf, isAmp, isFlat] = check_window_qc(segT, ampThreshUV, flatStdThreshUV)
    isNaNInf = any(~isfinite(segT(:)));
    isAmp    = max(abs(segT(:))) > ampThreshUV;
    chanStd  = std(segT, 0, 1);
    isFlat   = any(chanStd < flatStdThreshUV);

    reasons = {};
    if isNaNInf; reasons{end+1} = 'NaNInf'; end %#ok<AGROW>
    if isAmp;    reasons{end+1} = 'Amplitude'; end %#ok<AGROW>
    if isFlat;   reasons{end+1} = 'Flat'; end %#ok<AGROW>

    isBad = isNaNInf || isAmp || isFlat;
    reasonStr = strjoin(reasons, '+');
end

function bp = integrate_band_power(freqVec, Pxx, bandRange)
    mask = freqVec >= bandRange(1) & freqVec <= bandRange(2);
    if ~any(mask)
        bp = nan(1, size(Pxx, 2));
        return;
    end
    bp = trapz(freqVec(mask), Pxx(mask, :), 1);
end

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
