%% ========================================================================
% X04_aggregate_subject_level.m
%
% Band power per recording, channel and band from the X03 spectra: absolute power,
% ln(absolute power) and relative power (band power / 1-45 Hz power, ratio of means),
% plus whole-scalp (global) values. Bands: delta 1-4, theta 4-8, alpha 8-13, beta 13-30,
% exploratory beta 13-20 and 20-30, exploratory gamma 30-45, 55-80 and 80-100 Hz.
% Inputs : 02_metadata, 03_psd.
% Outputs: 04_subject_level/phase3_subject_level_bandpower.mat (BandPowerT,
%          GlobalBandPowerT, PARAMS3), recording summaries, 09_reports/.
% ========================================================================

clear; clc;

%% ------------------------------------------------------------------

% >>> Must match the ANALYSIS_ROOT used in Phases 1 and 2 <<<
config;
ANALYSIS_ROOT = OUT_DIR;

PSD_DIR     = fullfile(ANALYSIS_ROOT, '03_psd');
SUBJ_DIR    = fullfile(ANALYSIS_ROOT, '04_subject_level');
REPORT_DIR  = fullfile(ANALYSIS_ROOT, '09_reports');
if ~exist(SUBJ_DIR, 'dir');   mkdir(SUBJ_DIR);   end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

PSD_MAT_FILE = fullfile(PSD_DIR, 'intervention_psd.mat');
LOG_MAT_FILE = fullfile(PSD_DIR, 'phase2_processing_log.mat');
META_MAT_FILE = fullfile(ANALYSIS_ROOT, '02_metadata', 'intervention_metadata.mat');

fprintf('========================================================\n');
fprintf('PHASE 3 -- SUBJECT-LEVEL AGGREGATION\n');
fprintf('========================================================\n');
fprintf('ANALYSIS_ROOT : %s\n', ANALYSIS_ROOT);
fprintf('OUTPUT        : %s\n', SUBJ_DIR);
fprintf('--------------------------------------------------------\n\n');

if ~exist(PSD_MAT_FILE, 'file')
    error(['Phase 2 PSD file not found: %s\nRun ' ...
           '03_compute_intervention_psd.m first.'], PSD_MAT_FILE);
end
if ~exist(LOG_MAT_FILE, 'file')
    error(['Phase 2 processing log not found: %s\nRun ' ...
           '03_compute_intervention_psd.m first.'], LOG_MAT_FILE);
end
if ~exist(META_MAT_FILE, 'file')
    error(['Phase 1 metadata not found: %s\nRun ' ...
           '02_build_intervention_metadata.m first.'], META_MAT_FILE);
end

%% ------------------------------------------------------------------
%  1. LOAD PHASE 2 OUTPUT (source of truth; PSD is never recomputed)
% ------------------------------------------------------------------

psdFile = load(PSD_MAT_FILE, 'psd_data', 'FrequencyVectorHz', 'ChannelLabelsReference', 'PARAMS');
psd_data   = psdFile.psd_data;
freqVecRef = psdFile.FrequencyVectorHz(:);
chanLabelsRef = psdFile.ChannelLabelsReference(:);
PARAMS2 = psdFile.PARAMS;   % Phase 2's saved parameters -- source of truth for bands/denominator

logFile = load(LOG_MAT_FILE, 'ProcessingLogT');
ProcessingLogT = logFile.ProcessingLogT;

if ~ismember('RejectedBoundary', ProcessingLogT.Properties.VariableNames)
    error(['The Phase 2 output in %s was written by the OLD 03 script (no boundary ' ...
           'handling). Rerun 03_compute_intervention_psd.m (v2) successfully first.'], PSD_DIR);
end

nRec = numel(psd_data);
bandNames = fieldnames(PARAMS2.BandDefinitionsHz);
nBands = numel(bandNames);
nChanRef = numel(chanLabelsRef);
nFreqRef = numel(freqVecRef);

fprintf('Loaded Phase 2 output: %d recordings, %d bands, %d channels, %d frequency bins.\n\n', ...
    nRec, nBands, nChanRef, nFreqRef);

%% ------------------------------------------------------------------
%  2. PHASE 3 AGGREGATION PARAMETERS (methodological choices, saved
%     verbatim for reproducibility; see accompanying explanation)
% ------------------------------------------------------------------

PARAMS3 = struct();
PARAMS3.ScriptName  = '04_aggregate_subject_level.m';
PARAMS3.GeneratedOn = datestr(now);
PARAMS3.SourcePhase2File = PSD_MAT_FILE;

PARAMS3.WindowAggregationStatistic = 'arithmetic_mean';

% --- Absolute-power log transformation ---
% LogAbsolutePower = natural log (ln) of the already-aggregated,
% recording-level AbsolutePower (i.e. ln(mean(window power)), NOT
% mean(ln(window power)) -- a simple, monotonic transform of the FINAL
% per-recording value, not a different window-aggregation method). This
% keeps LogAbsolutePower and AbsolutePower numerically consistent by
% construction (one is just ln() of the other) and keeps Structure A
% (frequency-resolved) and Structure B (band power) using the identical
% aggregation convention throughout.
% Rationale: absolute EEG power is classically log-normally distributed
% ACROSS SUBJECTS/RECORDINGS (individual differences in skull thickness,
% electrode impedance, cortical folding etc. act multiplicatively on
% overall amplitude) -- a well-established reason to log-transform before
% later parametric statistics, and a DIFFERENT concern from within-
% recording window-level skew (addressed empirically by Diagnostic 2/
% skewness below, not assumed).
% In the X03 diagnostic output the SD of absolute power is at least as
% large as the mean in every band, consistent with a right-skewed,
% approximately log-normal distribution.
PARAMS3.AbsolutePowerLogTransform = true;
PARAMS3.LogTransformConvention = 'natural log (ln), applied to the recording-level AbsolutePower AFTER window-averaging';
PARAMS3.RelativePowerLogTransform = false;
PARAMS3.RelativePowerDefinition = 'ratio of means: mean(window band power) / mean(window 1-45 Hz total power), per channel; whole-scalp = mean across channels';
PARAMS3.RelativePowerLogTransformRationale = ['Relative power is a bounded proportion, typically far less skewed ' ...
    'than absolute power. No log-transform is applied here; if Phase 4 diagnostics show it is still ' ...
    'problematically skewed, a logit transform can be added then rather than assumed now.'];

% --- Global scalp aggregation ---
% Arithmetic mean ACROSS the 62 channels, computed from the already
% window-aggregated per-channel values (not a separate re-aggregation
% from raw windows). Rationale: preprocessing already interpolated any
% bad channels back into a complete, cleaned 62-channel montage, so
% (unlike windows) there is no strong a priori reason to expect outlier
% channels requiring a robust summary. A median-across-channels variant is
% ALSO computed and saved alongside, as a secondary cross-check, not as
% primary.
PARAMS3.GlobalScalpAggregation = 'arithmetic_mean_across_channels';
PARAMS3.GlobalScalpSecondary   = 'median_across_channels (saved alongside, not primary)';

% --- Channel aggregation (band power / frequency-resolved) ---
PARAMS3.ChannelHandling = 'Channels are NEVER averaged away at the per-channel product level (Structures A and B). A separate whole-scalp summary (Structure C) is derived FROM the per-channel values, not a replacement for them.';

% --- Frequency bands & relative-power denominator: inherited verbatim ---
PARAMS3.BandDefinitionsHz = PARAMS2.BandDefinitionsHz;
PARAMS3.PrimaryBands = PARAMS2.PrimaryBands;
PARAMS3.ExploratoryBands = PARAMS2.ExploratoryBands;
PARAMS3.RelativePowerDenominatorRangeHz = PARAMS2.RelativePowerDenominatorRangeHz;
PARAMS3.NotchGuardBandExcludedHz = PARAMS2.NotchGuardBandExcludedHz;
PARAMS3.BandEdgeConvention = PARAMS2.BandEdgeConvention;

% --- Valid-window handling ---
PARAMS3.ValidWindowHandling = ['Only windows with WindowValid=true (per Phase 2 QC) contribute to any ' ...
    'aggregation. Each recording contributes exactly ONE row per channel per band (and one whole-scalp ' ...
    'row per band) regardless of how many valid windows it had -- unequal window counts across ' ...
    'recordings are therefore never a source of unequal weighting in this or any later phase.'];

PARAMS3.ExtremeValueRule = 'robust z = (x - median)/(1.4826*MAD) on LogAbsolutePower, computed per band WITHIN each condition; |z| > 3 flagged; flag only, never exclude';

save(fullfile(SUBJ_DIR, 'phase3_parameters.mat'), 'PARAMS3');

fprintf('Window-aggregation statistic : %s\n', PARAMS3.WindowAggregationStatistic);
fprintf('Absolute-power log transform : %d (%s)\n', PARAMS3.AbsolutePowerLogTransform, PARAMS3.LogTransformConvention);
fprintf('Global scalp aggregation     : %s\n\n', PARAMS3.GlobalScalpAggregation);

%% ------------------------------------------------------------------
%  3. MAIN LOOP: AGGREGATE EACH RECORDING
% ------------------------------------------------------------------

subject_psd = struct([]);   % Structure A: one entry per recording

bandPowerCols = {'Subject', 'InterventionNumber', 'Condition', 'Channel', 'Band', ...
    'AbsolutePower', 'LogAbsolutePower', 'RelativePower', 'RelativePower_MeanOfWindowRatios', ...
    'WindowMedian_Absolute', 'WindowSD_Absolute', 'WindowIQR_Absolute', 'WindowSkewness_Absolute', ...
    'WindowMedian_Relative', 'WindowSD_Relative', 'WindowIQR_Relative', 'WindowSkewness_Relative'};
bandPowerRows = cell(0, numel(bandPowerCols));

globalCols = {'Subject', 'InterventionNumber', 'Condition', 'Band', ...
    'AbsolutePower', 'LogAbsolutePower', 'RelativePower', 'RelativePower_MeanOfWindowRatios', ...
    'AbsolutePower_MedianAcrossChannels', 'RelativePower_MedianAcrossChannels'};
globalRows = cell(0, numel(globalCols));

qcCols = {'Subject', 'InterventionNumber', 'Condition', 'ValidWindows', 'RejectedWindows', ...
    'CandidateWindows', 'PercentValid', 'DurationSec', ...
    'RejectedNaNInf', 'RejectedAmplitude', 'RejectedFlat', 'RejectedBoundary', 'NumBoundaryEvents'};
qcRows = cell(0, numel(qcCols));

integrityIssues = {};   % free-text log of anything that fails a check below

for i = 1:nRec
    rec = psd_data(i);
    subjStr = rec.Subject;
    intNum  = rec.InterventionNumber;
    condStr = rec.Condition;

    fprintf('[%d/%d] Subject %s, Intervention %g, Condition %s ... ', i, nRec, subjStr, intNum, condStr);

    % ---- Consistency checks against the reference channel/frequency sets ----
    if numel(rec.ChannelLabels) ~= nChanRef || ~isequal(rec.ChannelLabels, chanLabelsRef)
        integrityIssues{end+1} = sprintf('Sub%s Int%g: channel labels differ from the reference set.', subjStr, intNum); %#ok<AGROW>
    end
    if numel(rec.FrequenciesHz) ~= nFreqRef || max(abs(rec.FrequenciesHz - freqVecRef)) > 1e-9
        integrityIssues{end+1} = sprintf('Sub%s Int%g: frequency vector differs from the reference vector.', subjStr, intNum); %#ok<AGROW>
    end

    valid = rec.WindowValid;
    nValidThis = sum(valid);

    % ---- Structure D: recording-level QC, pulled from Phase 2's own log.
    %      Populated FIRST and unconditionally, so that even a recording
    %      with 0 valid windows (which cannot be aggregated into
    %      Structures A/B/C below) still appears in the QC table rather
    %      than silently vanishing from the output entirely. ----
    logIdx = find(strcmp(ProcessingLogT.Subject, subjStr) & ProcessingLogT.InterventionNumber == intNum);
    if numel(logIdx) ~= 1
        integrityIssues{end+1} = sprintf('Sub%s Int%g: expected exactly 1 matching Phase 2 processing-log row, found %d.', ...
            subjStr, intNum, numel(logIdx)); %#ok<AGROW>
        candW = rec.NumCandidateWindows; validW = rec.NumValidWindows; pctV = NaN; durS = NaN;
        rNaN = NaN; rAmp = NaN; rFlat = NaN; rBnd = NaN; nBnd = NaN;
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
        if ismember('RejectedBoundary', ProcessingLogT.Properties.VariableNames)
            rBnd = lr.RejectedBoundary; nBnd = lr.NumBoundaryEvents;
        else
            rBnd = NaN; nBnd = NaN;
        end
    end
    qcRows(end+1, :) = {subjStr, intNum, condStr, validW, candW - validW, candW, pctV, durS, rNaN, rAmp, rFlat, rBnd, nBnd}; %#ok<AGROW>

    if nValidThis == 0
        integrityIssues{end+1} = sprintf('Sub%s Int%g: 0 valid windows -- cannot aggregate into Structures A/B/C; recorded in the QC table only.', subjStr, intNum); %#ok<AGROW>
        fprintf('SKIPPED for aggregation (0 valid windows; QC row still recorded)\n');
        continue;
    end

    % ---- Structure A: frequency-resolved aggregation (reuses Phase 2's
    %      own arithmetic-mean MeanPSD; relative PSD normalizes by the
    %      SAME TotalPower_1to45Hz field Phase 2 already defined, averaged
    %      the same way, so numerator and denominator share one
    %      aggregation convention throughout) ----
    aggAbsPSD = rec.MeanPSD;                                    % nChan x nFreq, already arithmetic-mean-of-valid-windows
    meanTotalPower = mean(rec.TotalPower_1to45Hz(valid, :), 1); % 1 x nChan
    aggRelPSD = aggAbsPSD ./ meanTotalPower(:);                 % nChan x nFreq (broadcast over frequency)

    entryA = struct();
    entryA.Subject = subjStr;
    entryA.InterventionNumber = intNum;
    entryA.Condition = condStr;
    entryA.ChannelLabels = chanLabelsRef;
    entryA.FrequenciesHz = freqVecRef;
    entryA.AggregatedAbsolutePSD = aggAbsPSD;   % nChan x nFreq, uV^2/Hz
    entryA.AggregatedRelativePSD = aggRelPSD;   % nChan x nFreq, fraction of the Phase-2-defined 1-45 Hz total
    entryA.NumValidWindows = nValidThis;
    subject_psd = concat_struct_array(subject_psd, entryA);

    % ---- Structures B & C: band power, per channel and whole-scalp ----
    for b = 1:nBands
        band = bandNames{b};
        absVals = rec.AbsoluteBandPower.(band)(valid, :);   % nValid x nChan
        relVals = rec.RelativeBandPower.(band)(valid, :);   % nValid x nChan (per-window ratios)

        % --- per-channel (Structure B) ---
        absMean = mean(absVals, 1);
        relMean = absMean ./ meanTotalPower;
        relMeanOfRatios = mean(relVals, 1);
        absLog  = log(absMean);   % ln() of the recording-level mean; NaN/complex only if absMean<=0 (flagged below, never masked)
        if any(absMean <= 0)
            integrityIssues{end+1} = sprintf('Sub%s Int%g Band %s: non-positive AbsolutePower found in %d channel(s) -- LogAbsolutePower set to NaN there, not silently patched.', ...
                subjStr, intNum, band, sum(absMean <= 0)); %#ok<AGROW>
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
            bandPowerRows(end+1, :) = {subjStr, intNum, condStr, chanLabelsRef{ch}, band, ...
                absMean(ch), absLog(ch), relMean(ch), relMeanOfRatios(ch), ...
                absMedian(ch), absSD(ch), absIQR(ch), absSkew(ch), ...
                relMedian(ch), relSD(ch), relIQR(ch), relSkew(ch)}; %#ok<AGROW>
        end

        % --- whole-scalp (Structure C): mean of the per-channel means
        %     (mathematically identical to a flat pooled mean here, since
        %     every channel in a recording shares the same valid-window
        %     set -- documented in the explanation, not just asserted) ---
        globalAbs = mean(absMean);
        globalRel = mean(relMean);
        globalRelMoR = mean(relMeanOfRatios);
        globalLog = log(globalAbs);
        if globalAbs <= 0
            integrityIssues{end+1} = sprintf('Sub%s Int%g Band %s: non-positive whole-scalp AbsolutePower -- LogAbsolutePower set to NaN.', subjStr, intNum, band); %#ok<AGROW>
            globalLog = NaN;
        end
        globalAbsMedian = median(absMean);   % median-of-channel-means (secondary, not a flat pooled median)
        globalRelMedian = median(relMean);

        globalRows(end+1, :) = {subjStr, intNum, condStr, band, ...
            globalAbs, globalLog, globalRel, globalRelMoR, globalAbsMedian, globalRelMedian}; %#ok<AGROW>
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

fprintf('--- Data integrity checks ---\n');

% (a) Subject x InterventionNumber unique
subjIntKey = strcat(RecordingQCT.Subject, '_int', string(RecordingQCT.InterventionNumber));
if numel(unique(subjIntKey)) ~= height(RecordingQCT)
    integrityIssues{end+1} = 'DUPLICATE Subject x InterventionNumber combinations found in the subject-level output.';
end

% (b) exactly one Condition per recording (guaranteed by construction here,
%     and checked explicitly)
for i = 1:numel(subject_psd)
    if isempty(subject_psd(i).Condition) || ~ischar(subject_psd(i).Condition)
        integrityIssues{end+1} = sprintf('Sub%s Int%g: Condition is empty or non-char.', subject_psd(i).Subject, subject_psd(i).InterventionNumber); %#ok<AGROW>
    end
end

% (c) no subject receives two different conditions for the SAME intervention number
[uSubjInt, ~, uIdx] = unique(strcat(RecordingQCT.Subject, '_int', string(RecordingQCT.InterventionNumber)));
for k = 1:numel(uSubjInt)
    condsHere = unique(RecordingQCT.Condition(uIdx == k));
    if numel(condsHere) > 1
        integrityIssues{end+1} = sprintf('%s: multiple different Condition values found (%s) for the same Subject x InterventionNumber.', ...
            uSubjInt{k}, strjoin(condsHere, ', ')); %#ok<AGROW>
    end
end

metaFile = load(META_MAT_FILE, 'T');
MetaT = metaFile.T;
expN    = height(MetaT);
expCond = struct('E', sum(strcmp(MetaT.Condition, 'E')), 'N', sum(strcmp(MetaT.Condition, 'N')), ...
                 'R', sum(strcmp(MetaT.Condition, 'R')));
expInt  = struct('Int1', sum(MetaT.InterventionNumber == 1), 'Int2', sum(MetaT.InterventionNumber == 2));

nSubjLevel = numel(subject_psd);
if nSubjLevel ~= expN
    integrityIssues{end+1} = sprintf('Expected %d subject-level recordings (Phase 1), found %d.', expN, nSubjLevel);
end
condCountsCheck = struct('E', sum(strcmp(RecordingQCT.Condition,'E')), ...
                          'N', sum(strcmp(RecordingQCT.Condition,'N')), ...
                          'R', sum(strcmp(RecordingQCT.Condition,'R')));
intCountsCheck = struct('Int1', sum(RecordingQCT.InterventionNumber==1), 'Int2', sum(RecordingQCT.InterventionNumber==2));
if ~isequal(condCountsCheck, expCond)
    integrityIssues{end+1} = sprintf('Condition counts do not match Phase 1 (expected E=%d,N=%d,R=%d; found E=%d,N=%d,R=%d).', ...
        expCond.E, expCond.N, expCond.R, condCountsCheck.E, condCountsCheck.N, condCountsCheck.R);
end
if ~isequal(intCountsCheck, expInt)
    integrityIssues{end+1} = sprintf('Intervention-position counts do not match Phase 1 (expected %d/%d; found %d/%d).', ...
        expInt.Int1, expInt.Int2, intCountsCheck.Int1, intCountsCheck.Int2);
end

nNaNInf_Abs = sum(~isfinite(BandPowerT.AbsolutePower));
nNaNInf_Rel = sum(~isfinite(BandPowerT.RelativePower));
nNaNInf_Log = sum(~isfinite(BandPowerT.LogAbsolutePower));
if nNaNInf_Abs > 0
    integrityIssues{end+1} = sprintf('%d non-finite AbsolutePower values found in BandPowerT.', nNaNInf_Abs);
end
if nNaNInf_Rel > 0
    integrityIssues{end+1} = sprintf('%d non-finite RelativePower values found in BandPowerT.', nNaNInf_Rel);
end

% (g) valid-window counts preserved (spot check: RecordingQCT vs Phase 2 log totals)
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

% --- Diagnostic 1: distribution of participant-level (whole-scalp) band power ---
try
    fig1 = figure('Visible', 'off', 'Position', [100 100 1400 700]);
    nPrimary = numel(PARAMS3.PrimaryBands);
    for b = 1:nPrimary
        subplot(2, ceil(nPrimary/2), b);
        band = PARAMS3.PrimaryBands{b};
        vals = GlobalBandPowerT.AbsolutePower(strcmp(GlobalBandPowerT.Band, band));
        histogram(vals);
        title(sprintf('%s (whole-scalp AbsolutePower)', band), 'Interpreter', 'none');
        xlabel('uV^2'); ylabel('Recordings');
    end
    sgtitle('Diagnostic 1: Distribution of subject-level band power (primary bands)');
    saveas(fig1, fullfile(SUBJ_DIR, 'diagnostic1_subject_bandpower_distributions.png'));
    close(fig1);
catch ME1
    warning('Could not generate Diagnostic 1 figure: %s', ME1.message);
end

% --- Diagnostic 2: mean vs median ACROSS CHANNELS, whole-scalp, all
%     recordings x all bands (uses GlobalBandPowerT's mean vs
%     median-of-channel-means, i.e. the two central-tendency estimates
%     already computed above) ---
try
    fig2 = figure('Visible', 'off');
    scatter(GlobalBandPowerT.AbsolutePower_MedianAcrossChannels, GlobalBandPowerT.AbsolutePower, 15, 'filled', 'MarkerFaceAlpha', 0.4);
    hold on;
    lims = [0, max([GlobalBandPowerT.AbsolutePower; GlobalBandPowerT.AbsolutePower_MedianAcrossChannels])*1.05];
    plot(lims, lims, 'k--');
    hold off;
    xlabel('Median-across-channels AbsolutePower (uV^2)');
    ylabel('Mean-across-channels AbsolutePower (uV^2)');
    title('Diagnostic 2: Mean vs median ACROSS CHANNELS (whole-scalp, all recordings x bands)');
    axis(repmat(lims,1,2));
    saveas(fig2, fullfile(SUBJ_DIR, 'diagnostic2_mean_vs_median_across_channels.png'));
    close(fig2);
catch ME2
    warning('Could not generate Diagnostic 2 figure: %s', ME2.message);
end

% --- Diagnostic 3: valid windows per recording, by condition ---
try
    fig3 = figure('Visible', 'off');
    conds = {'E', 'N', 'R'};
    hold on;
    for c = 1:numel(conds)
        vals = RecordingQCT.ValidWindows(strcmp(RecordingQCT.Condition, conds{c}));
        jitterX = c + (rand(size(vals))-0.5)*0.3;
        scatter(jitterX, vals, 25, 'filled');
        errorbar(c, mean(vals), std(vals), 'k', 'LineWidth', 1.5, 'CapSize', 12);
    end
    hold off;
    xticks(1:numel(conds)); xticklabels(conds); xlim([0.5 numel(conds)+0.5]);
    ylabel('Valid windows per recording');
    title('Diagnostic 3: Valid-window counts by condition (descriptive, not a test)');
    saveas(fig3, fullfile(SUBJ_DIR, 'diagnostic3_valid_windows_by_condition.png'));
    close(fig3);
catch ME3
    warning('Could not generate Diagnostic 3 figure: %s', ME3.message);
end

% --- Diagnostic 4: participant-level alpha power by condition (descriptive
%     sanity check only -- R expected > E/N; NOT a hypothesis test) ---
try
    fig4 = figure('Visible', 'off');
    conds = {'E', 'N', 'R'};
    hold on;
    for c = 1:numel(conds)
        vals = GlobalBandPowerT.AbsolutePower(strcmp(GlobalBandPowerT.Band, 'Alpha') & strcmp(GlobalBandPowerT.Condition, conds{c}));
        jitterX = c + (rand(size(vals))-0.5)*0.3;
        scatter(jitterX, vals, 25, 'filled');
        errorbar(c, mean(vals), std(vals)/sqrt(numel(vals)), 'k', 'LineWidth', 1.5, 'CapSize', 12);
    end
    hold off;
    xticks(1:numel(conds)); xticklabels(conds); xlim([0.5 numel(conds)+0.5]);
    ylabel('Whole-scalp Alpha AbsolutePower (uV^2)');
    title({'Diagnostic 4: Alpha power by condition -- DESCRIPTIVE ONLY', 'error bars = mean +/- SEM; no statistical test performed'});
    saveas(fig4, fullfile(SUBJ_DIR, 'diagnostic4_alpha_power_by_condition.png'));
    close(fig4);
catch ME4
    warning('Could not generate Diagnostic 4 figure: %s', ME4.message);
end

extremeRows = {};
condsForFlags = {'E', 'N', 'R'};
for b = 1:nBands
    band = bandNames{b};
    for c = 1:numel(condsForFlags)
        sel   = strcmp(GlobalBandPowerT.Band, band) & strcmp(GlobalBandPowerT.Condition, condsForFlags{c});
        vals  = GlobalBandPowerT.LogAbsolutePower(sel);
        rawV  = GlobalBandPowerT.AbsolutePower(sel);
        subjs = GlobalBandPowerT.Subject(sel);
        ints  = GlobalBandPowerT.InterventionNumber(sel);
        conds = GlobalBandPowerT.Condition(sel);
        med = median(vals);
        madVal = median(abs(vals - med));
        robustSD = 1.4826 * madVal;
        if robustSD == 0 || isnan(robustSD)
            continue;
        end
        z = (vals - med) / robustSD;
        flagged = find(abs(z) > 3);
        for f = flagged(:)'
            extremeRows(end+1, :) = {subjs{f}, ints(f), conds{f}, band, rawV(f), vals(f), z(f)}; %#ok<AGROW>
        end
    end
end
if isempty(extremeRows)
    ExtremeT = cell2table(cell(0,7), 'VariableNames', {'Subject','InterventionNumber','Condition','Band','AbsolutePower','LogAbsolutePower','RobustZ'});
else
    ExtremeT = cell2table(extremeRows, 'VariableNames', {'Subject','InterventionNumber','Condition','Band','AbsolutePower','LogAbsolutePower','RobustZ'});
end
ExtremeT = normalize_text_columns(ExtremeT);
writetable(ExtremeT, fullfile(SUBJ_DIR, 'diagnostic5_extreme_value_flags.csv'));

%% ------------------------------------------------------------------
%  8. FINAL PHASE 3 REPORT
% ------------------------------------------------------------------

report_lines = {};
report_lines{end+1} = '========================================================';
report_lines{end+1} = 'PHASE 3 -- SUBJECT-LEVEL AGGREGATION REPORT';
report_lines{end+1} = '========================================================';
report_lines{end+1} = '';
report_lines{end+1} = sprintf('Input recordings          : %d', nRec);
report_lines{end+1} = sprintf('Subject-level recordings  : %d', nSubjLevel);
report_lines{end+1} = '';
report_lines{end+1} = sprintf('E : %d', condCountsCheck.E);
report_lines{end+1} = sprintf('N : %d', condCountsCheck.N);
report_lines{end+1} = sprintf('R : %d', condCountsCheck.R);
report_lines{end+1} = '';
report_lines{end+1} = sprintf('Intervention 1 : %d', intCountsCheck.Int1);
report_lines{end+1} = sprintf('Intervention 2 : %d', intCountsCheck.Int2);
report_lines{end+1} = '';
report_lines{end+1} = sprintf('Aggregation method: %s (window-level); LogAbsolutePower = ln(recording-level AbsolutePower)', ...
    PARAMS3.WindowAggregationStatistic);
report_lines{end+1} = '';
report_lines{end+1} = sprintf('Absolute power transformation: %s -- %s', ...
    iif(PARAMS3.AbsolutePowerLogTransform, 'ln (natural log) applied', 'none applied'), PARAMS3.LogTransformConvention);
report_lines{end+1} = '';
report_lines{end+1} = 'Valid-window summary (per recording):';
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
report_lines{end+1} = sprintf('Relative power definition: %s', PARAMS3.RelativePowerDefinition);
report_lines{end+1} = '';
report_lines{end+1} = sprintf('Extreme value flags (Diagnostic 5, flagged not removed): %d', height(ExtremeT));
report_lines{end+1} = sprintf('  Rule: %s', PARAMS3.ExtremeValueRule);
report_lines{end+1} = sprintf('  Flags by condition: E=%d, N=%d, R=%d', sum(strcmp(ExtremeT.Condition,'E')), ...
    sum(strcmp(ExtremeT.Condition,'N')), sum(strcmp(ExtremeT.Condition,'R')));
report_lines{end+1} = sprintf('Windows rejected for boundary (all recordings): %d', sum(RecordingQCT.RejectedBoundary(~isnan(RecordingQCT.RejectedBoundary))));
report_lines{end+1} = '';
report_lines{end+1} = 'Output files:';
report_lines{end+1} = sprintf('  %s', fullfile(SUBJ_DIR, 'phase3_subject_level.mat'));
report_lines{end+1} = sprintf('  %s', fullfile(SUBJ_DIR, 'phase3_subject_level_bandpower.mat'));
report_lines{end+1} = sprintf('  %s', fullfile(SUBJ_DIR, 'phase3_recording_summary.xlsx'));
report_lines{end+1} = sprintf('  %s', fullfile(SUBJ_DIR, 'phase3_parameters.mat'));
report_lines{end+1} = sprintf('  %s', fullfile(SUBJ_DIR, 'diagnostic1_subject_bandpower_distributions.png'));
report_lines{end+1} = sprintf('  %s', fullfile(SUBJ_DIR, 'diagnostic2_mean_vs_median_across_channels.png'));
report_lines{end+1} = sprintf('  %s', fullfile(SUBJ_DIR, 'diagnostic3_valid_windows_by_condition.png'));
report_lines{end+1} = sprintf('  %s', fullfile(SUBJ_DIR, 'diagnostic4_alpha_power_by_condition.png'));
report_lines{end+1} = sprintf('  %s', fullfile(SUBJ_DIR, 'diagnostic5_extreme_value_flags.csv'));
report_lines{end+1} = sprintf('  %s', fullfile(REPORT_DIR, 'phase3_aggregation_report.txt'));
report_lines{end+1} = '========================================================';

fullReport = strjoin(report_lines, '\n');
fprintf('%s\n', fullReport);

fid = fopen(fullfile(REPORT_DIR, 'phase3_aggregation_report.txt'), 'w');
fprintf(fid, '%s\n', fullReport);
fclose(fid);

fprintf('\nPhase 3 complete. Do not proceed to Phase 4 until this output has been reviewed.\n');


%% ========================================================================
%  LOCAL FUNCTIONS
%  (must follow ALL top-level executable statements per MATLAB's rule for
%  local functions in a script file)
% ========================================================================

function outArr = concat_struct_array(inArr, entry)
% Appends 'entry' (a scalar struct) to 'inArr' (a struct array, possibly
% struct([]) i.e. 0x0 with no fields yet).
    if isempty(inArr)
        outArr = entry;
    else
        outArr = inArr;
        outArr(end+1) = entry; %#ok<AGROW>
    end
end

function T = normalize_text_columns(T)
% Forces every char-typed table variable to cellstr, guarding against the
% char-matrix collapse discovered in Phase 0/1/2 (see those scripts).
    varNames = T.Properties.VariableNames;
    for v = 1:numel(varNames)
        col = T.(varNames{v});
        if ischar(col)
            T.(varNames{v}) = cellstr(col);
        end
    end
end

function v = local_iqr(X)
% Column-wise interquartile range (75th - 25th percentile), computed via
% linear-interpolation percentiles (MATLAB's default "R-7"-equivalent
% method) WITHOUT requiring the Statistics and Machine Learning Toolbox.
% X: nObs x nCols.
    v = local_prctile(X, 75) - local_prctile(X, 25);
end

function p = local_prctile(X, pct)
% Column-wise percentile via linear interpolation on the sorted data
% (equivalent to MATLAB's default prctile method), implemented locally to
% avoid a Statistics and Machine Learning Toolbox dependency.
% X: nObs x nCols. pct: scalar in [0,100].
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
    % MATLAB's default percentile method: linear interpolation between
    % order statistics at positions (i - 0.5)/n.
    posPct = (pct/100) * n + 0.5;
    posPct = min(max(posPct, 1), n);
    lo = floor(posPct);
    hi = ceil(posPct);
    frac = posPct - lo;
    p = Xs(lo, :) + frac .* (Xs(hi, :) - Xs(lo, :));
end

function s = local_skewness(X)
% Column-wise population (biased) Fisher-Pearson skewness coefficient,
% implemented locally to avoid a Statistics and Machine Learning Toolbox
% dependency. X: nObs x nCols.
    n = size(X, 1);
    mu = mean(X, 1);
    sigma = std(X, 1, 1);   % population (N-normalized) std, matching the classic skewness definition
    sigma(sigma == 0) = NaN;
    s = mean((X - mu).^3, 1) ./ (sigma.^3);
end

function out = iif(cond, a, b)
% Inline conditional helper for cleaner report-string assembly.
    if cond
        out = a;
    else
        out = b;
    end
end
