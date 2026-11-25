%% ========================================================================
% X03_compute_intervention_psd.m
%
% Window segmentation, window rejection and power spectra for every recording.
% 4-s Hamming windows (2048 samples) with 50% overlap; NFFT = 2048 (0.25 Hz resolution).
% A window is rejected if it contains a non-finite value, |amplitude| > 150 uV on any
% channel, a channel with SD < 0.01 uV, or a 'boundary' event (removed data).
% Retained frequencies: 1-45 Hz and 55-100 Hz (45-55 Hz guard band around the notch).
% Inputs : 02_metadata, DATA_DIR recordings.
% Outputs: 03_psd/intervention_psd.mat (per-window and mean PSD), processing log, 09_reports/.
% ========================================================================

clear; clc;

%% ------------------------------------------------------------------

% >>> Must match the ANALYSIS_ROOT produced by 02_build_intervention_metadata.m <<<
config;
ANALYSIS_ROOT = OUT_DIR;

META_DIR   = fullfile(ANALYSIS_ROOT, '02_metadata');
PSD_DIR    = fullfile(ANALYSIS_ROOT, '03_psd');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');
if ~exist(PSD_DIR, 'dir');    mkdir(PSD_DIR);    end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

META_MAT_FILE = fullfile(META_DIR, 'intervention_metadata.mat');

fprintf('========================================================\n');
fprintf('PHASE 2 -- CONTINUOUS EEG -> PSD\n');
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
           '02_build_intervention_metadata.m first.'], META_MAT_FILE);
end

%% ------------------------------------------------------------------
%  1. SPECTRAL-ANALYSIS PARAMETERS (saved verbatim for reproducibility)
% ------------------------------------------------------------------

PARAMS = struct();
PARAMS.ScriptName          = '03_compute_intervention_psd.m';
PARAMS.GeneratedOn         = datestr(now);

PARAMS.ExpectedSamplingRateHz = 512;   % cross-checked per recording
PARAMS.ExpectedNumChannels    = 62;    % cross-checked per recording

PARAMS.WindowLengthSec        = 4;                 % primary analysis window
PARAMS.WindowOverlapFraction  = 0.5;               % 50% overlap
PARAMS.WindowLengthSamples    = round(PARAMS.WindowLengthSec * PARAMS.ExpectedSamplingRateHz);
PARAMS.WindowStepSamples      = round(PARAMS.WindowLengthSamples * (1 - PARAMS.WindowOverlapFraction));

PARAMS.TaperType              = 'hamming';          % applied to the full window
PARAMS.NFFT                   = PARAMS.WindowLengthSamples;  % 2048, already a power of 2
PARAMS.FrequencyResolutionHz  = PARAMS.ExpectedSamplingRateHz / PARAMS.NFFT;  % 0.25 Hz

PARAMS.PrimaryAnalysisFreqRangeHz = [1 45];   % used for the relative-power denominator (unchanged rationale below)
PARAMS.HighGammaFreqRangeHz       = [55 100];
PARAMS.NotchGuardBandExcludedHz   = [45 55];  % EXCLUDED entirely: wider than the actual 47-53 Hz notch, giving a
PARAMS.RetainedFreqRangesHz = {PARAMS.PrimaryAnalysisFreqRangeHz, PARAMS.HighGammaFreqRangeHz};
% Rationale: 1 Hz matches preprocessing's own high-pass edge (no content
% below 1 Hz exists in these recordings); 45 Hz matches the freqrange
% used by preprocessing's own pop_rejchan step and gives a safety margin
% before the 47-53 Hz notch; 100 Hz matches preprocessing's own low-pass
% edge (no higher content exists in these recordings). The 45-55 Hz gap
% is deliberately wider than the notch itself and is fully excluded from
% the retained frequency-resolved PSD, from every band, and from the
% relative-power denominator.

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
PARAMS.PrimaryBands    = {'Delta', 'Theta', 'Alpha', 'Beta'};
PARAMS.ExploratoryBands = {'Beta_low_exploratory', 'Beta_high_exploratory', ...
    'Gamma_low_exploratory', 'Gamma_high1_exploratory', 'Gamma_high2_exploratory'};
PARAMS.BandEdgeConvention = 'both edges inclusive at every boundary; adjacent bands share their boundary frequency sample (exact under trapz, not a double-count)';

PARAMS.RelativePowerDenominatorRangeHz = PARAMS.PrimaryAnalysisFreqRangeHz;

PARAMS.WindowQC = struct( ...
    'RejectIfNaNInf',            true, ...
    'AmplitudeThresholdUV',      150, ...   % reject window if any sample, any channel, |amplitude| exceeds this
    'FlatChannelStdThresholdUV', 0.01);     % reject window if any channel's std within the window is below this (dead/flat channel)

PARAMS.WindowQC.RejectIfSpansBoundary = true;
PARAMS.WindowQC.BoundaryEventTypes = {'boundary'};   % case-insensitive match on event type
PARAMS.WindowQC.MaxFractionRejectedByBoundary = 0.50;
PARAMS.DurationMismatchToleranceSec = 1 / PARAMS.ExpectedSamplingRateHz;

PARAMS.AssumedDataUnits = 'microvolts (uV), per standard EEGLAB convention; not independently re-verified against amplifier calibration records.';

PARAMS.DurationHandlingPolicy = [ ...
    'Use ALL available analysis windows from each recording''s full ' ...
    'duration; do NOT truncate recordings to a common/shortest duration. ' ...
    'Unequal per-recording window counts are NOT allowed to bias later ' ...
    'comparisons because Phase 3 aggregates every subject x condition ' ...
    'to exactly one value before any group comparison; a participant ' ...
    'with more windows contributes more PRECISION to their own subject- ' ...
    'level estimate, never more WEIGHT in the between-subject analysis.'];

PARAMS.SaveFullWindowLevelPSD = false;
% If true, the full (window x channel x frequency) PSD array is ALSO
% saved per recording (large: for this dataset, roughly 1-1.5 GB total
% across all recordings). Default false: only per-window BAND POWER
% (small, and sufficient for the primary analysis, ROI, and topography)
% and a per-recording MEAN frequency-resolved PSD (also small, and
% sufficient for spectral-shape diagnostics/exploratory work) are saved.
% Set to true and re-run if full window-level frequency resolution is
% needed later (e.g. for individual alpha peak frequency estimation).

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
ref_labels = loadedMeta.ref_labels(:);   % 62x1 cellstr, Phase 1 reference channel order
nRefChan = numel(ref_labels);

nRec = height(T);
fprintf('Loaded Phase 1 metadata: %d recordings, %d reference channels.\n\n', nRec, nRefChan);

%% ------------------------------------------------------------------
%  3. SHARED, PRECOMPUTED QUANTITIES (identical for every recording,
%     since Fs is validated == 512 Hz for all recordings in Phase 1)
% ------------------------------------------------------------------

winSamples  = PARAMS.WindowLengthSamples;
stepSamples = PARAMS.WindowStepSamples;
NFFT        = PARAMS.NFFT;
taperWin    = hamming(winSamples);   % column vector, winSamples x 1

% One-sided frequency vector exactly as pwelch would produce it for real
% input of this NFFT/Fs (documented formula, avoids an extra dummy call).
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

psd_data = struct([]);   % one entry per SUCCESSFULLY processed recording
kSuccess = 0;

logCols = {'Subject', 'InterventionNumber', 'Condition', 'DurationSec', ...
           'NumCandidateWindows', 'NumValidWindows', 'PercentValid', ...
           'AnalysisDurationUsedSec', 'RejectedNaNInf', 'RejectedAmplitude', ...
           'RejectedFlat', 'RejectedBoundary', 'NumBoundaryEvents', 'DurationLoadedSec', ...
           'EventTypes', 'ProcessingStatus', 'ErrorMessage'};
log_rows = cell(0, numel(logCols));

for i = 1:nRec

    subjStr   = getstr(T.Subject{i});
    intNum    = T.InterventionNumber(i);
    condStr   = getstr(T.Condition{i});
    fpath     = getstr(T.FullPath{i});
    durPhase1 = T.DurationSec(i);

    fprintf('[%d/%d] Subject %s, Intervention %g, Condition %s ... ', i, nRec, subjStr, intNum, condStr);

    status = 'OK';
    errMsg = '';
    nCandidate = 0; nValid = 0; pctValid = NaN; analysisDurUsed = NaN;
    nRejNaNInf = 0; nRejAmp = 0; nRejFlat = 0;
    nRejBoundary = 0; nBoundary = NaN; durLoaded = NaN; evtSummary = '';

    try
        % ---- Load the already-preprocessed continuous recording ----
        [fdir, fname, fext] = fileparts(fpath);
        EEG = pop_loadset('filename', [fname fext], 'filepath', fdir);

        Fs = EEG.srate;
        if Fs ~= PARAMS.ExpectedSamplingRateHz
            error('Sampling rate mismatch: expected %d Hz, found %g Hz.', ...
                PARAMS.ExpectedSamplingRateHz, Fs);
        end

        durLoaded = EEG.pnts / Fs;
        if ~isnan(durPhase1) && abs(durLoaded - durPhase1) > PARAMS.DurationMismatchToleranceSec
            error(['Duration mismatch vs Phase 1 metadata (loaded %.3f s, Phase 1 %.3f s). ' ...
                   'The data folder changed after Phase 1 -- rerun 02_build_intervention_metadata.m.'], ...
                   durLoaded, durPhase1);
        end

        evtSummary = summarize_event_types(EEG);
        bLat = get_boundary_latencies(EEG, PARAMS.WindowQC.BoundaryEventTypes);
        nBoundary = numel(bLat);

        chanLabels = {EEG.chanlocs.labels};
        if numel(chanLabels) ~= nRefChan
            error('Channel count mismatch: expected %d, found %d.', nRefChan, numel(chanLabels));
        end

        % ---- Reorder channels (bookkeeping only, no value change) to
        %      match the Phase 1 reference channel order, so the channel
        %      dimension is aligned identically across every recording. ----
        reorderIdx = match_channel_order(chanLabels, ref_labels);

        data = double(EEG.data(reorderIdx, :));   % nChan x nSamples, uV, reordered
        nSamplesRec = size(data, 2);
        nChan = size(data, 1);

        if nSamplesRec < winSamples
            error('Recording shorter than one analysis window (%d samples < %d).', ...
                nSamplesRec, winSamples);
        end

        % ---- Candidate window start indices (sample space; the ~1.7 s
        %      artificial preprocessing events are NOT consulted here) ----
        starts = 0:stepSamples:(nSamplesRec - winSamples);
        nCandidate = numel(starts);

        % ---- Preallocate per-window outputs ----
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

        sumPSD = zeros(nFreqSel, nChan);   % accumulator for the recording-level mean PSD

        for w = 1:nCandidate
            idx = (starts(w) + 1):(starts(w) + winSamples);
            seg = data(:, idx);          % nChan x winSamples
            segT = seg.';                % winSamples x nChan

            [isBad, reasonStr, isNaNInf, isAmp, isFlat] = check_window_qc( ...
                segT, PARAMS.WindowQC.AmplitudeThresholdUV, PARAMS.WindowQC.FlatChannelStdThresholdUV);

            isBoundary = PARAMS.WindowQC.RejectIfSpansBoundary && ~isempty(bLat) && ...
                any(bLat >= starts(w) + 1 & bLat <= starts(w) + winSamples);
            if isBoundary
                nRejBoundary = nRejBoundary + 1;
                isBad = true;
                if isempty(reasonStr); reasonStr = 'Boundary'; else; reasonStr = [reasonStr '+Boundary']; end
            end

            if isNaNInf; nRejNaNInf = nRejNaNInf + 1; end
            if isAmp;    nRejAmp    = nRejAmp    + 1; end
            if isFlat;   nRejFlat   = nRejFlat   + 1; end

            windowRejectReason{w} = reasonStr;

            if isBad
                windowValid(w) = false;
                continue;
            end
            windowValid(w) = true;

            % ---- Welch/periodogram PSD for this window, all channels ----
            [Pxx, ~] = pwelch(segT, taperWin, 0, NFFT, Fs);   % nFreqFull x nChan, one-sided, uV^2/Hz
            PxxSel = Pxx(freqMask, :);                        % nFreqSel x nChan

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
        analysisDurUsed = nValid * PARAMS.WindowLengthSec;  % sum of valid-window lengths; overlap double-counted by design (documented)

        if nValid > 0
            meanPSD = (sumPSD / nValid).';   % nChan x nFreqSel (transposed for storage convention below)
        else
            meanPSD = nan(nChan, nFreqSel);
            warning('Recording %s (sub%s int%g): 0 valid windows out of %d candidates.', ...
                fname, subjStr, intNum, nCandidate);
        end

        % ---- Assemble this recording's entry ----
        entry = struct();
        entry.Subject = subjStr;
        entry.InterventionNumber = intNum;
        entry.Condition = condStr;
        entry.ConditionLetterRaw = getstr(T.ConditionLetterRaw{i});
        entry.FullPath = fpath;
        entry.DurationSec = durLoaded;
        entry.NumBoundaryEvents = nBoundary;
        entry.NumWindowsRejectedBoundary = nRejBoundary;
        entry.ChannelLabels = ref_labels;              % nChan x 1 cellstr, common order
        entry.FrequenciesHz = freqVecSel;               % nFreqSel x 1
        entry.WindowStartSec = windowStartSec;
        entry.WindowEndSec = windowEndSec;
        entry.WindowValid = windowValid;
        entry.WindowRejectReason = windowRejectReason;
        entry.AbsoluteBandPower = absBandPower;          % struct of (nWin x nChan), uV^2
        entry.RelativeBandPower = relBandPower;          % struct of (nWin x nChan), unitless fraction
        entry.TotalPower_1to45Hz = totalPower;           % nWin x nChan, uV^2 (the relative-power denominator, saved explicitly)
        entry.MeanPSD = meanPSD;                         % nChan x nFreqSel, uV^2/Hz, mean over VALID windows only
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

    log_rows(end+1, :) = {subjStr, intNum, condStr, durPhase1, nCandidate, nValid, ...
        pctValid, analysisDurUsed, nRejNaNInf, nRejAmp, nRejFlat, nRejBoundary, nBoundary, durLoaded, ...
        evtSummary, status, errMsg}; %#ok<AGROW>
end

kSuccess = numel(psd_data);   % re-derived from the actual array, not just the running counter

fracBoundary = cellfun(@(nb, nc) nb / max(nc, 1), log_rows(:, 12), log_rows(:, 5));
highBoundary = fracBoundary > PARAMS.WindowQC.MaxFractionRejectedByBoundary;
if median(fracBoundary) > PARAMS.WindowQC.MaxFractionRejectedByBoundary
    error(['The MEDIAN recording lost %.0f%% of its windows to boundary rejection (limit %.0f%%). ' ...
           'This dataset-wide pattern means boundaries sit at regular epoch junctions ' ...
           '(e.g. eeg_epoch2continuous), not at isolated splices. Nothing was saved.'], ...
           100 * median(fracBoundary), 100 * PARAMS.WindowQC.MaxFractionRejectedByBoundary);
end
for k = find(highBoundary(:))'
    fprintf(2, '  FLAG (kept, not excluded): sub%s int%g lost %d/%d windows to boundaries (%d boundary events).\n', ...
        log_rows{k, 1}, log_rows{k, 2}, log_rows{k, 12}, log_rows{k, 5}, log_rows{k, 13});
end

%% ------------------------------------------------------------------
%  5. SAVE OUTPUTS
% ------------------------------------------------------------------

save(fullfile(PSD_DIR, 'phase2_parameters.mat'), 'PARAMS');

FrequencyVectorHz = freqVecSel;
ChannelLabelsReference = ref_labels;
save(fullfile(PSD_DIR, 'intervention_psd.mat'), ...
    'psd_data', 'FrequencyVectorHz', 'ChannelLabelsReference', 'PARAMS', '-v7.3');

ProcessingLogT = cell2table(log_rows, 'VariableNames', logCols);
ProcessingLogT = normalize_text_columns(ProcessingLogT);
ProcessingLogT.FractionRejectedBoundary = fracBoundary(:);
ProcessingLogT.HighBoundaryRejection = highBoundary(:);
writetable(ProcessingLogT, fullfile(PSD_DIR, 'phase2_processing_log.csv'));
writetable(ProcessingLogT, fullfile(PSD_DIR, 'phase2_processing_log.xlsx'));
save(fullfile(PSD_DIR, 'phase2_processing_log.mat'), 'ProcessingLogT');

%% ------------------------------------------------------------------
%  6. DIAGNOSTICS
% ------------------------------------------------------------------

nConsistencyFailures = sum(strcmp(ProcessingLogT.ProcessingStatus, 'FAILED') & ...
    contains(ProcessingLogT.ErrorMessage, {'Sampling rate mismatch', 'Channel count mismatch', 'Channel labels'}));

% --- Diagnostic 2: distribution of valid-window counts across recordings ---
validCounts = arrayfun(@(x) x.NumValidWindows, psd_data);
try
    fig2 = figure('Visible', 'off');
    histogram(validCounts);
    xlabel('Number of valid windows per recording');
    ylabel('Number of recordings');
    title('Diagnostic 2: Distribution of valid-window counts');
    saveas(fig2, fullfile(PSD_DIR, 'diagnostic2_valid_window_count_hist.png'));
    close(fig2);
catch ME2
    warning('Could not generate Diagnostic 2 figure: %s', ME2.message);
end

% --- Diagnostic 1: example PSD curves, one representative recording per
%     condition, a few representative channels ---
try
    exampleChans = {'Fz', 'Cz', 'Oz'};
    exampleChanIdx = find(ismember(ref_labels, exampleChans));
    if isempty(exampleChanIdx)
        exampleChanIdx = 1:min(3, numel(ref_labels));
    end
    condsToPlot = {'E', 'N', 'R'};
    fig1 = figure('Visible', 'off', 'Position', [100 100 1200 400]);
    for c = 1:numel(condsToPlot)
        subplot(1, numel(condsToPlot), c);
        condMask = arrayfun(@(x) strcmp(x.Condition, condsToPlot{c}), psd_data);
        condIdxList = find(condMask);
        if isempty(condIdxList)
            title(sprintf('%s (no successful recordings)', condsToPlot{c}));
            continue;
        end
        exRec = psd_data(condIdxList(1));
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
            plot(plotFreq_g, plotPSD, 'LineWidth', 1.2);
        end
        hold off;
        xlabel('Frequency (Hz)');
        ylabel('PSD (dB, 10*log10 uV^2/Hz)');
        title(sprintf('%s: Subject %s, Int %g', condsToPlot{c}, exRec.Subject, exRec.InterventionNumber));
        legend(ref_labels(exampleChanIdx), 'Location', 'northeast');
        grid on;
    end
    sgtitle('Diagnostic 1: Example mean PSD curves (representative recordings/channels)');
    saveas(fig1, fullfile(PSD_DIR, 'diagnostic1_example_psd_curves.png'));
    close(fig1);
catch ME1
    warning('Could not generate Diagnostic 1 figure: %s', ME1.message);
end

% --- Diagnostic 4: basic descriptive summary of absolute/relative band
%     power. POOLED and broken out by condition for sanity-checking only
%     -- purely descriptive, NOT an inferential comparison, NOT a figure. ---
diag4_rows = {};
condsForSummary = {'ALL', 'E', 'N', 'R'};
for c = 1:numel(condsForSummary)
    cLabel = condsForSummary{c};
    if strcmp(cLabel, 'ALL')
        recIdx = 1:numel(psd_data);
    else
        recIdx = find(arrayfun(@(x) strcmp(x.Condition, cLabel), psd_data));
    end
    for b = 1:nBands
        allAbs = []; allRel = [];
        for r = recIdx
            v = psd_data(r).WindowValid;
            absVals = psd_data(r).AbsoluteBandPower.(bandNames{b})(v, :);
            relVals = psd_data(r).RelativeBandPower.(bandNames{b})(v, :);
            allAbs = [allAbs; absVals(:)]; %#ok<AGROW>
            allRel = [allRel; relVals(:)]; %#ok<AGROW>
        end
        diag4_rows(end+1, :) = {cLabel, bandNames{b}, mean(allAbs, 'omitnan'), std(allAbs, 'omitnan'), ...
            mean(allRel, 'omitnan'), std(allRel, 'omitnan'), numel(allAbs)}; %#ok<AGROW>
    end
end
Diagnostic4T = cell2table(diag4_rows, 'VariableNames', ...
    {'ConditionGroup', 'Band', 'MeanAbsolutePower_uV2', 'SDAbsolutePower_uV2', ...
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

condCounts = struct('E', 0, 'N', 0, 'R', 0, 'Unknown', 0);
for r = 1:numel(psd_data)
    c = psd_data(r).Condition;
    if ~isempty(c) && isfield(condCounts, c)
        condCounts.(c) = condCounts.(c) + 1;
    else
        condCounts.Unknown = condCounts.Unknown + 1;
    end
end

report_lines = {};
report_lines{end+1} = '========================================================';
report_lines{end+1} = 'PHASE 2 -- CONTINUOUS EEG -> PSD REPORT';
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
report_lines{end+1} = sprintf('Relative-power denominator: %g-%g Hz only (unchanged; high-gamma bands are extra fractions of this same baseline)', ...
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
report_lines{end+1} = sprintf('  Total windows rejected for boundary  : %d  (windows containing a ''boundary'' splice event)', sum(ProcessingLogT.RejectedBoundary));
report_lines{end+1} = sprintf('  Recordings containing >=1 boundary event: %d / %d (total boundary events: %d)', ...
    sum(ProcessingLogT.NumBoundaryEvents > 0), height(ProcessingLogT), sum(ProcessingLogT.NumBoundaryEvents(~isnan(ProcessingLogT.NumBoundaryEvents))));
report_lines{end+1} = '  NOTE: a window can be rejected for more than one reason, so these counts can overlap.';
report_lines{end+1} = sprintf('  Median fraction of windows lost to boundaries per recording: %.1f %%', 100 * median(ProcessingLogT.FractionRejectedBoundary));
hb = find(ProcessingLogT.HighBoundaryRejection);
report_lines{end+1} = sprintf('  Recordings flagged HighBoundaryRejection (> %.0f %% of windows lost; KEPT, not excluded): %d', ...
    100 * PARAMS.WindowQC.MaxFractionRejectedByBoundary, numel(hb));
for k = hb(:)'
    report_lines{end+1} = sprintf('    sub%s int%g (%s): %d valid of %d windows, %d boundary events', ProcessingLogT.Subject{k}, ...
        ProcessingLogT.InterventionNumber(k), ProcessingLogT.Condition{k}, ProcessingLogT.NumValidWindows(k), ...
        ProcessingLogT.NumCandidateWindows(k), ProcessingLogT.NumBoundaryEvents(k)); %#ok<AGROW>
end
report_lines{end+1} = '';
report_lines{end+1} = 'Event types found per recording (for inspection):';
for r = 1:height(ProcessingLogT)
    report_lines{end+1} = sprintf('  sub%s int%g (%s): %s', ProcessingLogT.Subject{r}, ProcessingLogT.InterventionNumber(r), ...
        ProcessingLogT.Condition{r}, ProcessingLogT.EventTypes{r}); %#ok<AGROW>
end
report_lines{end+1} = '';
report_lines{end+1} = 'Condition counts processed:';
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
report_lines{end+1} = sprintf('  %s', fullfile(PSD_DIR, 'intervention_psd.mat'));
report_lines{end+1} = sprintf('  %s', fullfile(PSD_DIR, 'phase2_processing_log.csv'));
report_lines{end+1} = sprintf('  %s', fullfile(PSD_DIR, 'phase2_processing_log.xlsx'));
report_lines{end+1} = sprintf('  %s', fullfile(PSD_DIR, 'diagnostic1_example_psd_curves.png'));
report_lines{end+1} = sprintf('  %s', fullfile(PSD_DIR, 'diagnostic2_valid_window_count_hist.png'));
report_lines{end+1} = sprintf('  %s', fullfile(PSD_DIR, 'diagnostic4_band_power_summary.csv'));
report_lines{end+1} = sprintf('  %s', fullfile(REPORT_DIR, 'phase2_psd_report.txt'));
report_lines{end+1} = '========================================================';

fullReport = strjoin(report_lines, '\n');
fprintf('%s\n', fullReport);

fid = fopen(fullfile(REPORT_DIR, 'phase2_psd_report.txt'), 'w');
fprintf(fid, '%s\n', fullReport);
fclose(fid);

fprintf('\nPhase 2 complete. Do not proceed to Phase 3 until this output has been reviewed.\n');


%% ========================================================================
%  LOCAL FUNCTIONS
%  (must follow ALL top-level executable statements per MATLAB's rule for
%  local functions in a script file)
% ========================================================================

function s = getstr(x)
% Defensive unwrap: returns a plain char regardless of whether the input
% arrived as char, a 1x1 cell containing char, a string scalar, or numeric.
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
% Appends 'entry' (a scalar struct) to 'inArr' (a struct array, possibly
% struct([]) i.e. 0x0 with no fields yet). Avoids MATLAB's "dissimilar
% structures" error when the first element defines the field set.
    if isempty(inArr)
        outArr = entry;
    else
        outArr = inArr;
        outArr(end+1) = entry; %#ok<AGROW>
    end
end

function T = normalize_text_columns(T)
% Forces every char-typed table variable to cellstr, guarding against the
% char-matrix collapse that occurs when every value in a column happens
% to have identical length in a given run (see Phase 0/1 scripts for the
% original discovery of this MATLAB behavior).
    varNames = T.Properties.VariableNames;
    for v = 1:numel(varNames)
        col = T.(varNames{v});
        if ischar(col)
            T.(varNames{v}) = cellstr(col);
        end
    end
end

function bLat = get_boundary_latencies(EEG, boundaryTypes)
    bLat = zeros(0, 1);
    if ~isfield(EEG, 'event') || isempty(EEG.event) || ~isfield(EEG.event, 'type')
        return;
    end
    hasDur = isfield(EEG.event, 'duration');
    for k = 1:numel(EEG.event)
        t = EEG.event(k).type;
        if isnumeric(t) || islogical(t); t = num2str(t); end
        if isstring(t); t = char(t); end
        if ~ischar(t) || ~any(strcmpi(strtrim(t), boundaryTypes))
            continue;
        end
        if hasDur
            d = EEG.event(k).duration;
            if ~isempty(d) && isnumeric(d) && isscalar(d) && ~isnan(d) && d == 0
                continue;   % explicitly zero-length -> not a real splice
            end
        end
        L = double(EEG.event(k).latency);
        if L > 1 && L < EEG.pnts
            bLat(end+1, 1) = L; %#ok<AGROW>
        end
    end
end

function s = summarize_event_types(EEG)
    s = '(no events)';
    if ~isfield(EEG, 'event') || isempty(EEG.event) || ~isfield(EEG.event, 'type')
        return;
    end
    types = cell(numel(EEG.event), 1);
    for k = 1:numel(EEG.event)
        t = EEG.event(k).type;
        if isnumeric(t) || islogical(t); t = num2str(t); end
        if isstring(t); t = char(t); end
        if ~ischar(t); t = 'UNKNOWN_TYPE'; end
        types{k} = t;
    end
    [u, ~, ic] = unique(types);
    cnt = accumarray(ic, 1);
    parts = arrayfun(@(k) sprintf('%s:%d', u{k}, cnt(k)), 1:numel(u), 'UniformOutput', false);
    s = strjoin(parts, '; ');
end
