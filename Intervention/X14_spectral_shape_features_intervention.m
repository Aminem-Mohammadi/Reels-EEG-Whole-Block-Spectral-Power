%% ========================================================================
% X14_spectral_shape_features_intervention.m
%
% Spectral-shape measures: individual alpha frequency (IAF, peak frequency within 8-13 Hz
% of each channel's mean spectrum) and the 95% spectral edge frequency (SEF95, 1-45 Hz).
% Global, regional and electrode level.
% Inputs : 03_psd, 07_roi_electrode_definitions.
% Outputs: 14_spectral_shape_features/<RUN_MODE>/, 09_reports/.
% Run twice, with RUN_MODE = 'full_sample' (N = 30) and RUN_MODE = 'sensitivity_excl_noisy'
% (cleaned sample, N = 28; excluded subjects listed in EXCLUSION_CONFIG_intervention.m).
% Outputs are written to <OUT_DIR>/<folder>/<RUN_MODE>/ and reports to <OUT_DIR>/09_reports/.
% ========================================================================

clear; clc;

%% Section 0: Paths and existence checks
config;
DATA_ROOT_FOR_SIBLING = DATA_DIR;
[DATA_ROOT_PARENT, ~] = fileparts(DATA_ROOT_FOR_SIBLING);
ANALYSIS_ROOT = OUT_DIR;

THIS_SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(THIS_SCRIPT_DIR);
run(fullfile(THIS_SCRIPT_DIR, 'EXCLUSION_CONFIG_intervention.m'));
RUN_MODE = 'full_sample';   % 'full_sample' (N = 30) or 'sensitivity_excl_noisy' (N = 28)
if ~ismember(RUN_MODE, RUN_MODES_intervention)
    error('RunMode:Bad', 'RUN_MODE must be one of: %s', strjoin(RUN_MODES_intervention, ', '));
end
RUN_MODE_LABEL = RUN_MODES_DISPLAY_intervention.(RUN_MODE);
RUN_MODE_EXCL_SUBJ = [];   % filled in after the data are loaded
fprintf('RUN_MODE = %s\n', RUN_MODE_LABEL);
if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
    fprintf('Excluding subject(s) {%s} -- see EXCLUSION_CONFIG_intervention.m for the reasons.\n', ...
        strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','));
end

PSD_DIR    = fullfile(ANALYSIS_ROOT, '03_psd');
DEF_DIR    = fullfile(ANALYSIS_ROOT, '07_roi_electrode_definitions');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');

psdFile    = fullfile(PSD_DIR, 'intervention_psd.mat');
roiDefFile = fullfile(DEF_DIR, 'phase5_channel_roi_definitions.mat');

if ~exist(psdFile, 'file')
    error('Phase9X14:MissingInput', ['X03 PSD file not found: %s\nRun ' ...
        'X03_compute_intervention_psd.m first.'], psdFile);
end
if ~exist(roiDefFile, 'file')
    error('Phase9X14:MissingInput', ['X07 ROI/electrode definition file not found: %s\nRun ' ...
        'X07_roi_electrode_definitions.m first. This script will NOT ' ...
        'reconstruct or guess ROI/electrode membership on its own.'], roiDefFile);
end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end
if isempty(ver('stats'))
    error('Phase9X14:MissingToolbox', ...
        'Statistics and Machine Learning Toolbox not detected. This script requires fitlme, anova(LinearMixedModel), coefTest, and qqplot.');
end

fprintf('========================================================\n');
fprintf('PHASE 9 (X14) -- SPECTRAL-SHAPE FEATURES (IAF, SEF95), INTERVENTION DATA\n');
fprintf('Single intervention analysis (conditions E / N / R).\n');
fprintf('NEW, PURELY ADDITIVE exploratory extension -- never reads or modifies X05/X06/X08/X09/X10 outputs.\n');
fprintf('========================================================\n\n');

%% Section 1: Load X03 (PSD) and X07 (ROI/electrode definitions) -- both
S3 = intervention_load(psdFile, 'psd_data', 'FrequencyVectorHz', 'ChannelLabelsReference', 'PARAMS');
psd_data_ALL           = S3.psd_data;
FrequencyVectorHz      = S3.FrequencyVectorHz(:);
ChannelLabelsReference = S3.ChannelLabelsReference(:);
PSD_PARAMS             = S3.PARAMS;
clear S3;

S7 = load(roiDefFile, 'ChannelDefT', 'ROI_DEFINITION', 'ROI_NAMES', 'UNASSIGNED_FOR_ROI', 'PARAMS5D');
ChannelDefT        = S7.ChannelDefT;
ROI_DEFINITION     = S7.ROI_DEFINITION;
ROI_NAMES          = S7.ROI_NAMES(:);
UNASSIGNED_FOR_ROI = S7.UNASSIGNED_FOR_ROI;
PARAMS5D            = S7.PARAMS5D; %#ok<NASGU>
clear S7;

nChan = numel(ChannelLabelsReference);
nROI  = numel(ROI_NAMES);
if nChan ~= 62
    error('Phase9X14:UnexpectedChannelCount', ...
        'ChannelLabelsReference has %d entries; this project is defined for a 62-channel montage. Stopping.', nChan);
end
if ~isequal(upper(ChannelDefT.Label), upper(ChannelLabelsReference))
    error('Phase9X14:ReferenceOrderMismatch', ...
        'X03''s ChannelLabelsReference and X07''s ChannelDefT.Label are not identical. Stopping -- investigate before proceeding.');
end

% Spot-check every psd_data entry uses the SAME frequency grid and channel
% order as the shared references loaded above -- a structural invariant,
% never assumed.
chanMismatch = false(numel(psd_data_ALL), 1);
freqMismatch = false(numel(psd_data_ALL), 1);
for i = 1:numel(psd_data_ALL)
    chanMismatch(i) = ~isequal(upper(psd_data_ALL(i).ChannelLabels(:)), upper(ChannelLabelsReference));
    freqMismatch(i) = numel(psd_data_ALL(i).FrequenciesHz) ~= numel(FrequencyVectorHz) || ...
        max(abs(psd_data_ALL(i).FrequenciesHz(:) - FrequencyVectorHz)) > 1e-9;
end
if any(chanMismatch)
    error('Phase9X14:ChannelOrderMismatch', ...
        '%d psd_data entries have a ChannelLabels order different from the shared ChannelLabelsReference. Stopping.', sum(chanMismatch));
end
if any(freqMismatch)
    error('Phase9X14:FrequencyGridMismatch', ...
        '%d psd_data entries have a FrequenciesHz grid different from the shared FrequencyVectorHz. Stopping.', sum(freqMismatch));
end
fprintf('Loaded %d psd_data recordings; channel order and frequency grid cross-checked against X03/X07 shared references (%d channels, %d frequency bins).\n\n', ...
    numel(psd_data_ALL), nChan, numel(FrequencyVectorHz));

%% Section 2: Frequency-range setup (Task-invariant, data-driven, done once)
%  IAF sub-range: 8-13 Hz inclusive.
IAF_RANGE_HZ = [8 13];
iafMask = FrequencyVectorHz >= IAF_RANGE_HZ(1) & FrequencyVectorHz <= IAF_RANGE_HZ(2);
if ~any(iafMask)
    error('Phase9X14:NoAlphaBins', 'No FrequencyVectorHz bins fall inside the 8-13 Hz IAF sub-range. Stopping.');
end

%  Low-frequency (1-45 Hz) contiguous retained block: located data-drivenly
%  as the first run of FrequencyVectorHz before the large gap to the
%  55-100 Hz portion, using the SAME gap-detection logic X03's own
%  Diagnostic 1 uses (diff > 2x the nominal frequency resolution) --
%  never hard-coded to a literal [1 45].
if isfield(PSD_PARAMS, 'FrequencyResolutionHz') && ~isempty(PSD_PARAMS.FrequencyResolutionHz)
    freqResHz = PSD_PARAMS.FrequencyResolutionHz;
else
    freqResHz = median(diff(FrequencyVectorHz));
end
freqDiffs = diff(FrequencyVectorHz);
gapPos = find(freqDiffs > 2 * freqResHz, 1, 'first');
if isempty(gapPos)
    error('Phase9X14:NoGapFound', ...
        'No gap (>2x frequency resolution) found in FrequencyVectorHz -- expected a notch-guard gap between the 1-45 Hz and 55-100 Hz retained ranges. Stopping rather than guessing the 1-45 Hz sub-range.');
end
lowFreqMask = false(size(FrequencyVectorHz));
lowFreqMask(1:gapPos) = true;
LowFreqBoundaryHz  = [FrequencyVectorHz(1), FrequencyVectorHz(gapPos)];
HighFreqBoundaryHz = [FrequencyVectorHz(gapPos + 1), FrequencyVectorHz(end)];

if ~all(lowFreqMask(iafMask))   % iafMask must be a subset of lowFreqMask
    error('Phase9X14:AlphaOutsideLowBlock', ...
        'The 8-13 Hz IAF sub-range is not entirely contained in the detected 1-45 Hz low-frequency block. Stopping.');
end

fprintf('Section 2: frequency ranges (data-driven from FrequencyVectorHz, resolution=%.4f Hz):\n', freqResHz);
fprintf('  IAF sub-range               : %g-%g Hz (%d bins)\n', IAF_RANGE_HZ(1), IAF_RANGE_HZ(2), sum(iafMask));
fprintf('  SEF95 low-frequency block   : %.2f-%.2f Hz (%d bins, contiguous, gap detected after this range)\n', ...
    LowFreqBoundaryHz(1), LowFreqBoundaryHz(2), sum(lowFreqMask));
fprintf('  Excluded high-frequency block (never used for SEF95): %.2f-%.2f Hz (%d bins)\n\n', ...
    HighFreqBoundaryHz(1), HighFreqBoundaryHz(2), sum(~lowFreqMask));

%% Section 3: MedianPSDBased variant auto-detection (Task-invariant, done once)
hasFullWindowField = isfield(psd_data_ALL, 'PSD_WindowLevel');
hasNonEmptyFullWindowData = false;
if hasFullWindowField
    for i = 1:numel(psd_data_ALL)
        if ~isempty(psd_data_ALL(i).PSD_WindowLevel)
            hasNonEmptyFullWindowData = true;
            break;
        end
    end
end
MEDIAN_VARIANT_AVAILABLE = hasFullWindowField && hasNonEmptyFullWindowData;

if isfield(PSD_PARAMS, 'SaveFullWindowLevelPSD')
    saveFlagValue = PSD_PARAMS.SaveFullWindowLevelPSD;
else
    saveFlagValue = NaN;
end

if MEDIAN_VARIANT_AVAILABLE
    VARIANTS = {'MeanPSDBased', 'MedianPSDBased'};
    variantSkipMessage = 'MedianPSDBased variant: AVAILABLE and RUN -- psd_data(i).PSD_WindowLevel (per-window, full-resolution PSD) was found in the loaded X03 output.';
    fprintf('Section 3: %s\n\n', variantSkipMessage);
else
    VARIANTS = {'MeanPSDBased'};
    variantSkipMessage = sprintf(['MedianPSDBased variant unavailable -- X03 was run with SaveFullWindowLevelPSD=%s; ' ...
        're-run X03 with that flag set to true to unlock full epoch-level PSD data for this variant. ' ...
        'Proceeding with MeanPSDBased only (this is a clean, expected, disclosed skip, not an error).'], ...
        local_yn_flag(saveFlagValue));
    fprintf('Section 3: %s\n\n', variantSkipMessage);
end

METRICS = {'IAF', 'SEF95'};
CONDITION_ORDER = {'E', 'N', 'R'};
REF_LEVEL = CONDITION_ORDER{1};
CONTRASTS = {'E','N','E - N'; 'E','R','E - R'; 'N','R','N - R'};   % {level2, level1, label}; estimate = level2 - level1
FDR_ALPHA = 0.05;

TASK_TYPES = {'intervention'};
TASK_DISPLAY = struct('intervention', 'Intervention');

%% ========================================================================
for taskI = 1:numel(TASK_TYPES)
    TASK_KEY = TASK_TYPES{taskI};
    TASK_LABEL = TASK_DISPLAY.(TASK_KEY);

    STATS_DIR = fullfile(ANALYSIS_ROOT, '14_spectral_shape_features', RUN_MODE);
    if ~exist(STATS_DIR, 'dir'); mkdir(STATS_DIR); end

    fprintf('\n########################################################\n');
    fprintf('### TASK = %s | RUN_MODE = %s\n', TASK_LABEL, RUN_MODE_LABEL);
    fprintf('########################################################\n\n');

    integrityIssues = {};

    %% Section 4a: filter to this Task, then to Condition in {E,N,R}
    taskMaskIdx = arrayfun(@(x) strcmp(x.TaskType, TASK_KEY), psd_data_ALL);
    thisTaskEntries = psd_data_ALL(taskMaskIdx);
    nBaselineExcluded = sum(arrayfun(@(x) strcmp(x.Condition, 'Baseline'), thisTaskEntries));
    condKeepMask = arrayfun(@(x) ismember(x.Condition, CONDITION_ORDER), thisTaskEntries);
    entries = thisTaskEntries(condKeepMask);
    [entries, RUN_MODE_N_EXCL_ROWS, RUN_MODE_EXCL_SUBJ] = apply_run_mode_filter(entries, RUN_MODE, EXCLUDE_NOISY_SUBJECTS);
    if RUN_MODE_N_EXCL_ROWS > 0; fprintf('RUN_MODE=%s: removed %d recording(s) of subject(s) {%s}.\n', RUN_MODE_LABEL, RUN_MODE_N_EXCL_ROWS, num2str(RUN_MODE_EXCL_SUBJ)); end
    fprintf('Task=%s: %d total recordings loaded, %d Baseline recordings excluded, %d E/N/R recordings retained.\n', ...
        TASK_LABEL, numel(thisTaskEntries), nBaselineExcluded, numel(entries));

    if isempty(entries)
        warning('Task=%s: no E/N/R recordings found after filtering -- skipping this Task entirely.', TASK_LABEL);
        continue;
    end

    %% Section 4b: data-driven design check (same check as X05 Section 2)
    entrySubj  = {entries.Subject}';
    entryBlock = [entries.BlockNumber]';
    entryCond  = {entries.Condition}';
    comboT = unique(table(entrySubj, entryBlock, entryCond, 'VariableNames', {'Subject','BlockNumber','Condition'}), 'rows');
    nRecordings = height(comboT);

    subjNumeric = local_subject_numeric(comboT.Subject);
    uniqSubj = unique(subjNumeric);
    nSubjects = numel(uniqSubj);

    subjPairs = cell(nSubjects, 1);
    pairTypeCount = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for i = 1:nSubjects
        rowsI = subjNumeric == uniqSubj(i);
        condsI = cellstr(comboT.Condition(rowsI));
        if numel(condsI) ~= 2
            integrityIssues{end+1} = sprintf('Subject %g has %d recordings (expected 2) for Task=%s.', uniqSubj(i), numel(condsI), TASK_LABEL); %#ok<AGROW>
        elseif strcmp(condsI{1}, condsI{2})
            integrityIssues{end+1} = sprintf('Subject %g has the SAME condition (%s) in both recordings for Task=%s -- not a valid pair.', uniqSubj(i), condsI{1}, TASK_LABEL); %#ok<AGROW>
        else
            pairKey = strjoin(sort(condsI), '_');
            if isKey(pairTypeCount, pairKey)
                pairTypeCount(pairKey) = pairTypeCount(pairKey) + 1;
            else
                pairTypeCount(pairKey) = 1;
            end
        end
        subjPairs{i} = condsI;
    end
    fprintf('Condition-pair coverage across %d subjects (Task=%s):\n', nSubjects, TASK_LABEL);
    pairKeys = keys(pairTypeCount);
    for k = 1:numel(pairKeys)
        fprintf('  %s : %d subjects\n', strrep(pairKeys{k}, '_', ' & '), pairTypeCount(pairKeys{k}));
    end

    allSubjectsComplete3 = false;
    condsPerSubjCounts = cellfun(@numel, subjPairs);
    if all(condsPerSubjCounts >= 3)
        allSubjectsComplete3 = true;
        integrityIssues{end+1} = sprintf(['Task=%s: every subject appears to have 3+ recordings/conditions -- the incomplete-design ' ...
            'justification for using an LMM instead of RM-ANOVA no longer applies. STOP and re-examine.'], TASK_LABEL);
    end
    if ~isempty(integrityIssues)
        fprintf(2, '\n*** DESIGN VERIFICATION FLAGGED %d ISSUE(S) (Task=%s) -- see report. ***\n\n', numel(integrityIssues), TASK_LABEL);
        for i = 1:numel(integrityIssues)
            fprintf(2, '  - %s\n', integrityIssues{i});
        end
    end
    if allSubjectsComplete3
        error('Phase9X14:DesignAssumptionViolated', ...
            'Task=%s: design verification found a complete 3-condition-per-subject structure. Stopping before fitting any model -- see integrityIssues.', TASK_LABEL);
    end

    %% Section 4c: build the Electrode-level feature-value long table
    %  (Subject, BlockNumber, Condition, Unit=Channel, Metric, Variant, Value)
    fprintf('\nSection 4c: extracting IAF/SEF95 per recording x channel x variant (Task=%s)...\n', TASK_LABEL);
    elecRows = cell(0, 7);
    nSkippedAllNaNCurves = 0;
    for e = 1:numel(entries)
        entry = entries(e);
        for vv = 1:numel(VARIANTS)
            variantName = VARIANTS{vv};
            if strcmp(variantName, 'MeanPSDBased')
                curveMat = entry.MeanPSD;   % nChan x nFreqSel
            else
                curveMat = local_median_psd_curve(entry);   % nChan x nFreqSel
            end

            % IAF: max within the 8-13 Hz sub-range, per channel (vectorized).
            pxxAlphaMat = curveMat(:, iafMask);
            freqAlphaVec = FrequencyVectorHz(iafMask);
            [maxVals, maxIdx] = max(pxxAlphaMat, [], 2);
            iafVec = freqAlphaVec(maxIdx);
            iafVec(any(isnan(pxxAlphaMat), 2) | isnan(maxVals)) = NaN;

            % SEF95: per channel, cumulative-power crossing with linear
            % interpolation, over the 1-45 Hz contiguous block only.
            freqLowVec = FrequencyVectorHz(lowFreqMask);
            sef95Vec = nan(nChan, 1);
            for ch = 1:nChan
                pxxLowRow = curveMat(ch, lowFreqMask)';
                sef95Vec(ch) = local_sef95_from_curve(freqLowVec, pxxLowRow);
            end

            if all(isnan(iafVec)) && all(isnan(sef95Vec))
                nSkippedAllNaNCurves = nSkippedAllNaNCurves + 1;
            end

            for ch = 1:nChan
                elecRows(end+1, :) = {entry.Subject, entry.BlockNumber, entry.Condition, ...
                    ChannelLabelsReference{ch}, 'IAF', variantName, iafVec(ch)}; %#ok<AGROW>
                elecRows(end+1, :) = {entry.Subject, entry.BlockNumber, entry.Condition, ...
                    ChannelLabelsReference{ch}, 'SEF95', variantName, sef95Vec(ch)}; %#ok<AGROW>
            end
        end
    end
    if nSkippedAllNaNCurves > 0
        integrityIssues{end+1} = sprintf(['%d recording x Variant combination(s) had an all-NaN PSD curve (0 valid windows) -- ' ...
            'IAF/SEF95 set to NaN for every channel there, never imputed; fitlme excludes these rows list-wise (Task=%s).'], ...
            nSkippedAllNaNCurves, TASK_LABEL); %#ok<AGROW>
    end

    ElectrodeFeatureT = local_build_table(elecRows, [1 3 4 5 6], [2 7], ...
        {'Subject','BlockNumber','Condition','Unit','Metric','Variant','Value'});
    fprintf('  Electrode-level feature table built: %d rows (%d recordings x %d channels x %d metrics x %d variant(s)).\n', ...
        height(ElectrodeFeatureT), numel(entries), nChan, numel(METRICS), numel(VARIANTS));

    ElectrodeFeatureT.ChannelUC = upper(ElectrodeFeatureT.Unit);

    %% Section 4d: ROI-level aggregation (arithmetic mean across ROI member channels)
    RoiFeatureT = table();
    for r = 1:nROI
        roiName = ROI_NAMES{r};
        roiChansUC = upper(ROI_DEFINITION.(roiName));
        subT = ElectrodeFeatureT(ismember(ElectrodeFeatureT.ChannelUC, roiChansUC), :);
        if isempty(subT); continue; end
        [G, keyT] = findgroups(subT(:, {'Subject','BlockNumber','Condition','Metric','Variant'}));
        meanVals = splitapply(@(v) mean(v, 'omitnan'), subT.Value, G);
        nMembers = splitapply(@numel, subT.Value, G);
        expectedN = numel(roiChansUC);
        if any(nMembers ~= expectedN)
            integrityIssues{end+1} = sprintf('ROI=%s (Task=%s): at least one Subject x Block x Condition x Metric x Variant group has %d member channel(s), expected %d.', ...
                roiName, TASK_LABEL, min(nMembers), expectedN); %#ok<AGROW>
        end
        keyT.Unit = repmat({roiName}, height(keyT), 1);
        keyT.Value = meanVals;
        RoiFeatureT = [RoiFeatureT; keyT]; %#ok<AGROW>
    end
    fprintf('  ROI-level feature table built: %d rows (%d recordings x %d ROIs x %d metrics x %d variant(s)).\n', ...
        height(RoiFeatureT), numel(entries), nROI, numel(METRICS), numel(VARIANTS));

    %% Section 4e: Global-level aggregation (arithmetic mean across all 62 channels)
    [G, keyT] = findgroups(ElectrodeFeatureT(:, {'Subject','BlockNumber','Condition','Metric','Variant'}));
    meanVals = splitapply(@(v) mean(v, 'omitnan'), ElectrodeFeatureT.Value, G);
    nMembers = splitapply(@numel, ElectrodeFeatureT.Value, G);
    if any(nMembers ~= nChan)
        integrityIssues{end+1} = sprintf('Global level (Task=%s): at least one Subject x Block x Condition x Metric x Variant group has %d member channel(s), expected %d.', ...
            TASK_LABEL, min(nMembers), nChan); %#ok<AGROW>
    end
    keyT.Unit = repmat({'Global'}, height(keyT), 1);
    keyT.Value = meanVals;
    GlobalFeatureT = keyT;
    fprintf('  Global-level feature table built: %d rows (%d recordings x %d metrics x %d variant(s)).\n\n', ...
        height(GlobalFeatureT), numel(entries), numel(METRICS), numel(VARIANTS));

    %% Section 4f: fit LMMs at all three levels, all Metrics, all Variants
    fprintf('--------------------------------------------------------\n');
    fprintf('SECTION 4f: MODEL FITTING (Task=%s, Electrode/ROI/Global x Metric x Variant)\n', TASK_LABEL);
    fprintf('--------------------------------------------------------\n');

    LEVELS = {'Electrode', 'ROI', 'Global'};
    LevelTables = containers.Map(LEVELS, {ElectrodeFeatureT, RoiFeatureT, GlobalFeatureT});
    LevelUnits  = containers.Map(LEVELS, {ChannelLabelsReference, ROI_NAMES, {'Global'}});

    omnibusRows = cell(0, 8);   % Level,Unit,Metric,Variant,FStat,DF1,DF2,pValue
    lmeMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    convergenceWarnings = {};
    nSkippedTooFewObs = 0;

    for lv = 1:numel(LEVELS)
        levelName = LEVELS{lv};
        LT = LevelTables(levelName);
        units = LevelUnits(levelName);
        for m = 1:numel(METRICS)
            metricName = METRICS{m};
            for vv = 1:numel(VARIANTS)
                variantName = VARIANTS{vv};
                for u = 1:numel(units)
                    unitName = units{u};
                    subT = LT(strcmp(LT.Metric, metricName) & strcmp(LT.Variant, variantName) & strcmp(LT.Unit, unitName), :);
                    subT = subT(~isnan(subT.Value), :);
                    if height(subT) < 4
                        nSkippedTooFewObs = nSkippedTooFewObs + 1;
                        continue;
                    end
                    subT.Subject = categorical(local_subject_numeric(subT.Subject));
                    subT.Condition = categorical(cellstr(subT.Condition), CONDITION_ORDER);
                    formula = 'Value ~ Condition + (1|Subject)';
                    try
                        lastwarn('');
                        lme = fitlme(subT, formula, 'FitMethod', 'REML');
                        [warnMsg, ~] = lastwarn();
                        if ~isempty(warnMsg)
                            convergenceWarnings{end+1} = sprintf('%s / %s / %s / %s: %s', levelName, unitName, metricName, variantName, warnMsg); %#ok<AGROW>
                        end
                    catch ME
                        integrityIssues{end+1} = sprintf('fitlme FAILED for %s / %s / %s / %s (Task=%s): %s', ...
                            levelName, unitName, metricName, variantName, TASK_LABEL, ME.message); %#ok<AGROW>
                        continue;
                    end
                    aov = anova(lme, 'DFMethod', 'satterthwaite');
                    rowIdx = find(strcmp(aov.Term, 'Condition'), 1);
                    if isempty(rowIdx)
                        integrityIssues{end+1} = sprintf('anova(lme) has no Condition row for %s / %s / %s / %s (Task=%s).', ...
                            levelName, unitName, metricName, variantName, TASK_LABEL); %#ok<AGROW>
                        continue;
                    end
                    key = sprintf('%s|%s|%s|%s', levelName, unitName, metricName, variantName);
                    lmeMap(key) = lme;
                    omnibusRows(end+1, :) = {levelName, unitName, metricName, variantName, ...
                        aov.FStat(rowIdx), aov.DF1(rowIdx), aov.DF2(rowIdx), aov.pValue(rowIdx)}; %#ok<AGROW>
                end
            end
        end
    end
    if nSkippedTooFewObs > 0
        integrityIssues{end+1} = sprintf('%d (Level x Unit x Metric x Variant) combination(s) had fewer than 4 non-missing observations and were skipped (Task=%s).', nSkippedTooFewObs, TASK_LABEL); %#ok<AGROW>
    end

    OmnibusT = local_build_table(omnibusRows, [1 2 3 4], [5 6 7 8], ...
        {'Level','Unit','Metric','Variant','FStat','DF1','DF2','pValue'});
    fprintf('Model fitting complete (Task=%s): %d omnibus tests fit successfully.\n', TASK_LABEL, height(OmnibusT));

    %% Section 4g: BH-FDR within each (Level x Metric x Variant) family
    OmnibusT.pFDR = nan(height(OmnibusT), 1);
    OmnibusT.SignificantFDR = false(height(OmnibusT), 1);
    OmnibusT.FamilySize = nan(height(OmnibusT), 1);

    for lv = 1:numel(LEVELS)
        levelName = LEVELS{lv};
        for m = 1:numel(METRICS)
            metricName = METRICS{m};
            for vv = 1:numel(VARIANTS)
                variantName = VARIANTS{vv};
                idx = strcmp(OmnibusT.Level, levelName) & strcmp(OmnibusT.Metric, metricName) & strcmp(OmnibusT.Variant, variantName);
                nInFamily = sum(idx);
                if nInFamily == 0; continue; end
                OmnibusT.FamilySize(idx) = nInFamily;
                pvals = OmnibusT.pValue(idx);
                if nInFamily == 1
                    OmnibusT.pFDR(idx) = pvals;
                    OmnibusT.SignificantFDR(idx) = pvals <= FDR_ALPHA;
                else
                    [adjP, reject] = local_fdr_bh(pvals, FDR_ALPHA);
                    OmnibusT.pFDR(idx) = adjP;
                    OmnibusT.SignificantFDR(idx) = reject;
                end
            end
        end
    end
    fprintf('\n--- Spectral-shape omnibus results (Task=%s) ---\n', TASK_LABEL);
    disp(OmnibusT);

    %% Section 4h: pairwise contrasts (E-N, E-R, N-R), gated on FDR, Holm-corrected
    pairwiseRows = cell(0, 13);   % Level,Unit,Metric,Variant,Contrast,Estimate,SE,DF,tStat,pValue,pHolm,CI_Low,CI_High
    sigMask = OmnibusT.SignificantFDR;
    for row = find(sigMask)'
        levelName = OmnibusT.Level{row};
        unitName = OmnibusT.Unit{row};
        metricName = OmnibusT.Metric{row};
        variantName = OmnibusT.Variant{row};
        key = sprintf('%s|%s|%s|%s', levelName, unitName, metricName, variantName);
        if ~isKey(lmeMap, key); continue; end
        lme = lmeMap(key);
        coeffNames = lme.CoefficientNames;
        beta = lme.Coefficients.Estimate;
        covB = lme.CoefficientCovariance;

        contrastP = nan(3, 1);
        contrastData = cell(3, 8);
        for c = 1:3
            lvl2 = CONTRASTS{c, 1}; lvl1 = CONTRASTS{c, 2}; label = CONTRASTS{c, 3};
            C = local_condition_contrast(coeffNames, 'Condition', REF_LEVEL, lvl1, lvl2);
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
        [adjP, ~] = local_holm_bonferroni(contrastP, FDR_ALPHA);
        for c = 1:3
            pairwiseRows(end+1, :) = [{levelName, unitName, metricName, variantName}, contrastData(c, 1:6), {adjP(c)}, contrastData(c, 7:8)]; %#ok<AGROW>
        end
    end
    PairwiseT = local_build_table(pairwiseRows, [1 2 3 4 5], [6 7 8 9 10 11 12 13], ...
        {'Level','Unit','Metric','Variant','Contrast','Estimate','SE','DF','tStat','pValue','pHolm','CI_Low','CI_High'});
    fprintf('\n--- Pairwise contrasts (Task=%s; Level x Unit x Metric x Variant significant under FDR) ---\n', TASK_LABEL);
    disp(PairwiseT);

    %% Section 4i: residual diagnostics -- Global level ONLY
    fprintf('\n--------------------------------------------------------\n');
    fprintf('SECTION 4i: MODEL DIAGNOSTICS (Task=%s, Global level only)\n', TASK_LABEL);
    fprintf('--------------------------------------------------------\n');
    nDiagSaved = 0;
    for m = 1:numel(METRICS)
        metricName = METRICS{m};
        for vv = 1:numel(VARIANTS)
            variantName = VARIANTS{vv};
            key = sprintf('Global|Global|%s|%s', metricName, variantName);
            if ~isKey(lmeMap, key); continue; end
            lme = lmeMap(key);
            try
                res = residuals(lme, 'ResidualType', 'Pearson');
                fig = figure('Visible', 'off', 'Position', [100 100 900 400]);
                subplot(1,2,1);
                histogram(res, 15);
                title(sprintf('Residuals: Global %s (%s), Task=%s', metricName, variantName, TASK_LABEL), 'Interpreter', 'none');
                xlabel('Pearson residual'); ylabel('Count');
                subplot(1,2,2);
                qqplot(res);
                title('Normal Q-Q');
                outPng = fullfile(STATS_DIR, sprintf('diagnostic_residuals_Global_%s_%s.png', metricName, variantName));
                saveas(fig, outPng);
                close(fig);
                nDiagSaved = nDiagSaved + 1;
            catch ME
                integrityIssues{end+1} = sprintf('Residual diagnostic failed for Global %s / %s (Task=%s): %s', metricName, variantName, TASK_LABEL, ME.message); %#ok<AGROW>
            end
        end
    end
    fprintf('Saved %d Global-level residual diagnostic plots (Task=%s).\n', nDiagSaved, TASK_LABEL);

    %% Section 4j: PARAMS14 and outputs
    PARAMS14 = struct();
    PARAMS14.ScriptName = 'X14_spectral_shape_features_intervention.m';
    PARAMS14.Task = TASK_LABEL;
    PARAMS14.GeneratedOn = datestr(now); %#ok<TNOW1,DATST>
    PARAMS14.AnalysisLevel = sprintf('Spectral-shape features (IAF, SEF95) from X03 PSD curves, Task=%s only, Condition in {E,N,R} only (Baseline excluded)', TASK_LABEL);
    PARAMS14.Model = 'LinearMixedEffects_RandomIntercept';
    PARAMS14.FixedEffectsPrimary = 'Condition (3-level: E, N, R; reference = E)';
    PARAMS14.RandomEffects = '(1|Subject)';
    PARAMS14.FitMethod = 'REML';
    PARAMS14.DFMethod = 'Satterthwaite';
    PARAMS14.Metrics = METRICS;
    PARAMS14.Variants = VARIANTS;
    PARAMS14.MedianVariantAvailable = MEDIAN_VARIANT_AVAILABLE;
    PARAMS14.MedianVariantNote = variantSkipMessage;
    PARAMS14.IAFRangeHz = IAF_RANGE_HZ;
    PARAMS14.SEF95LowFreqBoundaryHz = LowFreqBoundaryHz;
    PARAMS14.SEF95ExcludedHighFreqBoundaryHz = HighFreqBoundaryHz;
    PARAMS14.SEF95Threshold = 0.95;
    PARAMS14.Levels = LEVELS;
    PARAMS14.ROINames = ROI_NAMES;
    PARAMS14.FDRAlpha = FDR_ALPHA;
    PARAMS14.FDRFamilyScheme = 'Separate BH-FDR family per (Level x Metric x Variant); Electrode family size=62, ROI family size=5, Global family size=1 (reported as-is, no correction). Own families, never pooled with locked X05/X06/X08/X09 metrics or other X11-X15 families.';
    PARAMS14.PairwiseContrastMethod = 'Holm-Bonferroni across the 3 pairwise contrasts (E-N, E-R, N-R -- standardized on the X05/X06 convention), computed only for Level x Unit x Metric x Variant combinations whose FDR family test survives BH-FDR.';
    PARAMS14.ROIGlobalAggregation = 'Arithmetic mean of per-channel IAF/SEF95 (frequency values) across ROI member channels (ROI level) or all 62 reference channels (Global level); no transform-order ambiguity since these are frequency, not power, values.';
    PARAMS14.BaselineExclusionNote = sprintf('%d Baseline recording(s) excluded before any modeling; Baseline is processed by X03 but never enters this or any other inferential model in this pipeline.', nBaselineExcluded);
    PARAMS14.SourcePsdFile = psdFile;
    PARAMS14.SourceRoiDefFile = roiDefFile;
    PARAMS14.DesignVerification = struct('nRecordings', nRecordings, 'nSubjects', nSubjects, ...
        'PairTypeCounts', {pairTypeCount.keys}, 'IntegrityIssues', {integrityIssues});
    PARAMS14.RunDateTime = datestr(now); %#ok<TNOW1,DATST>

    save(fullfile(STATS_DIR, 'phase14_parameters.mat'), 'PARAMS14');

    FeatureValueT = [ElectrodeFeatureT(:, {'Unit','Subject','BlockNumber','Condition','Metric','Variant','Value'}); ...
                     RoiFeatureT(:, {'Unit','Subject','BlockNumber','Condition','Metric','Variant','Value'}); ...
                     GlobalFeatureT(:, {'Unit','Subject','BlockNumber','Condition','Metric','Variant','Value'})];
    levelTagElec = repmat({'Electrode'}, height(ElectrodeFeatureT), 1);
    levelTagRoi  = repmat({'ROI'}, height(RoiFeatureT), 1);
    levelTagGlob = repmat({'Global'}, height(GlobalFeatureT), 1);
    FeatureValueT.Level = [levelTagElec; levelTagRoi; levelTagGlob];
    FeatureValueT = FeatureValueT(:, {'Level','Unit','Subject','BlockNumber','Condition','Metric','Variant','Value'});

    save(fullfile(STATS_DIR, 'phase14_feature_values.mat'), 'FeatureValueT');
    writetable(FeatureValueT, fullfile(STATS_DIR, 'phase14_feature_values.xlsx'));
    writetable(FeatureValueT, fullfile(STATS_DIR, 'phase14_feature_values.csv'));

    save(fullfile(STATS_DIR, 'phase14_omnibus_results.mat'), 'OmnibusT');
    writetable(OmnibusT, fullfile(STATS_DIR, 'phase14_omnibus_results.xlsx'));
    writetable(OmnibusT, fullfile(STATS_DIR, 'phase14_omnibus_results.csv'));

    save(fullfile(STATS_DIR, 'phase14_pairwise_contrasts.mat'), 'PairwiseT');
    if ~isempty(PairwiseT)
        writetable(PairwiseT, fullfile(STATS_DIR, 'phase14_pairwise_contrasts.xlsx'));
        writetable(PairwiseT, fullfile(STATS_DIR, 'phase14_pairwise_contrasts.csv'));
    end

    fprintf('\nSaved output to: %s\n', STATS_DIR);

    %% Section 4k: Report
    reportLines = {};
    reportLines{end+1} = repmat('=', 1, 60);
    reportLines{end+1} = sprintf('PHASE 9 (X14) -- SPECTRAL-SHAPE FEATURES REPORT (Task=%s)', TASK_LABEL);
    reportLines{end+1} = repmat('=', 1, 60);
    reportLines{end+1} = '';
    reportLines{end+1} = 'Exploratory extension (X11-X15 family).';
    reportLines{end+1} = 'Never reads or modifies X05/X06/X08/X09/X10 outputs.';
    reportLines{end+1} = '';
    reportLines{end+1} = sprintf('Baseline recordings excluded before modeling: %d', nBaselineExcluded);
    reportLines{end+1} = sprintf('Recordings verified: %d | Subjects: %d', nRecordings, nSubjects);
    reportLines{end+1} = sprintf('Condition-pair coverage: %s', strjoin(cellfun(@(k) sprintf('%s=%d', strrep(k,'_',' & '), pairTypeCount(k)), pairKeys, 'UniformOutput', false), ', '));
    reportLines{end+1} = '';
    reportLines{end+1} = 'Model: Linear mixed-effects, Condition (fixed, 3-level) + (1|Subject), REML, Satterthwaite DF.';
    reportLines{end+1} = 'Metrics: IAF (peak alpha frequency, 8-13 Hz) and SEF95 (95% spectral edge frequency, 1-45 Hz only), co-exploratory, parallel, never pooled for FDR.';
    reportLines{end+1} = sprintf('IAF sub-range: %g-%g Hz.', IAF_RANGE_HZ(1), IAF_RANGE_HZ(2));
    reportLines{end+1} = sprintf('SEF95 low-frequency integration block (data-driven from FrequencyVectorHz): %.2f-%.2f Hz.', LowFreqBoundaryHz(1), LowFreqBoundaryHz(2));
    reportLines{end+1} = sprintf('Excluded high-frequency block (never used for SEF95): %.2f-%.2f Hz.', HighFreqBoundaryHz(1), HighFreqBoundaryHz(2));
    reportLines{end+1} = '';
    reportLines{end+1} = '--- MedianPSDBased variant status ---';
    reportLines{end+1} = ['  ' variantSkipMessage];
    reportLines{end+1} = sprintf('  Variants run this Task: %s', strjoin(VARIANTS, ', '));
    reportLines{end+1} = '';
    reportLines{end+1} = 'FDR family scheme: separate BH-FDR per (Level x Metric x Variant); Electrode family size=62, ROI family size=5, Global family size=1 (reported as-is).';
    reportLines{end+1} = 'ROI/Global aggregation: arithmetic mean of per-channel IAF/SEF95 across ROI member channels or all 62 channels (no power-metric transform-order ambiguity -- see script header).';
    reportLines{end+1} = 'Pairwise contrast convention: E-N, E-R, N-R (standardized on X05/X06; new families never inherit X08/X09''s E-N/E-R/N-R convention).';
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Omnibus results (all Level x Metric x Variant x Unit tests) ---';
    for r = 1:height(OmnibusT)
        reportLines{end+1} = sprintf('  [%-10s|%-8s] %-6s (%-14s) F(%g,%.1f)=%.3f p=%.4f pFDR=%.4f family_n=%d %s', ...
            OmnibusT.Level{r}, OmnibusT.Unit{r}, OmnibusT.Metric{r}, OmnibusT.Variant{r}, OmnibusT.DF1(r), OmnibusT.DF2(r), ...
            OmnibusT.FStat(r), OmnibusT.pValue(r), OmnibusT.pFDR(r), OmnibusT.FamilySize(r), ...
            local_yn(OmnibusT.SignificantFDR(r))); %#ok<AGROW>
    end
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Pairwise contrasts (Level x Unit x Metric x Variant significant under FDR; Holm-corrected within unit) ---';
    if isempty(PairwiseT)
        reportLines{end+1} = '  (none -- no Level x Unit x Metric x Variant combination survived its FDR family correction)';
    else
        for r = 1:height(PairwiseT)
            reportLines{end+1} = sprintf('  [%-10s|%-8s] %-6s (%-14s) %-8s est=%.4f (95%% CI %.4f to %.4f) t=%.3f df=%.1f p=%.4f pHolm=%.4f', ...
                PairwiseT.Level{r}, PairwiseT.Unit{r}, PairwiseT.Metric{r}, PairwiseT.Variant{r}, PairwiseT.Contrast{r}, PairwiseT.Estimate(r), ...
                PairwiseT.CI_Low(r), PairwiseT.CI_High(r), PairwiseT.tStat(r), PairwiseT.DF(r), ...
                PairwiseT.pValue(r), PairwiseT.pHolm(r)); %#ok<AGROW>
        end
    end
    reportLines{end+1} = '';
    if ~isempty(integrityIssues)
        reportLines{end+1} = '--- INTEGRITY ISSUES / DISCLOSED SKIPS (reported, not auto-corrected) ---';
        for i = 1:numel(integrityIssues)
            reportLines{end+1} = ['  - ' integrityIssues{i}]; %#ok<AGROW>
        end
    else
        reportLines{end+1} = 'Integrity issues: none.';
    end
    if ~isempty(convergenceWarnings)
        reportLines{end+1} = '';
        reportLines{end+1} = '--- Model convergence warnings ---';
        for i = 1:numel(convergenceWarnings)
            reportLines{end+1} = ['  - ' convergenceWarnings{i}]; %#ok<AGROW>
        end
    end
    reportLines{end+1} = '';
    reportLines{end+1} = 'Output files written to:';
    reportLines{end+1} = ['  ' STATS_DIR];
    reportLines{end+1} = repmat('=', 1, 60);

    reportText = strjoin(reportLines, '\n');
    fprintf('%s\n', reportText);
    fid = fopen(fullfile(REPORT_DIR, sprintf('phase_spectral_shape_features_report_%s_%s.txt', TASK_LABEL, RUN_MODE_LABEL)), 'w');
    fprintf(fid, '%s\n', reportText);
    fclose(fid);

end % TASK_TYPES loop

fprintf('\n=== X14 complete for all three Tasks. Purely additive -- X05/X06/X08/X09/X10 outputs untouched. ===\n');

%% ======================= LOCAL FUNCTIONS (all after top-level code) =======================

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

function [adjP, reject] = local_fdr_bh(pvals, alphaFDR)
    pvals = pvals(:);
    m = numel(pvals);
    [sortedP, sortIdx] = sort(pvals);
    adjSorted = nan(m, 1);
    runningMin = 1;
    for i = m:-1:1
        val = sortedP(i) * m / i;
        runningMin = min(runningMin, val);
        adjSorted(i) = runningMin;
    end
    adjP = nan(m, 1);
    adjP(sortIdx) = adjSorted;
    adjP = min(adjP, 1);
    reject = adjP <= alphaFDR;
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
            error('Phase9X14:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name1);
        end
        cRow(idx1) = cRow(idx1) - 1;
    end
    if ~strcmp(level2, refLevel)
        name2 = sprintf('%s_%s', varName, level2);
        idx2 = find(strcmp(coeffNames, name2));
        if isempty(idx2)
            error('Phase9X14:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name2);
        end
        cRow(idx2) = cRow(idx2) + 1;
    end
end

function s = local_yn(tf)
    if tf
        s = 'SIGNIFICANT';
    else
        s = 'n.s.';
    end
end

function s = local_yn_flag(flagVal)
    if isnan(flagVal)
        s = '<unrecorded in loaded PARAMS>';
    elseif flagVal
        s = 'true';
    else
        s = 'false';
    end
end

function T = local_build_table(rows, textColIdx, numColIdx, varNames)
    nCols = numel(varNames);
    if size(rows, 2) ~= nCols
        error('Phase9X14:BuildTableColumnMismatch', ...
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
            error('Phase9X14:BuildTableColumnMismatch', 'Column %d not classified as text or numeric.', j);
        end
    end
    T = table(cols{:}, 'VariableNames', varNames);
end

function s = local_fmt(v)
    if isnan(v)
        s = 'NaN';
    else
        s = sprintf('%.4f', v);
    end
end

function medCurve = local_median_psd_curve(entry)
    % Computes the MEDIAN-across-valid-windows PSD curve (nChan x nFreqSel)
    % from entry.PSD_WindowLevel (nCandidateWindows x nChan x nFreqSel),
    % using only windows flagged valid by entry.WindowValid. Returns an
    % all-NaN curve (same shape as entry.MeanPSD) if no valid windows or
    % no per-window data are present -- never fabricates a value.
    nChanLocal = size(entry.MeanPSD, 1);
    nFreqLocal = size(entry.MeanPSD, 2);
    if ~isfield(entry, 'PSD_WindowLevel') || isempty(entry.PSD_WindowLevel) || ~any(entry.WindowValid)
        medCurve = nan(nChanLocal, nFreqLocal);
        return;
    end
    validData = entry.PSD_WindowLevel(entry.WindowValid, :, :);   % nValid x nChan x nFreqSel
    medCurve = squeeze(median(validData, 1, 'omitnan'));          % nChan x nFreqSel
end

function sef95 = local_sef95_from_curve(freqLowVec, pxxLowVec)
    % 95% Spectral Edge Frequency via cumulative trapezoidal integration
    % over the supplied (contiguous, low-frequency-only) freq/power pair,
    % with linear interpolation between the two grid points that straddle
    % the 95%-of-total-power threshold.
    if any(isnan(pxxLowVec)) || any(~isfinite(pxxLowVec))
        sef95 = NaN;
        return;
    end
    cumP = cumtrapz(freqLowVec, pxxLowVec);
    total = cumP(end);
    if ~isfinite(total) || total <= 0
        sef95 = NaN;
        return;
    end
    target = 0.95 * total;
    idx = find(cumP >= target, 1, 'first');
    if isempty(idx)
        sef95 = freqLowVec(end);
        return;
    end
    if idx == 1
        sef95 = freqLowVec(1);
        return;
    end
    c0 = cumP(idx - 1); c1 = cumP(idx);
    f0 = freqLowVec(idx - 1); f1 = freqLowVec(idx);
    if c1 == c0
        sef95 = f1;
    else
        frac = (target - c0) / (c1 - c0);
        sef95 = f0 + frac * (f1 - f0);
    end
end
