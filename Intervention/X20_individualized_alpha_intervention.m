%% ========================================================================
% X20_individualized_alpha_intervention.m
%
% Individualized alpha power: power in the band IAF +/- 2 Hz of each channel, with IAF
% taken from X14; global, regional and electrode level.
% Inputs : 03_psd, 14_spectral_shape_features.
% Outputs: 20_individualized_alpha/<RUN_MODE>/, 09_reports/.
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

PSD_DIR    = fullfile(ANALYSIS_ROOT, '03_psd');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');

psdFile = fullfile(PSD_DIR, 'intervention_psd.mat');

if ~exist(psdFile, 'file')
    error('PhaseX20:MissingInput', ['X03 PSD file not found: %s\nRun ' ...
        'X03_compute_intervention_psd.m first.'], psdFile);
end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

addpath(fileparts(mfilename('fullpath')));
run(fullfile(fileparts(mfilename('fullpath')), 'EXCLUSION_CONFIG_intervention.m'));

RUN_MODE = 'full_sample';   % 'full_sample' (N = 30) or 'sensitivity_excl_noisy' (N = 28)
if ~ismember(RUN_MODE, RUN_MODES_intervention)
    error('PhaseX20:BadRunMode', 'RUN_MODE must be one of: %s', strjoin(RUN_MODES_intervention, ', '));
end
RUN_MODE_LABEL = RUN_MODES_DISPLAY_intervention.(RUN_MODE);
fprintf('RUN_MODE = %s\n', RUN_MODE_LABEL);
if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
    fprintf('Excluding subject(s) {%s} from every model in this run -- see EXCLUSION_CONFIG_intervention.m for rationale.\n', ...
        strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','));
end
% ---- END NEW ----

if isempty(ver('stats'))
    error('PhaseX20:MissingToolbox', ...
        'Statistics and Machine Learning Toolbox not detected. This script requires fitlme, anova(LinearMixedModel), coefTest, and qqplot.');
end

fprintf('========================================================\n');
fprintf('X20 -- INDIVIDUALIZED ALPHA (IAF +/- 2 Hz) BAND-POWER FAMILY, INTERVENTION DATA\n');
fprintf('Single intervention analysis (conditions E / N / R).\n');
fprintf('NEW, PURELY ADDITIVE exploratory extension -- never reads or modifies any other script''s outputs.\n');
fprintf('IAF reused read-only from X14''s saved output; IAF peak-finding is NEVER recomputed here.\n');
fprintf('========================================================\n\n');

%% Section 1: Load X03 (PSD) -- shared, Task-invariant, read-only, loaded ONCE
S3 = intervention_load(psdFile, 'psd_data', 'FrequencyVectorHz', 'ChannelLabelsReference', 'PARAMS');
psd_data_ALL           = S3.psd_data;
FrequencyVectorHz      = S3.FrequencyVectorHz(:);
ChannelLabelsReference = S3.ChannelLabelsReference(:);
PSD_PARAMS             = S3.PARAMS;
clear S3;

nChan = numel(ChannelLabelsReference);
if nChan ~= 62
    error('PhaseX20:UnexpectedChannelCount', ...
        'ChannelLabelsReference has %d entries; this project is defined for a 62-channel montage. Stopping.', nChan);
end

% Spot-check every psd_data entry uses the SAME frequency grid and channel
% order as the shared references loaded above -- a structural invariant,
% never assumed (same check X14 already performs on this same file).
chanMismatch = false(numel(psd_data_ALL), 1);
freqMismatch = false(numel(psd_data_ALL), 1);
for i = 1:numel(psd_data_ALL)
    chanMismatch(i) = ~isequal(upper(psd_data_ALL(i).ChannelLabels(:)), upper(ChannelLabelsReference));
    freqMismatch(i) = numel(psd_data_ALL(i).FrequenciesHz) ~= numel(FrequencyVectorHz) || ...
        max(abs(psd_data_ALL(i).FrequenciesHz(:) - FrequencyVectorHz)) > 1e-9;
end
if any(chanMismatch)
    error('PhaseX20:ChannelOrderMismatch', ...
        '%d psd_data entries have a ChannelLabels order different from the shared ChannelLabelsReference. Stopping.', sum(chanMismatch));
end
if any(freqMismatch)
    error('PhaseX20:FrequencyGridMismatch', ...
        '%d psd_data entries have a FrequenciesHz grid different from the shared FrequencyVectorHz. Stopping.', sum(freqMismatch));
end
fprintf('Loaded %d psd_data recordings; channel order and frequency grid cross-checked against X03''s shared references (%d channels, %d frequency bins).\n\n', ...
    numel(psd_data_ALL), nChan, numel(FrequencyVectorHz));

if isfield(PSD_PARAMS, 'RelativePowerDenominatorRangeHz') && ~isempty(PSD_PARAMS.RelativePowerDenominatorRangeHz)
    DENOM_RANGE_HZ = PSD_PARAMS.RelativePowerDenominatorRangeHz;
else
    error('PhaseX20:MissingDenominatorRange', ...
        'PARAMS.RelativePowerDenominatorRangeHz not found in the loaded X03 PARAMS. This script uses that SAME range as the [1,45] Hz clipping boundary and RelativePower denominator range -- refusing to invent a different one.');
end

%% Section 2: Fixed configuration (Task-invariant)
IAF_HALF_WIDTH_HZ = 2;
IAF_VARIANT_USED  = 'MeanPSDBased';   % deliberate pairing choice -- see header
IAF_SEARCH_RANGE_HZ = [8 13];   % X14's own IAF search range (informational only, used for a sanity check below)

METRICS = {'LogAbsolutePower_IndividualizedAlpha', 'RelativePower_IndividualizedAlpha'};
CONDITION_ORDER = {'E', 'N', 'R'};
REF_LEVEL = CONDITION_ORDER{1};
CONTRASTS = {'E','N','E - N'; 'E','R','E - R'; 'N','R','N - R'};   % {level2, level1, label}; estimate = level2 - level1
FDR_ALPHA = 0.05;

TASK_TYPES = {'intervention'};
TASK_DISPLAY = struct('intervention', 'Intervention');

PriorTaskGranularity = '';

fprintf('Section 2: IAF half-width = %g Hz; IAF variant used = %s; RelativePower/clipping denominator range = [%g, %g] Hz.\n\n', ...
    IAF_HALF_WIDTH_HZ, IAF_VARIANT_USED, DENOM_RANGE_HZ(1), DENOM_RANGE_HZ(2));

%% ========================================================================
for taskI = 1:numel(TASK_TYPES)
    TASK_KEY = TASK_TYPES{taskI};
    TASK_LABEL = TASK_DISPLAY.(TASK_KEY);

    STATS_DIR = fullfile(ANALYSIS_ROOT, '20_individualized_alpha', RUN_MODE);
    if ~exist(STATS_DIR, 'dir'); mkdir(STATS_DIR); end

    fprintf('\n########################################################\n');
    fprintf('### TASK = %s | RUN_MODE = %s\n', TASK_LABEL, RUN_MODE_LABEL);
    fprintf('########################################################\n\n');

    integrityIssues = {};
    convergenceWarnings = {};

    %% Section 3: load this Task's X14 IAF output (per-Task file; must
    %  already exist -- X14 is a hard prerequisite, run first)
    iafFile = fullfile(ANALYSIS_ROOT, '14_spectral_shape_features', RUN_MODE, 'phase14_feature_values.mat')
    if ~exist(iafFile, 'file')
        error('PhaseX20:MissingInput', ['X14 IAF output not found for Task=%s: %s\nRun ' ...
            'X14_spectral_shape_features_intervention.m first (all three Tasks) before running X20. ' ...
            'This script REUSES X14''s already-computed IAF values and never recomputes peak-finding on its own.'], ...
            TASK_LABEL, iafFile);
    end
    SIAF = load(iafFile, 'FeatureValueT');
    IAFAllT = SIAF.FeatureValueT;
    clear SIAF;

    requiredIafCols = {'Level','Unit','Subject','BlockNumber','Condition','Metric','Variant','Value'};
    missingIafCols = requiredIafCols(~ismember(requiredIafCols, IAFAllT.Properties.VariableNames));
    if ~isempty(missingIafCols)
        error('PhaseX20:MissingColumns', 'X14''s FeatureValueT (Task=%s) is missing expected column(s): %s', ...
            TASK_LABEL, strjoin(missingIafCols, ', '));
    end

    iafOnlyT = IAFAllT(strcmp(IAFAllT.Metric, 'IAF'), :);
    availableVariants = unique(iafOnlyT.Variant);
    if ~ismember(IAF_VARIANT_USED, availableVariants)
        error('PhaseX20:MissingVariant', ...
            'X14''s IAF output (Task=%s) does not contain the "%s" variant this script requires. Available variant(s): %s.', ...
            TASK_LABEL, IAF_VARIANT_USED, strjoin(availableVariants, ', '));
    end
    iafOnlyT = iafOnlyT(strcmp(iafOnlyT.Variant, IAF_VARIANT_USED), :);

    availableLevels = unique(iafOnlyT.Level);
    hasElectrode = ismember('Electrode', availableLevels);
    hasGlobal    = ismember('Global', availableLevels);
    if hasElectrode
        GRANULARITY = 'PerChannel';
        IAF_LEVEL_KEY = 'Electrode';
        UNITS = ChannelLabelsReference;
    elseif hasGlobal
        GRANULARITY = 'WholeScalp';
        IAF_LEVEL_KEY = 'Global';
        UNITS = {'Global'};
    else
        error('PhaseX20:NoUsableIafLevel', ...
            'X14''s IAF output (Task=%s, Variant=%s) has neither an Electrode nor a Global level. Stopping -- cannot proceed without a usable IAF granularity.', ...
            TASK_LABEL, IAF_VARIANT_USED);
    end
    fprintf('Section 3 (Task=%s): IAF granularity detected from X14''s own output = %s (Level=''%s'' found).\n', TASK_LABEL, GRANULARITY, IAF_LEVEL_KEY);
    if strcmp(GRANULARITY, 'PerChannel')
        fprintf('  -> Confirmed case: one individualized band per recording x channel; FDR family per metric = across channels.\n\n');
    else
        fprintf('  -> DEFENSIVE FALLBACK case (not expected on inspection of X14''s code): one band per recording using the Global IAF; family size = 1, no FDR correction, pairwise gated on nominal p<=0.05 -- see script header.\n\n');
    end
    if ~isempty(PriorTaskGranularity) && ~strcmp(PriorTaskGranularity, GRANULARITY)
        error('PhaseX20:InconsistentGranularityAcrossTasks', ...
            'Task=%s detected IAF granularity "%s", but an earlier Task detected "%s". X14''s output should be structurally identical across Tasks -- stopping rather than silently mixing granularities.', ...
            TASK_LABEL, GRANULARITY, PriorTaskGranularity);
    end
    PriorTaskGranularity = GRANULARITY;

    iafSubT = iafOnlyT(strcmp(iafOnlyT.Level, IAF_LEVEL_KEY), :);
    iafMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
    nDuplicateIafKeys = 0;
    for i = 1:height(iafSubT)
        k = local_iaf_key(iafSubT.Subject{i}, iafSubT.BlockNumber(i), iafSubT.Condition{i}, iafSubT.Unit{i});
        if isKey(iafMap, k)
            nDuplicateIafKeys = nDuplicateIafKeys + 1;
            continue;
        end
        iafMap(k) = iafSubT.Value(i);
    end
    if nDuplicateIafKeys > 0
        integrityIssues{end+1} = sprintf(['%d duplicate (Subject x BlockNumber x Condition x Unit) key(s) found in X14''s IAF output ' ...
            '(Task=%s, Level=%s) -- first occurrence kept for each, later duplicate(s) ignored.'], nDuplicateIafKeys, TASK_LABEL, IAF_LEVEL_KEY); %#ok<AGROW>
    end

    % Sanity check on IAF's own search range (informational; X14's own
    % code enforces this, this is just a confirming spot-check here).
    iafValsCheck = iafSubT.Value(~isnan(iafSubT.Value));
    if ~isempty(iafValsCheck) && (min(iafValsCheck) < IAF_SEARCH_RANGE_HZ(1) - 1e-9 || max(iafValsCheck) > IAF_SEARCH_RANGE_HZ(2) + 1e-9)
        integrityIssues{end+1} = sprintf(['Task=%s: at least one non-NaN IAF value from X14 falls outside the expected %g-%g Hz search range ' ...
            '(min=%.4f, max=%.4f). The "no clipping expected" assumption in this script''s header may not hold -- check the clipping count below.'], ...
            TASK_LABEL, IAF_SEARCH_RANGE_HZ(1), IAF_SEARCH_RANGE_HZ(2), min(iafValsCheck), max(iafValsCheck)); %#ok<AGROW>
    end

    %% Section 4a: filter psd_data to this Task, then to Condition in {E,N,R}
    taskMaskIdx = arrayfun(@(x) strcmp(x.TaskType, TASK_KEY), psd_data_ALL);
    thisTaskEntries = psd_data_ALL(taskMaskIdx);
    nBaselineExcluded = sum(arrayfun(@(x) strcmp(x.Condition, 'Baseline'), thisTaskEntries));
    condKeepMask = arrayfun(@(x) ismember(x.Condition, CONDITION_ORDER), thisTaskEntries);
    entries = thisTaskEntries(condKeepMask);
    fprintf('Task=%s: %d total recordings loaded, %d Baseline recordings excluded, %d E/N/R recordings retained.\n', ...
        TASK_LABEL, numel(thisTaskEntries), nBaselineExcluded, numel(entries));

    if isempty(entries)
        warning('Task=%s: no E/N/R recordings found after filtering -- skipping this Task entirely.', TASK_LABEL);
        continue;
    end

    %% Section 4b: data-driven design check (same check as X05/X14; runs on
    %  the full E/N/R recording set, before per-channel IAF-availability
    %  filtering)
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
        fprintf(2, '\n*** DESIGN VERIFICATION FLAGGED %d ISSUE(S) SO FAR (Task=%s) -- see report. ***\n\n', numel(integrityIssues), TASK_LABEL);
        for i = 1:numel(integrityIssues)
            fprintf(2, '  - %s\n', integrityIssues{i});
        end
    end
    if allSubjectsComplete3
        error('PhaseX20:DesignAssumptionViolated', ...
            'Task=%s: design verification found a complete 3-condition-per-subject structure. Stopping before fitting any model -- see integrityIssues.', TASK_LABEL);
    end

    %% Section 4c: per recording x Unit -- look up IAF, define the band,
    %  integrate MeanPSD, build AbsolutePower/LogAbsolutePower/RelativePower.
    fprintf('\nSection 4c: extracting individualized-alpha band power per recording x %s (Task=%s)...\n', ...
        lower(iif(strcmp(GRANULARITY,'PerChannel'),'channel','unit')), TASK_LABEL);

    nUnits = numel(UNITS);
    wideRows = cell(0, 11);   % Subject,BlockNumber,Condition,Unit,IafUsedHz,BandLowHz,BandHighHz,BandClipped,AbsolutePower,LogAbsolutePower,RelativePower (NaN + ExclusionReason handled via a parallel cell)
    exclusionReasonRows = {};   % parallel to wideRows, kept separate for the long/model table decision below
    longRows = cell(0, 6);      % Subject,BlockNumber,Condition,Unit,Metric,Value

    nAttempted = 0;
    nExcludedNotFound = 0;
    nExcludedNaN = 0;
    nClippedBand = 0;
    nNonPositiveAbsPowerChannels = 0;
    nRecordingsAllUnitsExcluded = 0;

    for e = 1:numel(entries)
        entry = entries(e);
        subjStr = entry.Subject; blockNum = entry.BlockNumber; condStr = entry.Condition;

        usableCountThisEntry = 0;

        for u = 1:nUnits
            unitName = UNITS{u};
            nAttempted = nAttempted + 1;

            key = local_iaf_key(subjStr, blockNum, condStr, unitName);
            if ~isKey(iafMap, key)
                nExcludedNotFound = nExcludedNotFound + 1;
                wideRows(end+1, :) = {subjStr, blockNum, condStr, unitName, NaN, NaN, NaN, false, NaN, NaN, NaN}; %#ok<AGROW>
                exclusionReasonRows{end+1, 1} = 'NoIAFRowInX14Output'; %#ok<AGROW>
                continue;
            end
            iafVal = iafMap(key);
            if isnan(iafVal)
                nExcludedNaN = nExcludedNaN + 1;
                wideRows(end+1, :) = {subjStr, blockNum, condStr, unitName, NaN, NaN, NaN, false, NaN, NaN, NaN}; %#ok<AGROW>
                exclusionReasonRows{end+1, 1} = 'IAFIsNaNInX14Output'; %#ok<AGROW>
                continue;
            end

            lowRaw = iafVal - IAF_HALF_WIDTH_HZ;
            highRaw = iafVal + IAF_HALF_WIDTH_HZ;
            lowClip = max(lowRaw, DENOM_RANGE_HZ(1));
            highClip = min(highRaw, DENOM_RANGE_HZ(2));
            wasClipped = (lowClip ~= lowRaw) || (highClip ~= highRaw);
            if wasClipped; nClippedBand = nClippedBand + 1; end

            [pxxRow, meanTotalPowerScalar] = local_get_unit_curve(entry, unitName, ChannelLabelsReference);
            absPow = local_integrate_band_power_row(FrequencyVectorHz, pxxRow, [lowClip, highClip]);

            logAbsPow = NaN;
            if ~isnan(absPow)
                if absPow <= 0
                    nNonPositiveAbsPowerChannels = nNonPositiveAbsPowerChannels + 1;
                else
                    logAbsPow = log(absPow);
                end
            end
            relPow = absPow ./ meanTotalPowerScalar;

            wideRows(end+1, :) = {subjStr, blockNum, condStr, unitName, iafVal, lowClip, highClip, wasClipped, absPow, logAbsPow, relPow}; %#ok<AGROW>
            exclusionReasonRows{end+1, 1} = ''; %#ok<AGROW>

            longRows(end+1, :) = {subjStr, blockNum, condStr, unitName, 'LogAbsolutePower_IndividualizedAlpha', logAbsPow}; %#ok<AGROW>
            longRows(end+1, :) = {subjStr, blockNum, condStr, unitName, 'RelativePower_IndividualizedAlpha', relPow}; %#ok<AGROW>

            usableCountThisEntry = usableCountThisEntry + 1;
        end

        if usableCountThisEntry == 0
            nRecordingsAllUnitsExcluded = nRecordingsAllUnitsExcluded + 1;
        end
    end

    if nNonPositiveAbsPowerChannels > 0
        integrityIssues{end+1} = sprintf(['%d recording x unit combination(s) had a non-positive AbsolutePower_IndividualizedAlpha -- ' ...
            'LogAbsolutePower_IndividualizedAlpha set to NaN there, not silently patched (Task=%s).'], nNonPositiveAbsPowerChannels, TASK_LABEL); %#ok<AGROW>
    end
    if nExcludedNotFound > 0 || nExcludedNaN > 0
        integrityIssues{end+1} = sprintf(['Task=%s: %d/%d recording x unit combination(s) excluded from the individualized-alpha family ' ...
            '(no usable IAF): %d had no matching IAF row in X14''s output, %d had a NaN IAF value there. ' ...
            '%d/%d combinations retained as usable.'], TASK_LABEL, nExcludedNotFound + nExcludedNaN, nAttempted, ...
            nExcludedNotFound, nExcludedNaN, nAttempted - nExcludedNotFound - nExcludedNaN, nAttempted); %#ok<AGROW>
    end
    if nRecordingsAllUnitsExcluded > 0
        integrityIssues{end+1} = sprintf(['Task=%s: %d recording(s) had EVERY unit excluded (0 usable channels anywhere in that recording).'], ...
            TASK_LABEL, nRecordingsAllUnitsExcluded); %#ok<AGROW>
    end
    if nClippedBand > 0
        integrityIssues{end+1} = sprintf(['Task=%s: %d recording x unit combination(s) had their individualized-alpha band clipped to the ' ...
            '[%g, %g] Hz boundary (NOT expected given X14''s 8-13 Hz IAF search range -- see the IAF-search-range sanity check above).'], ...
            TASK_LABEL, nClippedBand, DENOM_RANGE_HZ(1), DENOM_RANGE_HZ(2)); %#ok<AGROW>
    end

    fprintf('  Attempted %d recording x unit combination(s); %d usable (retained), %d excluded (no usable IAF); %d recording(s) fully excluded; %d band(s) clipped.\n\n', ...
        nAttempted, nAttempted - nExcludedNotFound - nExcludedNaN, nExcludedNotFound + nExcludedNaN, nRecordingsAllUnitsExcluded, nClippedBand);

    WideFeatureT = local_build_table(wideRows, [1 3 4], [2 5 6 7 8 9 10 11], ...
        {'Subject','BlockNumber','Condition','Unit','IAF_Used_Hz','BandLowHz','BandHighHz','BandClipped','AbsolutePower_IndividualizedAlpha','LogAbsolutePower_IndividualizedAlpha','RelativePower_IndividualizedAlpha'});
    WideFeatureT.BandClipped = logical(WideFeatureT.BandClipped);
    WideFeatureT.ExclusionReason = exclusionReasonRows;

    LongFeatureT = local_build_table(longRows, [1 3 4 5], [2 6], ...
        {'Subject','BlockNumber','Condition','Unit','Metric','Value'});

    % ---- NEW: apply the shared noisy-subject exclusion list (see
    % EXCLUSION_CONFIG_intervention.m), only when RUN_MODE=sensitivity_excl_noisy.
    % Never applied to the full_sample run. Applied to LongFeatureT --
    % the table the LMMs below are actually fit on -- AFTER the IAF-based
    % individualized-alpha band power is built above and BEFORE any model
    % fitting, identically regardless of which IAF granularity branch
    % (PerChannel, confirmed; or the defensive WholeScalp fallback) is
    % active, so both branches remain correct if the fallback ever
    % triggers. WideFeatureT (the full-transparency audit table, including
    % excluded-for-no-usable-IAF rows) is intentionally left untouched. ----
    nExcludedRowsSensitivity = 0;
    if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
        subjNumForExcl = local_subject_numeric(LongFeatureT.Subject);
        excludeMaskSensitivity = ismember(subjNumForExcl, EXCLUDE_NOISY_SUBJECTS);
        nExcludedRowsSensitivity = sum(excludeMaskSensitivity);
        LongFeatureT = LongFeatureT(~excludeMaskSensitivity, :);
        fprintf('Task=%s: RUN_MODE=%s -- excluded %d row(s) belonging to subject(s) {%s}; %d rows remain for modeling (%s granularity).\n', ...
            TASK_LABEL, RUN_MODE_LABEL, nExcludedRowsSensitivity, ...
            strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','), height(LongFeatureT), GRANULARITY);
    end

    %% Section 4d: fit LMMs per Unit x Metric
    fprintf('--------------------------------------------------------\n');
    fprintf('SECTION 4d: MODEL FITTING (Task=%s, %s level, Unit x Metric)\n', TASK_LABEL, GRANULARITY);
    fprintf('--------------------------------------------------------\n');

    omnibusRows = cell(0, 6);   % Unit,Metric,FStat,DF1,DF2,pValue
    lmeMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    nSkippedTooFewObs = 0;

    for m = 1:numel(METRICS)
        metricName = METRICS{m};
        for u = 1:nUnits
            unitName = UNITS{u};
            subT = LongFeatureT(strcmp(LongFeatureT.Metric, metricName) & strcmp(LongFeatureT.Unit, unitName), :);
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
                    convergenceWarnings{end+1} = sprintf('%s / %s: %s', unitName, metricName, warnMsg); %#ok<AGROW>
                end
            catch ME
                integrityIssues{end+1} = sprintf('fitlme FAILED for %s / %s (Task=%s): %s', unitName, metricName, TASK_LABEL, ME.message); %#ok<AGROW>
                continue;
            end
            aov = anova(lme, 'DFMethod', 'satterthwaite');
            rowIdx = find(strcmp(aov.Term, 'Condition'), 1);
            if isempty(rowIdx)
                integrityIssues{end+1} = sprintf('anova(lme) has no Condition row for %s / %s (Task=%s).', unitName, metricName, TASK_LABEL); %#ok<AGROW>
                continue;
            end
            key = sprintf('%s|%s', unitName, metricName);
            lmeMap(key) = lme;
            omnibusRows(end+1, :) = {unitName, metricName, aov.FStat(rowIdx), aov.DF1(rowIdx), aov.DF2(rowIdx), aov.pValue(rowIdx)}; %#ok<AGROW>
        end
    end
    if nSkippedTooFewObs > 0
        integrityIssues{end+1} = sprintf('%d (Unit x Metric) combination(s) had fewer than 4 non-missing observations and were skipped (Task=%s).', nSkippedTooFewObs, TASK_LABEL); %#ok<AGROW>
    end

    OmnibusT = local_build_table(omnibusRows, [1 2], [3 4 5 6], {'Unit','Metric','FStat','DF1','DF2','pValue'});
    OmnibusT.Level = repmat({GRANULARITY}, height(OmnibusT), 1);
    fprintf('Model fitting complete (Task=%s): %d omnibus tests fit successfully.\n', TASK_LABEL, height(OmnibusT));

    %% Section 4e: BH-FDR within each Metric family (family = across Units);
    OmnibusT.pFDR = nan(height(OmnibusT), 1);
    OmnibusT.SignificantFDR = false(height(OmnibusT), 1);
    OmnibusT.FamilySize = nan(height(OmnibusT), 1);
    OmnibusT.SingleTestNotFDRCorrected = false(height(OmnibusT), 1);

    for m = 1:numel(METRICS)
        metricName = METRICS{m};
        idx = strcmp(OmnibusT.Metric, metricName);
        nInFamily = sum(idx);
        if nInFamily == 0; continue; end
        OmnibusT.FamilySize(idx) = nInFamily;
        pvals = OmnibusT.pValue(idx);
        if nInFamily == 1
            OmnibusT.pFDR(idx) = pvals;
            OmnibusT.SignificantFDR(idx) = pvals <= FDR_ALPHA;
            OmnibusT.SingleTestNotFDRCorrected(idx) = true;
        else
            [adjP, reject] = local_fdr_bh(pvals, FDR_ALPHA);
            OmnibusT.pFDR(idx) = adjP;
            OmnibusT.SignificantFDR(idx) = reject;
        end
    end
    fprintf('\n--- Individualized-alpha omnibus results (Task=%s, %s) ---\n', TASK_LABEL, GRANULARITY);
    disp(OmnibusT);

    %% Section 4f: pairwise contrasts (E-N, E-R, N-R), Holm-corrected.
    %  Gating: FDR survival when family size>1; nominal p<=0.05 when family
    %  size==1 (there is no FDR family to survive)
    pairwiseRows = cell(0, 11);   % Unit,Metric,Contrast,Estimate,SE,DF,tStat,pValue,pHolm,CI_Low,CI_High
    if any(OmnibusT.SingleTestNotFDRCorrected)
        gateMask = OmnibusT.pValue <= FDR_ALPHA;
    else
        gateMask = OmnibusT.SignificantFDR;
    end
    for row = find(gateMask)'
        unitName = OmnibusT.Unit{row};
        metricName = OmnibusT.Metric{row};
        key = sprintf('%s|%s', unitName, metricName);
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
            pairwiseRows(end+1, :) = [{unitName, metricName}, contrastData(c, 1:6), {adjP(c)}, contrastData(c, 7:8)]; %#ok<AGROW>
        end
    end
    PairwiseT = local_build_table(pairwiseRows, [1 2 3], [4 5 6 7 8 9 10 11], ...
        {'Unit','Metric','Contrast','Estimate','SE','DF','tStat','pValue','pHolm','CI_Low','CI_High'});
    fprintf('\n--- Pairwise contrasts (Task=%s; Unit x Metric gated as described above) ---\n', TASK_LABEL);
    disp(PairwiseT);

    %% Section 4g: residual diagnostics -- bounded to gated (significant)
    %  Unit x Metric combinations only
    fprintf('\n--------------------------------------------------------\n');
    fprintf('SECTION 4g: MODEL DIAGNOSTICS (Task=%s, bounded to significant results)\n', TASK_LABEL);
    fprintf('--------------------------------------------------------\n');
    nDiagSaved = 0;
    for row = find(gateMask)'
        unitName = OmnibusT.Unit{row};
        metricName = OmnibusT.Metric{row};
        key = sprintf('%s|%s', unitName, metricName);
        if ~isKey(lmeMap, key); continue; end
        lme = lmeMap(key);
        try
            res = residuals(lme, 'ResidualType', 'Pearson');
            fig = figure('Visible', 'off', 'Position', [100 100 900 400]);
            subplot(1,2,1);
            histogram(res, 15);
            title(sprintf('Residuals: %s %s, Task=%s', unitName, metricName, TASK_LABEL), 'Interpreter', 'none');
            xlabel('Pearson residual'); ylabel('Count');
            subplot(1,2,2);
            qqplot(res);
            title('Normal Q-Q');
            outPng = fullfile(STATS_DIR, sprintf('diagnostic_residuals_%s_%s.png', unitName, metricName));
            saveas(fig, outPng);
            close(fig);
            nDiagSaved = nDiagSaved + 1;
        catch ME
            integrityIssues{end+1} = sprintf('Residual diagnostic failed for %s / %s (Task=%s): %s', unitName, metricName, TASK_LABEL, ME.message); %#ok<AGROW>
        end
    end
    fprintf('Saved %d residual diagnostic plot(s) (Task=%s, bounded to gated/significant Unit x Metric combinations).\n', nDiagSaved, TASK_LABEL);

    %% Section 4h: PARAMS20 and outputs
    PARAMS20 = struct();
    PARAMS20.ScriptName = 'X20_individualized_alpha_intervention.m';
    PARAMS20.Task = TASK_LABEL;
    PARAMS20.GeneratedOn = datestr(now); %#ok<TNOW1,DATST>
    PARAMS20.AnalysisLevel = sprintf('Individualized-alpha (IAF +/- %g Hz) band power from X03 PSD curves, Task=%s only, Condition in {E,N,R} only (Baseline excluded)', IAF_HALF_WIDTH_HZ, TASK_LABEL);
    PARAMS20.Model = 'LinearMixedEffects_RandomIntercept';
    PARAMS20.FixedEffectsPrimary = 'Condition (3-level: E, N, R; reference = E)';
    PARAMS20.RandomEffects = '(1|Subject)';
    PARAMS20.FitMethod = 'REML';
    PARAMS20.DFMethod = 'Satterthwaite';
    PARAMS20.Metrics = METRICS;
    PARAMS20.IAFHalfWidthHz = IAF_HALF_WIDTH_HZ;
    PARAMS20.IAFHalfWidthConfirmation = 'Default half-width of the individualized alpha band.';
    PARAMS20.IAFSourceFile = iafFile;
    PARAMS20.IAFVariantUsed = IAF_VARIANT_USED;
    PARAMS20.IAFVariantUsedRationale = 'Paired with MeanPSD (the SAME curve type used here for band-power integration); MedianPSDBased IAF (if present in X14''s output) is deliberately not used, to avoid mixing two different underlying curves'' worth of information into one metric.';
    PARAMS20.GranularityDetected = GRANULARITY;
    PARAMS20.GranularityDetectionMethod = 'Data-driven, per Task: inspected which Level values are present in X14''s FeatureValueT for Metric=IAF, Variant=MeanPSDBased. Electrode present -> PerChannel; else Global present -> WholeScalp (defensive fallback, not expected).';
    PARAMS20.ClippingRule = sprintf('IndividualAlphaBand = [IAF-%g, IAF+%g] Hz, each edge independently clipped to the [%g, %g] Hz RelativePowerDenominatorRangeHz boundary if it would otherwise fall outside it.', IAF_HALF_WIDTH_HZ, IAF_HALF_WIDTH_HZ, DENOM_RANGE_HZ(1), DENOM_RANGE_HZ(2));
    PARAMS20.RelativePowerDenominatorRangeHz = DENOM_RANGE_HZ;
    PARAMS20.BandEdgeConvention = 'both edges inclusive, same trapz-mask convention as X03''s integrate_band_power / PARAMS.BandEdgeConvention';
    PARAMS20.ExclusionPolicy = 'Recording x Unit combinations with no matching IAF row in X14''s output, or a NaN IAF value there, are excluded entirely from this family (never imputed). Counts always reported per Task.';
    PARAMS20.ExclusionCounts = struct('nAttempted', nAttempted, 'nExcludedNotFoundInX14', nExcludedNotFound, ...
        'nExcludedNaNIAF', nExcludedNaN, 'nUsable', nAttempted - nExcludedNotFound - nExcludedNaN, ...
        'nRecordingsAllUnitsExcluded', nRecordingsAllUnitsExcluded, 'nBandsClipped', nClippedBand);
    PARAMS20.FDRAlpha = FDR_ALPHA;
    PARAMS20.FDRFamilyScheme = 'Separate BH-FDR family per Metric, across Units (channels, or the single Global unit in the WholeScalp fallback). Family size 1 is reported as-is, no correction applied or needed (same convention as X05''s Primary tier / X14''s Global level). Own family, never pooled with any other script''s metrics.';
    PARAMS20.PairwiseContrastMethod = 'Holm-Bonferroni across the 3 pairwise contrasts (E-N, E-R, N-R -- standardized on the X05/X06/X11-X15 convention). Gated on BH-FDR survival when family size>1, or on nominal p<=0.05 when family size==1 (no FDR family exists to survive) -- see header for full disclosure of this distinction.';
    PARAMS20.ScopeDecisionNoRoiOrGlobal = 'This script produces ONLY the Electrode (per-channel) level (or the WholeScalp fallback, if ever detected). ROI- and Global-level individualized-alpha power are explicitly out of scope -- see script header for why X14''s frequency-value ROI/Global-averaging argument does not transfer to a power value integrated over a per-channel-specific band.';
    PARAMS20.BaselineExclusionNote = sprintf('%d Baseline recording(s) excluded before any modeling; Baseline is processed by X03 but never enters this or any other inferential model in this pipeline.', nBaselineExcluded);
    PARAMS20.SourcePsdFile = psdFile;
    PARAMS20.DesignVerification = struct('nRecordings', nRecordings, 'nSubjects', nSubjects, ...
        'PairTypeCounts', {pairTypeCount.keys}, 'IntegrityIssues', {integrityIssues});
    PARAMS20.RunDateTime = datestr(now); %#ok<TNOW1,DATST>
    PARAMS20.RunMode = RUN_MODE;
    if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
        PARAMS20.ExcludedSubjects = EXCLUDE_NOISY_SUBJECTS;
        PARAMS20.ExcludedSubjectsRationale = EXCLUSION_RATIONALE_intervention;
        PARAMS20.RowsExcludedForSensitivity = nExcludedRowsSensitivity;
    else
        PARAMS20.ExcludedSubjects = [];
    end

    save(fullfile(STATS_DIR, 'phase20_parameters.mat'), 'PARAMS20');

    save(fullfile(STATS_DIR, 'phase20_feature_values.mat'), 'WideFeatureT');
    writetable(WideFeatureT, fullfile(STATS_DIR, 'phase20_feature_values.xlsx'));
    writetable(WideFeatureT, fullfile(STATS_DIR, 'phase20_feature_values.csv'));

    save(fullfile(STATS_DIR, 'phase20_feature_values_long.mat'), 'LongFeatureT');
    writetable(LongFeatureT, fullfile(STATS_DIR, 'phase20_feature_values_long.xlsx'));
    writetable(LongFeatureT, fullfile(STATS_DIR, 'phase20_feature_values_long.csv'));

    save(fullfile(STATS_DIR, 'phase20_omnibus_results.mat'), 'OmnibusT');
    writetable(OmnibusT, fullfile(STATS_DIR, 'phase20_omnibus_results.xlsx'));
    writetable(OmnibusT, fullfile(STATS_DIR, 'phase20_omnibus_results.csv'));

    save(fullfile(STATS_DIR, 'phase20_pairwise_contrasts.mat'), 'PairwiseT');
    if ~isempty(PairwiseT)
        writetable(PairwiseT, fullfile(STATS_DIR, 'phase20_pairwise_contrasts.xlsx'));
        writetable(PairwiseT, fullfile(STATS_DIR, 'phase20_pairwise_contrasts.csv'));
    end

    fprintf('\nSaved output to: %s\n', STATS_DIR);

    %% Section 4i: Report
    reportLines = {};
    reportLines{end+1} = repmat('=', 1, 60);
    reportLines{end+1} = sprintf('X20 -- INDIVIDUALIZED ALPHA (IAF +/- %g Hz) BAND-POWER FAMILY REPORT (Task=%s, RUN_MODE=%s)', IAF_HALF_WIDTH_HZ, TASK_LABEL, RUN_MODE_LABEL);
    reportLines{end+1} = repmat('=', 1, 60);
    reportLines{end+1} = '';
    reportLines{end+1} = 'Exploratory family (IAF +/- 2 Hz width).';
    reportLines{end+1} = 'Never reads or modifies any other script''s outputs. Reuses X14''s already-computed IAF read-only.';
    reportLines{end+1} = '';
    reportLines{end+1} = sprintf('IAF granularity detected (data-driven, from X14''s own output): %s%s', GRANULARITY, iif(strcmp(GRANULARITY,'WholeScalp'), ' *** DEFENSIVE FALLBACK, NOT the expected case -- see header ***', ''));
    reportLines{end+1} = sprintf('IAF variant used: %s (paired deliberately with MeanPSD-based power integration -- see header).', IAF_VARIANT_USED);
    reportLines{end+1} = '';
    reportLines{end+1} = sprintf('Baseline recordings excluded before modeling: %d', nBaselineExcluded);
    if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
        reportLines{end+1} = sprintf('RUN_MODE=%s: %d additional row(s) excluded from LongFeatureT (modeling table) belonging to subject(s) {%s} (see EXCLUSION_CONFIG_intervention.m).', ...
            RUN_MODE_LABEL, nExcludedRowsSensitivity, strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','));
    else
        reportLines{end+1} = sprintf('RUN_MODE=%s: no subjects excluded.', RUN_MODE_LABEL);
    end
    reportLines{end+1} = sprintf('Recordings verified (design check): %d | Subjects: %d', nRecordings, nSubjects);
    reportLines{end+1} = sprintf('Condition-pair coverage: %s', strjoin(cellfun(@(k) sprintf('%s=%d', strrep(k,'_',' & '), pairTypeCount(k)), pairKeys, 'UniformOutput', false), ', '));
    reportLines{end+1} = '';
    reportLines{end+1} = '--- No-usable-IAF exclusion policy and counts ---';
    reportLines{end+1} = sprintf('  Recording x unit combinations attempted     : %d', nAttempted);
    reportLines{end+1} = sprintf('  Excluded -- no IAF row in X14 output        : %d', nExcludedNotFound);
    reportLines{end+1} = sprintf('  Excluded -- IAF NaN in X14 output           : %d', nExcludedNaN);
    reportLines{end+1} = sprintf('  Usable / retained                           : %d', nAttempted - nExcludedNotFound - nExcludedNaN);
    reportLines{end+1} = sprintf('  Whole recordings with ALL units excluded    : %d', nRecordingsAllUnitsExcluded);
    reportLines{end+1} = sprintf('  Bands clipped to [%g, %g] Hz boundary        : %d (expected 0, given IAF''s 8-13 Hz search range)', DENOM_RANGE_HZ(1), DENOM_RANGE_HZ(2), nClippedBand);
    reportLines{end+1} = '';
    reportLines{end+1} = 'Model: Linear mixed-effects, Condition (fixed, 3-level) + (1|Subject), REML, Satterthwaite DF.';
    reportLines{end+1} = 'Metrics: LogAbsolutePower_IndividualizedAlpha and RelativePower_IndividualizedAlpha, co-exploratory, never pooled for FDR.';
    familySizeForReport = NaN;
    if ~isempty(OmnibusT); familySizeForReport = max(OmnibusT.FamilySize); end
    reportLines{end+1} = sprintf('FDR family scheme: separate BH-FDR per Metric, across Units; family size 1 reported as-is, no correction (see header). This Task''s family size = %d.', familySizeForReport);
    reportLines{end+1} = 'Pairwise contrast convention: E-N, E-R, N-R (standardized on X05/X06/X11-X15).';
    if any(OmnibusT.SingleTestNotFDRCorrected)
        reportLines{end+1} = 'Pairwise gating this Task: NOMINAL significance (p<=0.05) -- family size is 1, so there is no FDR family to survive.';
    else
        reportLines{end+1} = 'Pairwise gating this Task: BH-FDR survival (SignificantFDR==true).';
    end
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Omnibus results (all Unit x Metric tests) ---';
    for r = 1:height(OmnibusT)
        sigTag = local_yn(OmnibusT.SignificantFDR(r));
        if OmnibusT.SingleTestNotFDRCorrected(r)
            sigTag = [sigTag ' (single test, NOT FDR-corrected)'];
        end
        reportLines{end+1} = sprintf('  [%-8s] %-40s F(%g,%.1f)=%.3f p=%.4f pFDR=%.4f family_n=%d %s', ...
            OmnibusT.Unit{r}, OmnibusT.Metric{r}, OmnibusT.DF1(r), OmnibusT.DF2(r), ...
            OmnibusT.FStat(r), OmnibusT.pValue(r), OmnibusT.pFDR(r), OmnibusT.FamilySize(r), sigTag); %#ok<AGROW>
    end
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Pairwise contrasts (Unit x Metric gated as described above; Holm-corrected within unit) ---';
    if isempty(PairwiseT)
        reportLines{end+1} = '  (none -- no Unit x Metric combination met the gating criterion above)';
    else
        for r = 1:height(PairwiseT)
            reportLines{end+1} = sprintf('  [%-8s] %-40s %-8s est=%.4f (95%% CI %.4f to %.4f) t=%.3f df=%.1f p=%.4f pHolm=%.4f', ...
                PairwiseT.Unit{r}, PairwiseT.Metric{r}, PairwiseT.Contrast{r}, PairwiseT.Estimate(r), ...
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
    fid = fopen(fullfile(REPORT_DIR, sprintf('X20_individualized_alpha_report_%s_%s.txt', TASK_LABEL, RUN_MODE_LABEL)), 'w');
    fprintf(fid, '%s\n', reportText);
    fclose(fid);

end % TASK_TYPES loop

fprintf('\n=== X20 complete for all three Tasks. Purely additive -- no other script''s outputs touched. ===\n');

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
            error('PhaseX20:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name1);
        end
        cRow(idx1) = cRow(idx1) - 1;
    end
    if ~strcmp(level2, refLevel)
        name2 = sprintf('%s_%s', varName, level2);
        idx2 = find(strcmp(coeffNames, name2));
        if isempty(idx2)
            error('PhaseX20:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name2);
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

function T = local_build_table(rows, textColIdx, numColIdx, varNames)
    nCols = numel(varNames);
    if size(rows, 2) ~= nCols
        error('PhaseX20:BuildTableColumnMismatch', ...
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
            error('PhaseX20:BuildTableColumnMismatch', 'Column %d not classified as text or numeric.', j);
        end
    end
    T = table(cols{:}, 'VariableNames', varNames);
end

function k = local_iaf_key(subjStr, blockNum, condStr, unitStr)
    % Builds a plain-string join key for the (Subject x BlockNumber x
    % Condition x Unit) lookup into X14's IAF table. Accepts either raw
    % char/numeric values (as stored on a psd_data entry) or 1x1 cell
    % contents (as returned by indexing a table's cell-backed column).
    if iscell(subjStr); subjStr = subjStr{1}; end
    if iscell(condStr); condStr = condStr{1}; end
    if iscell(unitStr); unitStr = unitStr{1}; end
    k = sprintf('%s|%g|%s|%s', char(subjStr), double(blockNum), char(condStr), upper(char(unitStr)));
end

function bp = local_integrate_band_power_row(freqVec, pxxRow, bandRange)
    % Integrates a single (1 x nFreq) power row over bandRange using the
    % IDENTICAL "both edges inclusive" mask-then-trapz logic as X03's own
    % integrate_band_power helper and PARAMS.BandEdgeConvention. X03's
    % helper operates on Pxx oriented nFreq x nChan (trapz along dim 1);
    % this helper is adapted only for the nChan x nFreq orientation
    % entry.MeanPSD is actually stored in (a single already-selected row,
    % i.e. one channel or one whole-scalp average) -- the mask and the
    % trapz arithmetic themselves are unchanged.
    freqVec = freqVec(:)';
    pxxRow = pxxRow(:)';
    mask = freqVec >= bandRange(1) & freqVec <= bandRange(2);
    if ~any(mask) || any(isnan(pxxRow(mask))) || any(~isfinite(pxxRow(mask)))
        bp = NaN;
        return;
    end
    bp = trapz(freqVec(mask), pxxRow(mask));
end

function [pxxRow, meanTotalPowerScalar] = local_get_unit_curve(entry, unitName, ChannelLabelsReference)
    validWin = entry.WindowValid;
    if strcmpi(unitName, 'Global')
        pxxRow = mean(entry.MeanPSD, 1, 'omitnan');
        if any(validWin)
            meanTotalPowerScalar = mean(mean(entry.TotalPower_1to45Hz(validWin, :), 1, 'omitnan'), 2, 'omitnan');
        else
            meanTotalPowerScalar = NaN;
        end
    else
        chIdx = find(strcmpi(ChannelLabelsReference, unitName), 1);
        if isempty(chIdx)
            pxxRow = nan(1, size(entry.MeanPSD, 2));
            meanTotalPowerScalar = NaN;
            return;
        end
        pxxRow = entry.MeanPSD(chIdx, :);
        if any(validWin)
            meanTotalPowerScalar = mean(entry.TotalPower_1to45Hz(validWin, chIdx), 1);
        else
            meanTotalPowerScalar = NaN;
        end
    end
end

function out = iif(cond, valTrue, valFalse)
    % Tiny inline-conditional helper used only for report-string wording
    % (never for any numeric/statistical decision in this script).
    if cond
        out = valTrue;
    else
        out = valFalse;
    end
end
