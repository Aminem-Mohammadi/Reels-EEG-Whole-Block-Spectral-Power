%% ========================================================================
% X17_half_split_condition_time_intervention.m
%
% Condition x recording-half model: band power recomputed for the first and second half
% of each recording and analysed with Metric ~ Condition*TimeHalf + (1|Subject).
% Two observations per recording enter the model without modelling their dependence, so
% results are supportive only.
% Inputs : 03_psd.
% Outputs: 17_half_split_condition_time/<RUN_MODE>/, 09_reports/.
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

psdMatFile = fullfile(PSD_DIR, 'intervention_psd.mat');

if ~exist(psdMatFile, 'file')
    error('Phase17:MissingInput', ['X03 output not found: %s\n' ...
        'Run X03_compute_intervention_psd.m first.'], psdMatFile);
end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

addpath(fileparts(mfilename('fullpath')));
run(fullfile(fileparts(mfilename('fullpath')), 'EXCLUSION_CONFIG_intervention.m'));

RUN_MODE = 'full_sample';   % 'full_sample' (N = 30) or 'sensitivity_excl_noisy' (N = 28)
if ~ismember(RUN_MODE, RUN_MODES_intervention)
    error('Phase17:BadRunMode', 'RUN_MODE must be one of: %s', strjoin(RUN_MODES_intervention, ', '));
end
RUN_MODE_LABEL = RUN_MODES_DISPLAY_intervention.(RUN_MODE);
fprintf('RUN_MODE = %s\n', RUN_MODE_LABEL);
if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
    fprintf('Excluding subject(s) {%s} from every model in this run -- see EXCLUSION_CONFIG_intervention.m for rationale.\n', ...
        strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','));
end
% ---- END NEW ----

if isempty(ver('stats'))
    error('Phase17:MissingToolbox', ...
        'Statistics and Machine Learning Toolbox not detected. This script requires fitlme, anova(LinearMixedModel), and coefTest.');
end

%% Section 1: Load X03 output (read-only; the ONLY input this script has)
S3 = intervention_load(psdMatFile, 'psd_data', 'ChannelLabelsReference', 'PARAMS');
psd_data = S3.psd_data;
ChannelLabelsReference = S3.ChannelLabelsReference(:);
PARAMS2 = S3.PARAMS;
clear S3;

if ~isfield(psd_data, 'WindowStartSec')
    error('Phase17:MissingField', 'psd_data is missing .WindowStartSec -- required to order windows chronologically before splitting.');
end

ALL_BANDS = fieldnames(PARAMS2.BandDefinitionsHz);
nBands = numel(ALL_BANDS);
nChanRef = numel(ChannelLabelsReference);

% --- Band tiers, IDENTICAL to X05/X09/X11/X18's own locked tier lists,
%     used for the GLOBAL-level FDR family boundary (Tier x Metric,
%     pooling across bands within a tier) and, at ELECTRODE level, as a
%     LABEL ONLY (never a family boundary there). ---
PRIMARY_TIER_BANDS     = {'Alpha'};
SECONDARY_TIER_BANDS   = {'Delta', 'Theta', 'Beta'};
EXPLORATORY_TIER_BANDS = {'Gamma_low_exploratory', 'Gamma_high1_exploratory', 'Gamma_high2_exploratory'};
EXPLORATORY_BETA_TIER_BANDS = {'Beta_low_exploratory', 'Beta_high_exploratory'};
expectedTierBands = [PRIMARY_TIER_BANDS, SECONDARY_TIER_BANDS, EXPLORATORY_TIER_BANDS, EXPLORATORY_BETA_TIER_BANDS];
if ~isempty(setxor(sort(expectedTierBands), sort(ALL_BANDS')))
    error('Phase17:BandMismatch', ...
        ['The hard-coded band tiers do not exactly match PARAMS2.BandDefinitionsHz. ' ...
         'Tiers: %s | Phase3 bands: %s'], strjoin(expectedTierBands, ','), strjoin(ALL_BANDS', ','));
end
tierDefs = {'Primary', PRIMARY_TIER_BANDS; 'Secondary', SECONDARY_TIER_BANDS; 'Exploratory', EXPLORATORY_TIER_BANDS; 'ExploratoryBeta', EXPLORATORY_BETA_TIER_BANDS};
tierOf = containers.Map();
for i = 1:numel(PRIMARY_TIER_BANDS);     tierOf(PRIMARY_TIER_BANDS{i}) = 'Primary'; end     %#ok<NOCOL>
for i = 1:numel(SECONDARY_TIER_BANDS);   tierOf(SECONDARY_TIER_BANDS{i}) = 'Secondary'; end %#ok<NOCOL>
for i = 1:numel(EXPLORATORY_TIER_BANDS); tierOf(EXPLORATORY_TIER_BANDS{i}) = 'Exploratory'; end %#ok<NOCOL>
for i = 1:numel(EXPLORATORY_BETA_TIER_BANDS); tierOf(EXPLORATORY_BETA_TIER_BANDS{i}) = 'ExploratoryBeta'; end

METRICS = {'LogAbsolutePower', 'RelativePower'};   % co-primary; NEVER pooled for FDR (see header)
CONDITION_ORDER = {'E', 'N', 'R'};
REF_COND = CONDITION_ORDER{1};
CONTRASTS = {'E','N','E - N'; 'E','R','E - R'; 'N','R','N - R'};   % {level2, level1, label}; estimate = level2 - level1

TIMEHALF_ORDER = {'1', '2'};
REF_TIMEHALF = TIMEHALF_ORDER{1};
NONREF_TIMEHALF = TIMEHALF_ORDER{2};

MIN_VALID_WINDOWS_PER_RECORDING = 4;

FDR_ALPHA = 0.05;

TASK_TYPES = {'intervention'};
TASK_DISPLAY = struct('intervention', 'Intervention');

fprintf('========================================================\n');
fprintf('PHASE 17 (X17) -- CONDITION x TIME-HALF EXPLORATORY ANALYSIS, INTERVENTION DATA\n');
fprintf('Single intervention analysis (conditions E / N / R).\n');
fprintf('========================================================\n\n');

%% ========================================================================
for taskI = 1:numel(TASK_TYPES)
    TASK_KEY = TASK_TYPES{taskI};
    TASK_LABEL = TASK_DISPLAY.(TASK_KEY);

    STATS_DIR = fullfile(ANALYSIS_ROOT, '17_half_split_condition_time', RUN_MODE);
    if ~exist(STATS_DIR, 'dir'); mkdir(STATS_DIR); end

    fprintf('\n########################################################\n');
    fprintf('### TASK = %s | RUN_MODE = %s\n', TASK_LABEL, RUN_MODE_LABEL);
    fprintf('########################################################\n\n');

    integrityIssues = {};

    %% Section 2: filter psd_data to this Task, then to Condition in {E,N,R}
    taskMaskPSD = strcmp({psd_data.TaskType}, TASK_KEY);
    taskEntries = psd_data(taskMaskPSD);
    nBaselineExcluded = sum(strcmp({taskEntries.Condition}, 'Baseline'));
    condMaskPSD = ismember({taskEntries.Condition}, CONDITION_ORDER);
    taskEntriesENR = taskEntries(condMaskPSD);
    fprintf('Task=%s: %d recordings loaded from X03, %d Baseline recordings excluded, %d E/N/R recordings retained.\n', ...
        TASK_LABEL, numel(taskEntries), nBaselineExcluded, numel(taskEntriesENR));

    if isempty(taskEntriesENR)
        warning('Task=%s: no E/N/R recordings found -- skipping this Task entirely.', TASK_LABEL);
        continue;
    end

    %% ====================================================================
    %  SECTION 3: HALF-SPLIT + AGGREGATE EACH RECORDING (X04's exact
    %  aggregation formula, applied once per half).
    % ====================================================================
    elecRows = cell(0, 10);
    elecCols = {'Subject','TaskType','BlockNumber','Condition','TimeHalf','Channel','Band', ...
        'AbsolutePower','LogAbsolutePower','RelativePower'};
    globRows = cell(0, 9);
    globCols = {'Subject','TaskType','BlockNumber','Condition','TimeHalf','Band', ...
        'AbsolutePower','LogAbsolutePower','RelativePower'};
    countRows = cell(0, 8);
    countCols = {'Subject','TaskType','BlockNumber','Condition','NumValidWindowsTotal', ...
        'Half1Count','Half2Count','Included'};

    nExcludedTooFewWindows = 0;
    nReordered = 0;
    nNonPositiveLogFlagged = 0;

    for i = 1:numel(taskEntriesENR)
        entry = taskEntriesENR(i);
        subjStr = entry.Subject; blockNum = entry.BlockNumber; condStr = entry.Condition;

        validIdx = find(entry.WindowValid);
        n = numel(validIdx);

        if n < MIN_VALID_WINDOWS_PER_RECORDING
            countRows(end+1, :) = {subjStr, TASK_KEY, blockNum, condStr, n, NaN, NaN, false}; %#ok<AGROW>
            integrityIssues{end+1} = sprintf(['Sub%s Block%g Cond%s (Task=%s): only %d valid window(s) ' ...
                '(< minimum of %d, see Judgment Call 1) -- EXCLUDED from this analysis entirely.'], ...
                subjStr, blockNum, condStr, TASK_LABEL, n, MIN_VALID_WINDOWS_PER_RECORDING); %#ok<AGROW>
            nExcludedTooFewWindows = nExcludedTooFewWindows + 1;
            continue;
        end

        startsHere = entry.WindowStartSec(validIdx);
        [sortedStarts, sortOrder] = sort(startsHere, 'ascend');
        if ~isequal(sortOrder(:), (1:n)')
            nReordered = nReordered + 1;
            integrityIssues{end+1} = sprintf(['Sub%s Block%g Cond%s (Task=%s): valid windows were NOT already in ' ...
                'chronological order by WindowStartSec -- reordered defensively before splitting.'], ...
                subjStr, blockNum, condStr, TASK_LABEL); %#ok<AGROW>
        end
        orderedValidIdx = validIdx(sortOrder);

        half1Count = ceil(n / 2);
        half2Count = n - half1Count;
        half1Idx = orderedValidIdx(1:half1Count);
        half2Idx = orderedValidIdx(half1Count+1:end);
        countRows(end+1, :) = {subjStr, TASK_KEY, blockNum, condStr, n, half1Count, half2Count, true}; %#ok<AGROW>

        halves = {half1Idx, half2Idx};
        for h = 1:2
            hIdx = halves{h};
            timeHalfLabel = TIMEHALF_ORDER{h};

            absMeanAllBands = nan(1, nChanRef);
            relMeanAllBands = nan(1, nChanRef);
            for b = 1:nBands
                bandName = ALL_BANDS{b};
                absVals = entry.AbsoluteBandPower.(bandName)(hIdx, :);

                % --- X04's EXACT aggregation formula, unchanged, applied
                %     to this half's windows only ---
                absMean = mean(absVals, 1);
                meanTotalHalf = mean(entry.TotalPower_1to45Hz(hIdx, :), 1);
                relMean = absMean ./ meanTotalHalf;
                absLog = log(absMean);
                if any(absMean <= 0)
                    nBad = sum(absMean <= 0);
                    nNonPositiveLogFlagged = nNonPositiveLogFlagged + nBad;
                    integrityIssues{end+1} = sprintf(['Sub%s Block%g Cond%s Band%s TimeHalf%s (Task=%s): non-positive ' ...
                        'AbsolutePower found in %d channel(s) -- LogAbsolutePower set to NaN there, not silently patched.'], ...
                        subjStr, blockNum, condStr, bandName, timeHalfLabel, TASK_LABEL, nBad); %#ok<AGROW>
                    absLog(absMean <= 0) = NaN;
                end

                for ch = 1:nChanRef
                    elecRows(end+1, :) = {subjStr, TASK_KEY, blockNum, condStr, timeHalfLabel, ...
                        ChannelLabelsReference{ch}, bandName, absMean(ch), absLog(ch), relMean(ch)}; %#ok<AGROW>
                end

                % --- Global = arithmetic mean across channels of the
                %     per-channel means (X04's exact Global formula) ---
                globalAbs = mean(absMean);
                globalRel = mean(relMean);
                globalLog = log(globalAbs);
                if globalAbs <= 0
                    integrityIssues{end+1} = sprintf(['Sub%s Block%g Cond%s Band%s TimeHalf%s (Task=%s): non-positive ' ...
                        'whole-scalp AbsolutePower -- LogAbsolutePower set to NaN.'], ...
                        subjStr, blockNum, condStr, bandName, timeHalfLabel, TASK_LABEL); %#ok<AGROW>
                    globalLog = NaN;
                end
                globRows(end+1, :) = {subjStr, TASK_KEY, blockNum, condStr, timeHalfLabel, bandName, globalAbs, globalLog, globalRel}; %#ok<AGROW>
            end
        end
    end

    ElectrodeT = local_build_table(elecRows, [1 2 4 5 6 7], [3 8 9 10], elecCols);
    GlobalT = local_build_table(globRows, [1 2 4 5 6], [3 7 8 9], globCols);
    RecordingHalfWindowCountsT = local_build_table(countRows, [1 2 4], [3 5 6 7 8], countCols);

    fprintf(['Task=%s: %d recordings excluded (< %d valid windows), %d reordered (windows not already ' ...
        'chronological). ElectrodeT: %d rows. GlobalT: %d rows.\n'], ...
        TASK_LABEL, nExcludedTooFewWindows, MIN_VALID_WINDOWS_PER_RECORDING, nReordered, height(ElectrodeT), height(GlobalT));

    if isempty(GlobalT)
        warning('Task=%s: no recordings survived the minimum-valid-windows filter -- skipping this Task entirely.', TASK_LABEL);
        continue;
    end

    nExcludedRowsSensitivityElec = 0;
    nExcludedRowsSensitivityGlobal = 0;
    if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
        subjNumForExclElec = local_subject_numeric(ElectrodeT.Subject);
        excludeMaskElec = ismember(subjNumForExclElec, EXCLUDE_NOISY_SUBJECTS);
        nExcludedRowsSensitivityElec = sum(excludeMaskElec);
        ElectrodeT = ElectrodeT(~excludeMaskElec, :);

        subjNumForExclGlobal = local_subject_numeric(GlobalT.Subject);
        excludeMaskGlobal = ismember(subjNumForExclGlobal, EXCLUDE_NOISY_SUBJECTS);
        nExcludedRowsSensitivityGlobal = sum(excludeMaskGlobal);
        GlobalT = GlobalT(~excludeMaskGlobal, :);

        fprintf('Task=%s: RUN_MODE=%s -- excluded %d Electrode-level row(s) and %d Global-level row(s) belonging to subject(s) {%s}; %d Global-level rows remain.\n', ...
            TASK_LABEL, RUN_MODE_LABEL, nExcludedRowsSensitivityElec, nExcludedRowsSensitivityGlobal, ...
            strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','), height(GlobalT));

        if isempty(GlobalT)
            warning('Task=%s: no recordings remain after the sensitivity exclusion -- skipping this Task entirely.', TASK_LABEL);
            continue;
        end
    end
    % ---- END NEW ----

    %% Section 4: data-driven design check (reported, never
    %  hard-coded against a literal -- same "report don't guess" style as
    %  every other script in this pipeline)
    comboT = unique(GlobalT(:, {'Subject','BlockNumber','Condition'}), 'rows');
    nRecordingsIncluded = height(comboT);
    subjNumeric = local_subject_numeric(comboT.Subject);
    uniqSubj = unique(subjNumeric);
    nSubjects = numel(uniqSubj);

    condsPerSubjCounts = nan(nSubjects, 1);
    for i = 1:nSubjects
        rowsI = subjNumeric == uniqSubj(i);
        condsPerSubjCounts(i) = numel(unique(cellstr(comboT.Condition(rowsI))));
    end
    allSubjectsComplete3 = all(condsPerSubjCounts >= 3);
    if allSubjectsComplete3
        integrityIssues{end+1} = sprintf(['Task=%s: every included subject has recordings for all 3 conditions -- ' ...
            'the incomplete-design justification for using an LMM no longer applies. STOP and re-examine.'], TASK_LABEL); %#ok<AGROW>
    end

    expectedGlobalRows = nRecordingsIncluded * nBands * 2;
    if height(GlobalT) ~= expectedGlobalRows
        integrityIssues{end+1} = sprintf(['Task=%s: GlobalT has %d rows, expected %d (%d recordings x %d bands x 2 ' ...
            'TimeHalf levels) -- a recording is missing a TimeHalf or Band row.'], ...
            TASK_LABEL, height(GlobalT), expectedGlobalRows, nRecordingsIncluded, nBands); %#ok<AGROW>
    end

    fprintf('Section 4: %d recordings included, %d subjects, design-completeness flag=%d.\n', ...
        nRecordingsIncluded, nSubjects, allSubjectsComplete3);

    if allSubjectsComplete3
        error('Phase17:DesignAssumptionViolated', ...
            'Task=%s: design verification found a complete 3-condition-per-subject structure. Stopping before fitting any model.', TASK_LABEL);
    end

    %% ====================================================================
    PARAMS17 = struct();
    PARAMS17.ScriptName = 'X17_half_split_condition_time_intervention.m';
    PARAMS17.Task = TASK_LABEL;
    PARAMS17.GeneratedOn = datestr(now); %#ok<TNOW1,DATST>
    PARAMS17.Purpose = ['Explicitly exploratory analysis testing whether the ' ...
        'Condition effect on band power changes across the first vs. second half of a recording (fatigue/practice ' ...
        'question), for the intervention recordings (E / N / R).'];
    PARAMS17.AnalysisLevels = {'Electrode', 'Global'};
    PARAMS17.ROILevelOmitted = ['ROI level is NOT computed by this script -- ROI granularity was not required ' ...
        'for this specific analysis.'];
    PARAMS17.MinValidWindowsPerRecording = MIN_VALID_WINDOWS_PER_RECORDING;
    PARAMS17.MinValidWindowsRationale = ['A recording with fewer than 4 valid windows is excluded entirely (both ' ...
        'TimeHalf halves, at every level). The mean-based half aggregation itself would technically work with 1 ' ...
        'window per half, but a 2-windows-per-half minimum is a reasonable, disclosed floor for the SD/variance-based ' ...
        'diagnostics used elsewhere in this pipeline; disclosed judgment call, see header item 1.'];
    PARAMS17.HalfSplitRule = ['Valid windows ordered chronologically by WindowStartSec ascending (re-sorted ' ...
        'defensively, reordering logged if needed); TimeHalf=1 = first ceil(n/2) valid windows, TimeHalf=2 = ' ...
        'remaining floor(n/2) valid windows, n = total valid windows for that recording.'];
    PARAMS17.AggregationFormula = ['IDENTICAL to X04''s formula, applied once per half: arithmetic mean across that ' ...
        'half''s valid windows, per band per channel; LogAbsolutePower = ln(mean(AbsolutePower)) [never ' ...
        'mean(ln(...))]; RelativePower = mean of window-level relative power; Global = arithmetic mean across ' ...
        'channels of the per-channel means.'];
    PARAMS17.Model = 'LinearMixedEffects_RandomIntercept_FullInteraction';
    PARAMS17.ModelFormula = 'Metric ~ Condition*TimeHalf + (1|Subject)';
    PARAMS17.FixedEffects = 'Condition (3-level: E,N,R; ref=E) * TimeHalf (2-level: 1,2; ref=1), full interaction';
    PARAMS17.RandomEffects = '(1|Subject)';
    PARAMS17.FitMethod = 'REML';
    PARAMS17.DFMethod = 'Satterthwaite';
    PARAMS17.Metrics = METRICS;
    PARAMS17.InteractionFirstRationale = ['One model per Band x Metric x Unit yields THREE separate omnibus F-tests ' ...
        '(Condition main, TimeHalf main, Condition:TimeHalf interaction) from a single anova(lme) call, ' ...
        'specifically to avoid diluting interpretability with a single mixed FDR family ' ...
        'across term-types.'];
    PARAMS17.BandTiers = struct('Primary', {PRIMARY_TIER_BANDS}, 'Secondary', {SECONDARY_TIER_BANDS}, 'Exploratory', {EXPLORATORY_TIER_BANDS}, 'ExploratoryBeta', {EXPLORATORY_BETA_TIER_BANDS});
    PARAMS17.FDRAlpha = FDR_ALPHA;
    PARAMS17.FDRFamilyScheme = ['THREE SEPARATE families per level (Global, Electrode), one per term-type ' ...
        '(Condition-main, TimeHalf-main, Interaction), never pooled with each other and never pooled across Metric. ' ...
        'Global: BH-FDR per (Term x Tier x Metric), pooling across bands within tier (Primary={Alpha}, size 1, ' ...
        'reported as-is/no correction; Secondary={Delta,Theta,Beta}, size 3; Exploratory={5 bands}, size 5) -- ' ...
        'IDENTICAL tiered structure to X05''s own Global scheme (verified by re-reading X05''s FDR loop directly), ' ...
        'just computed three times (once per term-type) instead of once. Electrode: BH-FDR per (Term x individual ' ...
        'Band x Metric), pooling across up to 62 channels for that single band -- IDENTICAL structure to X09''s own ' ...
        'Electrode scheme (verified by re-reading X09''s FDR loop directly); Tier is a label column only at this ' ...
        'level, never a family boundary. Electrode family size = nChannels (per band); Global family size = ' ...
        'nBandsInTier (1, 3, or 5).'];
    PARAMS17.PairwiseContrastMethod = ['Computed ONLY for (Unit x Metric x Band x Term) combinations whose own term-' ...
        'type FDR family survives BH-FDR, independently per term-type (a unit surviving more than one term-type ' ...
        'family gets ALL applicable contrast sets below): Condition-main survival -> marginal E-N/E-R/N-R contrasts ' ...
        '(averaged over TimeHalf), Holm-corrected, family size=3. Interaction survival -> simple-effects E-N/E-R/N-R ' ...
        'contrasts computed separately within TimeHalf=1 and TimeHalf=2 (6 contrasts total), Holm-corrected, family ' ...
        'size=6. TimeHalf-main survival -> a single marginal TimeHalf2-vs-TimeHalf1 contrast (averaged over ' ...
        'Condition), NO Holm correction (family size=1, a single comparison needs no multiplicity correction).'];
    PARAMS17.BaselineExclusionNote = sprintf('%d Baseline recording(s) excluded before any modeling for Task=%s.', nBaselineExcluded, TASK_LABEL);
    PARAMS17.JudgmentCalls = { ...
        sprintf(['Minimum valid windows per recording = %d (see header item 1 / MinValidWindowsRationale above); ' ...
        '%d recording(s) excluded under this rule for Task=%s.'], MIN_VALID_WINDOWS_PER_RECORDING, nExcludedTooFewWindows, TASK_LABEL), ...
        'ROI level is not computed by this script (not required for this analysis).', ...
        ['FDR family construction mirrors this pipeline''s TWO existing, level-specific conventions exactly (verified ' ...
         'by re-reading X05''s and X09''s actual FDR loops directly, not their headers): Global pools p-values across ' ...
         'bands WITHIN a tier (Primary/Secondary/Exploratory), per term-type x Metric, matching X05; Electrode pools ' ...
         'p-values across channels for one individual band at a time (Tier is a label only, not a family boundary), ' ...
         'per term-type x Metric, matching X09. See header item D / PARAMS17.FDRFamilyScheme.'], ...
        ['The marginal-Condition-contrast-across-TimeHalf arithmetic (Section 8) is NOT the naive "Condition-' ...
         'coefficient-only" contrast -- it includes half the Condition:TimeHalf interaction coefficient, because the ' ...
         'model includes that interaction term. This was derived and verified against the model''s own ' ...
         'CoefficientNames at runtime, not assumed; see the extended comment above Section 8.']};
    PARAMS17.SourcePhase2PSDFile = psdMatFile;
    PARAMS17.DesignVerification = struct('nRecordingsIncluded', nRecordingsIncluded, 'nSubjects', nSubjects, ...
        'nExcludedTooFewWindows', nExcludedTooFewWindows, 'nReordered', nReordered, ...
        'nNonPositiveLogFlagged', nNonPositiveLogFlagged, 'IntegrityIssues', {integrityIssues});
    PARAMS17.RunDateTime = datestr(now); %#ok<TNOW1,DATST>
    PARAMS17.RunMode = RUN_MODE;
    if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
        PARAMS17.ExcludedSubjects = EXCLUDE_NOISY_SUBJECTS;
        PARAMS17.ExcludedSubjectsRationale = EXCLUSION_RATIONALE_intervention;
        PARAMS17.RowsExcludedForSensitivity = struct( ...
            'ElectrodeLevel', nExcludedRowsSensitivityElec, ...
            'GlobalLevel', nExcludedRowsSensitivityGlobal);
    else
        PARAMS17.ExcludedSubjects = [];
    end

    save(fullfile(STATS_DIR, 'phase17_parameters.mat'), 'PARAMS17');
    save(fullfile(STATS_DIR, 'half_split_data.mat'), 'ElectrodeT', 'GlobalT', 'RecordingHalfWindowCountsT');

    %% ====================================================================
    %  SECTION 6: ONE-MODEL, THREE-TERM OMNIBUS FITS
    %  (Electrode: per Channel x Band x Metric. Global: per Band x Metric.)
    % ====================================================================
    TERM_NAMES = {'Condition', 'TimeHalf', 'Condition:TimeHalf'};
    TERM_KEYS  = {'ConditionMain', 'TimeHalfMain', 'Interaction'};   % used as output column values

    levelConfigs = struct('Name', {'Electrode', 'Global'}, ...
        'Table', {ElectrodeT, GlobalT}, ...
        'UnitCol', {'Channel', ''}, ...
        'Units', {ChannelLabelsReference, {'Global'}});

    omnibusRows = cell(0, 9);
    omnibusCols = {'Level','Unit','Metric','Band','Term','FStat','DF1','DF2','pValue'};
    lmeMap = containers.Map();   % key: 'Level|Unit|Metric|Band' -> fitted lme
    coeffNamesLogged = false;
    fitFailCount = 0;

    for lc = 1:numel(levelConfigs)
        levelName = levelConfigs(lc).Name;
        Tlevel = levelConfigs(lc).Table;
        unitCol = levelConfigs(lc).UnitCol;
        units = levelConfigs(lc).Units;
        for u = 1:numel(units)
            unitName = units{u};
            if strcmp(levelName, 'Electrode')
                unitMask = strcmp(Tlevel.(unitCol), unitName);
            else
                unitMask = true(height(Tlevel), 1);
            end
            for mI = 1:numel(METRICS)
                metricName = METRICS{mI};
                for b = 1:nBands
                    bandName = ALL_BANDS{b};
                    subT = Tlevel(unitMask & strcmp(Tlevel.Band, bandName), :);
                    fr = local_fit_two_factor_omnibus(subT, metricName, CONDITION_ORDER, TIMEHALF_ORDER);
                    if ~fr.ok
                        fitFailCount = fitFailCount + 1;
                        integrityIssues{end+1} = sprintf('%s %s: fitlme FAILED for %s / %s (Task=%s): %s', ...
                            levelName, unitName, metricName, bandName, TASK_LABEL, fr.errMsg); %#ok<AGROW>
                        continue;
                    end
                    if ~coeffNamesLogged
                        fprintf('Verified model CoefficientNames (first successful fit, Task=%s): %s\n', ...
                            TASK_LABEL, strjoin(fr.lme.CoefficientNames, ', '));
                        coeffNamesLogged = true;
                    end
                    key = sprintf('%s|%s|%s|%s', levelName, unitName, metricName, bandName);
                    lmeMap(key) = fr.lme;
                    for tI = 1:numel(TERM_NAMES)
                        omnibusRows(end+1, :) = {levelName, unitName, metricName, bandName, TERM_KEYS{tI}, ...
                            fr.F(tI), fr.DF1(tI), fr.DF2(tI), fr.p(tI)}; %#ok<AGROW>
                    end
                end
            end
        end
        fprintf('Section 6 (%s level, Task=%s): %d units x %d metrics x %d bands fit attempted; %d failure(s) so far.\n', ...
            levelName, TASK_LABEL, numel(units), numel(METRICS), nBands, fitFailCount);
    end

    OmnibusT = local_build_table(omnibusRows, [1 2 3 4 5], [6 7 8 9], omnibusCols);
    OmnibusT.Tier = cellfun(@(b) tierOf(b), OmnibusT.Band, 'UniformOutput', false);   % label at both levels; a family boundary ONLY at Global (see Section 7)
    fprintf('Section 6 complete (Task=%s): %d omnibus rows (3 terms x Electrode+Global units).\n', TASK_LABEL, height(OmnibusT));

    %% ====================================================================
    OmnibusT.pFDR = nan(height(OmnibusT), 1);
    OmnibusT.SignificantFDR = false(height(OmnibusT), 1);
    OmnibusT.FamilySize = nan(height(OmnibusT), 1);

    % --- GLOBAL: pool across bands WITHIN a tier, per Term x Metric ---
    for tI = 1:numel(TERM_KEYS)
        termKey = TERM_KEYS{tI};
        for t = 1:size(tierDefs, 1)
            tierBands = tierDefs{t, 2};
            for mI = 1:numel(METRICS)
                metricName = METRICS{mI};
                idx = strcmp(OmnibusT.Level, 'Global') & strcmp(OmnibusT.Term, termKey) & ...
                    strcmp(OmnibusT.Metric, metricName) & ismember(OmnibusT.Band, tierBands);
                if ~any(idx); continue; end
                pvals = OmnibusT.pValue(idx);
                OmnibusT.FamilySize(idx) = numel(pvals);
                [adjP, reject] = local_bh_fdr_or_raw(pvals, FDR_ALPHA);
                OmnibusT.pFDR(idx) = adjP;
                OmnibusT.SignificantFDR(idx) = reject;
            end
        end
    end

    % --- ELECTRODE: pool across channels for ONE individual band, per
    %     Term x Band x Metric (Tier NOT a family boundary here) ---
    for tI = 1:numel(TERM_KEYS)
        termKey = TERM_KEYS{tI};
        for b = 1:nBands
            bandName = ALL_BANDS{b};
            for mI = 1:numel(METRICS)
                metricName = METRICS{mI};
                idx = strcmp(OmnibusT.Level, 'Electrode') & strcmp(OmnibusT.Term, termKey) & ...
                    strcmp(OmnibusT.Metric, metricName) & strcmp(OmnibusT.Band, bandName);
                if ~any(idx); continue; end
                pvals = OmnibusT.pValue(idx);
                OmnibusT.FamilySize(idx) = numel(pvals);
                [adjP, reject] = local_bh_fdr_or_raw(pvals, FDR_ALPHA);
                OmnibusT.pFDR(idx) = adjP;
                OmnibusT.SignificantFDR(idx) = reject;
            end
        end
    end

    OmnibusCondition_T   = OmnibusT(strcmp(OmnibusT.Term, 'ConditionMain'), :);
    OmnibusTimeHalf_T    = OmnibusT(strcmp(OmnibusT.Term, 'TimeHalfMain'), :);
    OmnibusInteraction_T = OmnibusT(strcmp(OmnibusT.Term, 'Interaction'), :);

    fprintf(['Section 7 complete (Task=%s): SignificantFDR counts -- ConditionMain=%d/%d, TimeHalfMain=%d/%d, ' ...
        'Interaction=%d/%d.\n'], TASK_LABEL, sum(OmnibusCondition_T.SignificantFDR), height(OmnibusCondition_T), ...
        sum(OmnibusTimeHalf_T.SignificantFDR), height(OmnibusTimeHalf_T), ...
        sum(OmnibusInteraction_T.SignificantFDR), height(OmnibusInteraction_T));

    %% ====================================================================
    %  SECTION 8: PAIRWISE / SIMPLE-EFFECTS CONTRASTS
    %
    %  Under treatment (dummy) coding with Condition ref=E and TimeHalf
    %  ref='1', the fitted value at (Condition=c, TimeHalf=h) is:
    %     mu(c,h) = Intercept + Condition_c [0 if c=E] + TimeHalf_2*[h=2]
    %               + (Condition_c:TimeHalf_2)*[c~=E]*[h=2]
    %
    %  MARGINAL TIMEHALF CONTRAST (averaged equally over the 3 Condition
    %  levels -- used when the TIMEHALF MAIN EFFECT survives FDR):
    %     Marginal(h) = [mu(E,h) + mu(N,h) + mu(R,h)] / 3
    %     Marginal(2) - Marginal(1) = TimeHalf_2 + (1/3)*(Condition_N:TimeHalf_2 + Condition_R:TimeHalf_2)
    %  (E contributes 0 to the interaction sum because E is the Condition
    %  reference level and has no interaction dummy.)
    % ====================================================================
    pairwiseRows = cell(0, 14);
    pairwiseCols = {'Level','Unit','Metric','Band','SurvivingTerm','Contrast','Estimate','SE','DF','tStat','pValue','pHolm','CI_Low','CI_High'};

    % --- Condition main effect survivors: marginal contrasts, Holm(3) ---
    for r = find(OmnibusCondition_T.SignificantFDR)'
        levelName = OmnibusCondition_T.Level{r}; unitName = OmnibusCondition_T.Unit{r};
        metricName = OmnibusCondition_T.Metric{r}; bandName = OmnibusCondition_T.Band{r};
        key = sprintf('%s|%s|%s|%s', levelName, unitName, metricName, bandName);
        if ~isKey(lmeMap, key); continue; end
        lme = lmeMap(key);
        rows3 = local_run_pairwise_condition_marginal(lme, REF_COND, NONREF_TIMEHALF, CONTRASTS, FDR_ALPHA);
        for c = 1:3
            pairwiseRows(end+1, :) = [{levelName, unitName, metricName, bandName, 'ConditionMain'}, rows3(c, :)]; %#ok<AGROW>
        end
    end

    % --- Interaction survivors: simple effects within each TimeHalf, Holm(6) ---
    for r = find(OmnibusInteraction_T.SignificantFDR)'
        levelName = OmnibusInteraction_T.Level{r}; unitName = OmnibusInteraction_T.Unit{r};
        metricName = OmnibusInteraction_T.Metric{r}; bandName = OmnibusInteraction_T.Band{r};
        key = sprintf('%s|%s|%s|%s', levelName, unitName, metricName, bandName);
        if ~isKey(lmeMap, key); continue; end
        lme = lmeMap(key);
        rows6 = local_run_pairwise_condition_simple(lme, REF_COND, REF_TIMEHALF, NONREF_TIMEHALF, TIMEHALF_ORDER, CONTRASTS, FDR_ALPHA);
        for c = 1:6
            pairwiseRows(end+1, :) = [{levelName, unitName, metricName, bandName, 'Interaction'}, rows6(c, :)]; %#ok<AGROW>
        end
    end

    % --- TimeHalf main effect survivors: single marginal contrast, NO Holm ---
    for r = find(OmnibusTimeHalf_T.SignificantFDR)'
        levelName = OmnibusTimeHalf_T.Level{r}; unitName = OmnibusTimeHalf_T.Unit{r};
        metricName = OmnibusTimeHalf_T.Metric{r}; bandName = OmnibusTimeHalf_T.Band{r};
        key = sprintf('%s|%s|%s|%s', levelName, unitName, metricName, bandName);
        if ~isKey(lmeMap, key); continue; end
        lme = lmeMap(key);
        row1 = local_run_pairwise_timehalf_marginal(lme, REF_COND, CONDITION_ORDER, REF_TIMEHALF, NONREF_TIMEHALF);
        pairwiseRows(end+1, :) = [{levelName, unitName, metricName, bandName, 'TimeHalfMain'}, row1]; %#ok<AGROW>
    end

    PairwiseT = local_build_table(pairwiseRows, [1 2 3 4 5 6], [7 8 9 10 11 12 13 14], pairwiseCols);
    fprintf('Section 8 complete (Task=%s): %d pairwise/simple-effects contrast rows.\n', TASK_LABEL, height(PairwiseT));

    %% ====================================================================
    %  SECTION 9: SAVE OUTPUTS
    % ====================================================================
    save(fullfile(STATS_DIR, 'phase17_omnibus_results.mat'), 'OmnibusT', 'OmnibusCondition_T', 'OmnibusTimeHalf_T', 'OmnibusInteraction_T');
    save(fullfile(STATS_DIR, 'phase17_pairwise_contrasts.mat'), 'PairwiseT');
    writetable(OmnibusT, fullfile(STATS_DIR, 'phase17_omnibus_results.csv'));
    writetable(OmnibusCondition_T, fullfile(STATS_DIR, 'phase17_fdr_family_condition_main.csv'));
    writetable(OmnibusTimeHalf_T, fullfile(STATS_DIR, 'phase17_fdr_family_timehalf_main.csv'));
    writetable(OmnibusInteraction_T, fullfile(STATS_DIR, 'phase17_fdr_family_interaction.csv'));
    if ~isempty(PairwiseT)
        writetable(PairwiseT, fullfile(STATS_DIR, 'phase17_pairwise_contrasts.csv'));
    end
    writetable(RecordingHalfWindowCountsT, fullfile(STATS_DIR, 'phase17_recording_half_window_counts.csv'));

    fprintf('\nSaved output to: %s\n', STATS_DIR);

    %% ====================================================================
    %  SECTION 10: REPORT
    % ====================================================================
    reportLines = {};
    reportLines{end+1} = repmat('=', 1, 60);
    reportLines{end+1} = sprintf('PHASE 17 (X17) -- CONDITION x TIME-HALF EXPLORATORY ANALYSIS REPORT (Task=%s)', TASK_LABEL);
    reportLines{end+1} = repmat('=', 1, 60);
    reportLines{end+1} = '';
    reportLines{end+1} = 'Explicitly exploratory analysis (fatigue/practice';
    reportLines{end+1} = 'question). Full is out of scope for this half-split analysis.';
    reportLines{end+1} = '';
    reportLines{end+1} = sprintf('Baseline recordings excluded before modeling: %d', nBaselineExcluded);
    reportLines{end+1} = sprintf('Recordings excluded for < %d valid windows (Judgment Call 1): %d', MIN_VALID_WINDOWS_PER_RECORDING, nExcludedTooFewWindows);
    reportLines{end+1} = sprintf('Recordings with windows defensively reordered (not already chronological): %d', nReordered);
    reportLines{end+1} = sprintf('Recordings included in this analysis: %d | Subjects: %d', nRecordingsIncluded, nSubjects);
    reportLines{end+1} = sprintf('Non-finite LogAbsolutePower values flagged (non-positive AbsolutePower denominator): %d', nNonPositiveLogFlagged);
    reportLines{end+1} = '';
    reportLines{end+1} = 'Model: Metric ~ Condition*TimeHalf + (1|Subject), REML, Satterthwaite DF.';
    reportLines{end+1} = 'Analysis levels: Electrode (per channel) and Global (whole-scalp). ROI level NOT computed (Judgment Call 2).';
    reportLines{end+1} = 'Metrics: LogAbsolutePower and RelativePower, co-primary, parallel, never pooled for FDR.';
    reportLines{end+1} = '';
    reportLines{end+1} = 'FDR family scheme: THREE SEPARATE families per level (Condition-main, TimeHalf-main,';
    reportLines{end+1} = 'Interaction), never pooled across Metric. Global: BH-FDR per (Term x Tier x Metric),';
    reportLines{end+1} = 'pooling across bands WITHIN a tier (Primary=Alpha size1/no-correction, Secondary=3 bands,';
    reportLines{end+1} = 'Exploratory=5 bands) -- matching X05''s own Global scheme exactly. Electrode: BH-FDR per';
    reportLines{end+1} = '(Term x individual Band x Metric), pooling across up to 62 channels for that one band --';
    reportLines{end+1} = 'matching X09''s own Electrode scheme exactly (Tier is a label only there, never a family';
    reportLines{end+1} = 'boundary).';
    reportLines{end+1} = '';
    reportLines{end+1} = '--- FDR-surviving omnibus results, per term-type family ---';
    reportLines{end+1} = sprintf('Condition-main family : %d / %d survive FDR (alpha=%.2f)', sum(OmnibusCondition_T.SignificantFDR), height(OmnibusCondition_T), FDR_ALPHA);
    reportLines{end+1} = sprintf('TimeHalf-main family   : %d / %d survive FDR (alpha=%.2f)', sum(OmnibusTimeHalf_T.SignificantFDR), height(OmnibusTimeHalf_T), FDR_ALPHA);
    reportLines{end+1} = sprintf('Interaction family     : %d / %d survive FDR (alpha=%.2f)', sum(OmnibusInteraction_T.SignificantFDR), height(OmnibusInteraction_T), FDR_ALPHA);
    reportLines{end+1} = '';
    survivingCond = OmnibusCondition_T(OmnibusCondition_T.SignificantFDR, :);
    reportLines{end+1} = sprintf('--- Condition-main FDR survivors (%d) ---', height(survivingCond));
    for r = 1:height(survivingCond)
        reportLines{end+1} = sprintf('  [%s] Level=%-9s Unit=%-6s %-25s F(%g,%.1f)=%.3f p=%.4f pFDR=%.4f family_n=%d', ...
            survivingCond.Metric{r}, survivingCond.Level{r}, survivingCond.Unit{r}, survivingCond.Band{r}, ...
            survivingCond.DF1(r), survivingCond.DF2(r), survivingCond.FStat(r), survivingCond.pValue(r), ...
            survivingCond.pFDR(r), survivingCond.FamilySize(r)); %#ok<AGROW>
    end
    reportLines{end+1} = '';
    survivingTH = OmnibusTimeHalf_T(OmnibusTimeHalf_T.SignificantFDR, :);
    reportLines{end+1} = sprintf('--- TimeHalf-main FDR survivors (%d) ---', height(survivingTH));
    for r = 1:height(survivingTH)
        reportLines{end+1} = sprintf('  [%s] Level=%-9s Unit=%-6s %-25s F(%g,%.1f)=%.3f p=%.4f pFDR=%.4f family_n=%d', ...
            survivingTH.Metric{r}, survivingTH.Level{r}, survivingTH.Unit{r}, survivingTH.Band{r}, ...
            survivingTH.DF1(r), survivingTH.DF2(r), survivingTH.FStat(r), survivingTH.pValue(r), ...
            survivingTH.pFDR(r), survivingTH.FamilySize(r)); %#ok<AGROW>
    end
    reportLines{end+1} = '';
    survivingInt = OmnibusInteraction_T(OmnibusInteraction_T.SignificantFDR, :);
    reportLines{end+1} = sprintf('--- Interaction FDR survivors (%d) ---', height(survivingInt));
    for r = 1:height(survivingInt)
        reportLines{end+1} = sprintf('  [%s] Level=%-9s Unit=%-6s %-25s F(%g,%.1f)=%.3f p=%.4f pFDR=%.4f family_n=%d', ...
            survivingInt.Metric{r}, survivingInt.Level{r}, survivingInt.Unit{r}, survivingInt.Band{r}, ...
            survivingInt.DF1(r), survivingInt.DF2(r), survivingInt.FStat(r), survivingInt.pValue(r), ...
            survivingInt.pFDR(r), survivingInt.FamilySize(r)); %#ok<AGROW>
    end
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Pairwise / simple-effects contrasts (gated on each term-type''s own FDR survival) ---';
    if isempty(PairwiseT)
        reportLines{end+1} = '  (none -- no Unit x Metric x Band x Term combination survived its own FDR family)';
    else
        for r = 1:height(PairwiseT)
            reportLines{end+1} = sprintf('  [%s|%s] Level=%-9s Unit=%-6s %-25s Term=%-13s %-8s est=%.4f (95%% CI %.4f to %.4f) t=%.3f df=%.1f p=%.4f pHolm=%.4f', ...
                PairwiseT.Metric{r}, PairwiseT.SurvivingTerm{r}, PairwiseT.Level{r}, PairwiseT.Unit{r}, PairwiseT.Band{r}, ...
                PairwiseT.SurvivingTerm{r}, PairwiseT.Contrast{r}, PairwiseT.Estimate(r), PairwiseT.CI_Low(r), ...
                PairwiseT.CI_High(r), PairwiseT.tStat(r), PairwiseT.DF(r), PairwiseT.pValue(r), PairwiseT.pHolm(r)); %#ok<AGROW>
        end
    end
    reportLines{end+1} = '';
    reportLines{end+1} = 'Holm family sizes used for pairwise contrasts (disclosed, differ by surviving term-type):';
    reportLines{end+1} = '  ConditionMain survivors : 3 contrasts (E-N, E-R, N-R, averaged over TimeHalf), Holm family size = 3.';
    reportLines{end+1} = '  Interaction survivors   : 6 contrasts (E-N, E-R, N-R within EACH TimeHalf level), Holm family size = 6.';
    reportLines{end+1} = '  TimeHalfMain survivors  : 1 contrast (TimeHalf2 - TimeHalf1, averaged over Condition), NO Holm correction (family size = 1).';
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Judgment calls (disclosed) ---';
    for i = 1:numel(PARAMS17.JudgmentCalls)
        reportLines{end+1} = ['  - ' PARAMS17.JudgmentCalls{i}]; %#ok<AGROW>
    end
    reportLines{end+1} = '';
    if ~isempty(integrityIssues)
        reportLines{end+1} = sprintf('--- INTEGRITY ISSUES (%d, reported, not auto-corrected) ---', numel(integrityIssues));
        for i = 1:numel(integrityIssues)
            reportLines{end+1} = ['  - ' integrityIssues{i}]; %#ok<AGROW>
        end
    else
        reportLines{end+1} = 'Integrity issues: none.';
    end
    reportLines{end+1} = '';
    reportLines{end+1} = 'Output files written to:';
    reportLines{end+1} = ['  ' STATS_DIR];
    reportLines{end+1} = repmat('=', 1, 60);

    reportText = strjoin(reportLines, '\n');
    fprintf('%s\n', reportText);
    fid = fopen(fullfile(REPORT_DIR, sprintf('X17_half_split_report_%s_%s.txt', TASK_LABEL, RUN_MODE_LABEL)), 'w');
    fprintf(fid, '%s\n', reportText);
    fclose(fid);

end % TASK_TYPES loop

fprintf('\n=== X17 complete. This is a purely additive exploratory analysis -- no other script''s output is touched. ===\n');

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
            error('Phase17:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name1);
        end
        cRow(idx1) = cRow(idx1) - 1;
    end
    if ~strcmp(level2, refLevel)
        name2 = sprintf('%s_%s', varName, level2);
        idx2 = find(strcmp(coeffNames, name2));
        if isempty(idx2)
            error('Phase17:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name2);
        end
        cRow(idx2) = cRow(idx2) + 1;
    end
end

function T = local_build_table(rows, textColIdx, numColIdx, varNames)
    nCols = numel(varNames);
    if size(rows, 2) ~= nCols
        error('Phase17:BuildTableColumnMismatch', ...
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
            error('Phase17:BuildTableColumnMismatch', 'Column %d not classified as text or numeric.', j);
        end
    end
    T = table(cols{:}, 'VariableNames', varNames);
end

function [adjP, reject] = local_bh_fdr_or_raw(pvals, alphaFDR)
    if numel(pvals) == 1
        adjP = pvals;
        reject = pvals <= alphaFDR;
    else
        [adjP, reject] = local_fdr_bh(pvals, alphaFDR);
    end
end

function fr = local_fit_two_factor_omnibus(subT, metricName, conditionOrder, timeHalfOrder)
% Fits Metric ~ Condition*TimeHalf + (1|Subject), REML, and extracts the
% THREE omnibus Satterthwaite F-tests (Condition, TimeHalf,
% Condition:TimeHalf) from ONE anova(lme) call. Returns a struct rather
% than erroring out, so the calling loop can log a failure and continue,
% exactly like X05/X11's own try/catch pattern around fitlme. fr.F/DF1/
% DF2/p are each 1x3 vectors, in the order {Condition, TimeHalf,
% Condition:TimeHalf} (matching TERM_NAMES in the caller).
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
    name1 = sprintf('%s:%s', condDummyName, timeHalfDummyName);
    name2 = sprintf('%s:%s', timeHalfDummyName, condDummyName);
    idx = find(strcmp(coeffNames, name1), 1);
    if isempty(idx)
        idx = find(strcmp(coeffNames, name2), 1);
    end
    if isempty(idx)
        error('Phase17:InteractionCoeffNotFound', ...
            ['Could not find interaction coefficient "%s" or "%s" among the model''s CoefficientNames: %s. ' ...
             'The marginal/simple-effects contrast arithmetic in this script assumes one of these two exact ' ...
             'naming orders; if MATLAB''s actual naming differs, this must be fixed here before trusting any ' ...
             'contrast output.'], name1, name2, strjoin(coeffNames, ', '));
    end
end

function cRow = local_marginal_condition_contrast(coeffNames, refCond, nonRefTimeHalf, level1, level2)
% Marginal Condition contrast (level2 - level1), averaged EQUALLY over
% TimeHalf=1 and TimeHalf=2. See the extended derivation comment above
% Section 8 in the main script: this is the Condition-coefficient
% contrast PLUS HALF the corresponding Condition:TimeHalf interaction
% coefficient contrast -- NOT the Condition coefficient alone.
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
    cRow = local_condition_contrast(coeffNames, 'Condition', refCond, level1, level2);
    if strcmp(timeHalfLevel, refTimeHalf)
        return;   % reference TimeHalf level: no interaction term contributes
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
% Marginal TimeHalf contrast (TimeHalf=nonRefTimeHalf minus
% TimeHalf=refTimeHalf), averaged EQUALLY over all levels of Condition
% (including the reference level, which contributes 0 to the interaction
% sum since it has no interaction dummy). See the extended derivation
% comment above Section 8 in the main script: TimeHalf_2 coefficient PLUS
% (1/nCondLevels) times the sum of every non-reference Condition's
% interaction coefficient.
    cRow = zeros(1, numel(coeffNames));
    timeHalfDummyName = sprintf('TimeHalf_%s', nonRefTimeHalf);
    idxTH = find(strcmp(coeffNames, timeHalfDummyName), 1);
    if isempty(idxTH)
        error('Phase17:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', timeHalfDummyName);
    end
    cRow(idxTH) = cRow(idxTH) + 1;
    nCondLevels = numel(condLevels);
    for k = 1:nCondLevels
        lvl = condLevels{k};
        if strcmp(lvl, refCond); continue; end   % reference level: no interaction dummy, contributes 0
        condDummyName = sprintf('Condition_%s', lvl);
        idxInt = local_find_interaction_coeff_idx(coeffNames, condDummyName, timeHalfDummyName);
        cRow(idxInt) = cRow(idxInt) + (1 / nCondLevels);
    end
end

function rows3 = local_run_pairwise_condition_marginal(lme, refCond, nonRefTimeHalf, CONTRASTS, alphaFDR)
% Marginal Condition contrasts (E-N, E-R, N-R), averaged over TimeHalf,
% Holm-corrected (family size = 3). Returns a 3x9 cell array of
% {Contrast, Estimate, SE, DF, tStat, pValue, pHolm, CI_Low, CI_High}.
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
% Simple-effects Condition contrasts (E-N, E-R, N-R), computed SEPARATELY
% within EACH TimeHalf level (6 contrasts total), Holm-corrected across
% all 6 (family size = 6). Contrast labels are
% suffixed with "@ TimeHalf=<h>" so the level is unambiguous in output.
% Returns a 6x9 cell array of {Contrast, Estimate, SE, DF, tStat, pValue,
% pHolm, CI_Low, CI_High}.
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
% Single marginal TimeHalf contrast (TimeHalf2 - TimeHalf1), averaged over
% Condition. NO Holm correction applied (a single comparison needs no
% multiplicity correction) -- pHolm is set equal to the
% raw pValue for schema consistency with the other pairwise tables, not
% because a correction was computed. Returns a 1x9 cell row of {Contrast,
% Estimate, SE, DF, tStat, pValue, pHolm, CI_Low, CI_High}.
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
    row1 = {label, est, se, df2c, tstat, p, p, ciLow, ciHigh};   % pHolm := pValue (no correction; see comment above)
end
