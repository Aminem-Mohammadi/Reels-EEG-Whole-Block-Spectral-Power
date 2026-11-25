%% ========================================================================
% X08_roi_level_statistics_cognitive.m
%
% Purpose: predefined 5-ROI statistical analysis of band power. Runs as
% three separate, parallel analyses, one per Task (Full / N-back / SART);
% Task is never modeled as a factor and Task x Condition is never fit.
% Before ROI aggregation, the channel-level band-power table is filtered
% to Condition in {E, N, R} (Baseline rows are excluded; the excluded
% count is logged per Task).
%
% Method: LMM (Metric ~ Condition + (1|Subject), REML, Satterthwaite DF)
% on ROI-averaged band power. ROI membership (Frontal, Central, Temporal,
% Parietal, Occipital) is loaded read-only from X07's shared output
% (phase5_channel_roi_definitions.mat) and never redefined here; P9/P10
% are excluded from ROI averaging only. ROI AbsolutePower/RelativePower =
% arithmetic mean of channel-level values across the ROI's member
% channels; ROI LogAbsolutePower = log(ROI AbsolutePower) (the ROI's
% absolute power is averaged first, then log-transformed). Outcome tiers
% (Primary = Alpha; Secondary = Delta/Theta/Beta; Exploratory = 3 gamma
% bands) label the outputs only; BH-FDR families are Band x Metric across
% the 5 ROIs (family size 5), independently per Task -- 14 families x 5
% ROIs = 70 omnibus tests per Task, never pooled. The co-primary metrics
% (LogAbsolutePower, RelativePower) are never pooled for FDR. Within an
% ROI x Band x Metric whose omnibus test survives its FDR family, the 3
% pairwise contrasts are Holm-corrected (family size 3). This script's
% own pairwise-contrast label convention is E-N/E-R/N-R (E as minuend),
% which differs from X05/X06's N-E/R-E/R-N (E as subtrahend); see
% CODE_NOTES.md's "Naming differences from the paper" for how this is
% harmonized for display. Effect sizes are model-based estimated
% difference + 95% CI only; no standardized effect size is computed.
% Design-verification counts (recording/condition/pair-type coverage) are
% computed from the data and reported, never asserted against a
% hard-coded literal; the same structural stop-check as X05 applies (if
% every subject has all 3 conditions, the incomplete-design justification
% for an LMM no longer applies and the script stops before fitting).
% Residual diagnostics are produced for the Primary+Secondary bands only.
%
% Requirements: MATLAB with Statistics and Machine Learning Toolbox
% (fitlme, anova(LinearMixedModel), coefTest, qqplot). X04 output must
% exist (phase3_subject_level_bandpower.mat); X07 output must exist
% (phase5_channel_roi_definitions.mat) -- this script errors immediately
% if that file is missing and never reconstructs ROI membership itself.
%
% Inputs: X04's channel-level band-power table, X07's ROI definitions.
%
% Outputs (one complete set per Task per RUN_MODE):
%   <ANALYSIS_ROOT>/<Full|Nback|SART>/08_roi_stats/<RUN_MODE>/
%       phase5a_roi_bandpower.mat            (RoiBandPowerT + PARAMS5A)
%       phase5a_roi_omnibus_results.mat/.xlsx/.csv
%       phase5a_roi_pairwise_contrasts.mat/.xlsx/.csv
%       phase5a_parameters.mat
%       diagnostic_residuals_<ROI>_<band>_<metric>.png  (Primary+Secondary bands only)
%   <ANALYSIS_ROOT>/09_reports/
%       phase5a_roi_level_statistics_report_<Full|Nback|SART>_<RunModeLabel>.txt
%
% RUN_MODE: full_sample / sensitivity_excl_noisy (see
% EXCLUSION_CONFIG_cognitive.m); run this script once per mode. ROI
% membership does not depend on RUN_MODE.
% ========================================================================

clear; clc;

%% ------------------------------------------------------------------
%  0. SETTINGS, PATHS, EXISTENCE CHECKS (shared across all three Tasks)
% ------------------------------------------------------------------

config;
DATA_ROOT_FOR_SIBLING = DATA_DIR;
[DATA_ROOT_PARENT, ~] = fileparts(DATA_ROOT_FOR_SIBLING);
ANALYSIS_ROOT = OUT_DIR;

SUBJ_DIR   = fullfile(ANALYSIS_ROOT, '04_subject_level');
DEF_DIR    = fullfile(ANALYSIS_ROOT, '07_roi_electrode_definitions');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');

bandPowerFile = fullfile(SUBJ_DIR, 'phase3_subject_level_bandpower.mat');
roiDefFile    = fullfile(DEF_DIR, 'phase5_channel_roi_definitions.mat');

if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

% ---- Noisy-subject exclusion list ----
% Loads the single shared exclusion-list config (EXCLUDE_NOISY_SUBJECTS,
% RUN_MODES_cognitive) -- see EXCLUSION_CONFIG_cognitive.m in this same
% folder for the vector itself and each subject's rationale. To add or
% remove a subject from the sensitivity run, edit ONLY that file.
run(fullfile(fileparts(mfilename('fullpath')), 'EXCLUSION_CONFIG_cognitive.m'));

% USER-EDITABLE: run this script ONCE PER VALUE below to get both the
% full-sample and the sensitivity (noisy-subjects-excluded) results,
% saved side by side under separate subfolders -- nothing else in this
% script needs to change between the two runs.
RUN_MODE = 'sensitivity_excl_noisy';  % or: 'sensitivity_excl_noisy'
if ~ismember(RUN_MODE, RUN_MODES_cognitive)
    error('Phase5A:BadRunMode', 'RUN_MODE must be one of: %s', strjoin(RUN_MODES_cognitive, ', '));
end
RUN_MODE_LABEL = RUN_MODES_DISPLAY_cognitive.(RUN_MODE);
fprintf('RUN_MODE = %s\n', RUN_MODE_LABEL);
if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
    fprintf('Excluding subject(s) {%s} from every model in this run -- see EXCLUSION_CONFIG_cognitive.m for rationale.\n', ...
        strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','));
end
% ---- END NEW ----

fprintf('========================================================\n');
fprintf('PHASE 5A -- PREDEFINED 5-ROI STATISTICAL ANALYSIS, COGNITIVE-TASK DATA\n');
fprintf('Running as THREE INDEPENDENT, PARALLEL analyses (Full / Nback / SART).\n');
fprintf('========================================================\n');
fprintf('ANALYSIS_ROOT   : %s\n', ANALYSIS_ROOT);
fprintf('bandPowerFile    : %s\n', bandPowerFile);
fprintf('roiDefFile       : %s\n', roiDefFile);
fprintf('--------------------------------------------------------\n\n');

if ~exist(bandPowerFile, 'file')
    error('Phase5A:MissingInput', ['X04 band-power file not found: %s\n' ...
        'Run X04_aggregate_cognitive_subject_level.m first.'], bandPowerFile);
end
if ~exist(roiDefFile, 'file')
    error('Phase5A:MissingInput', ['X07 ROI definition file not found: %s\n' ...
        'Run X07_roi_electrode_definitions_cognitive.m first. This script will NOT ' ...
        'reconstruct or guess ROI membership on its own.'], roiDefFile);
end
if isempty(ver('stats'))
    error('Phase5A:MissingToolbox', ...
        'Statistics and Machine Learning Toolbox not detected. This script requires fitlme, anova(LinearMixedModel), coefTest, and qqplot.');
end

%% ------------------------------------------------------------------
%  1. LOAD X04 (channel-level BandPowerT) AND X07 (locked ROI definitions)
%     -- both shared, read-only, loaded ONCE before the Task loop.
% ------------------------------------------------------------------

S3 = load(bandPowerFile, 'BandPowerT', 'PARAMS3');
BandPowerT_ALL = S3.BandPowerT;
PARAMS3 = S3.PARAMS3;
clear S3;

requiredCols = {'Subject','TaskType','BlockNumber','Condition','Channel','Band', ...
                'AbsolutePower','LogAbsolutePower','RelativePower'};
missingCols = requiredCols(~ismember(requiredCols, BandPowerT_ALL.Properties.VariableNames));
if ~isempty(missingCols)
    error('Phase5A:MissingColumns', 'BandPowerT is missing expected column(s): %s', strjoin(missingCols, ', '));
end

S5D = load(roiDefFile, 'ChannelDefT', 'ROI_DEFINITION', 'ROI_NAMES', 'UNASSIGNED_FOR_ROI', 'PARAMS5D');
ChannelDefT = S5D.ChannelDefT;
ROI_DEFINITION = S5D.ROI_DEFINITION;
ROI_NAMES = S5D.ROI_NAMES;
UNASSIGNED_FOR_ROI = S5D.UNASSIGNED_FOR_ROI;
PARAMS5D = S5D.PARAMS5D;
clear S5D;

ALL_BANDS = fieldnames(PARAMS3.BandDefinitionsHz);
METRICS   = {'LogAbsolutePower', 'RelativePower'};
CONDITION_ORDER = {'E', 'N', 'R'};   % Baseline excluded from this analysis (see header item 2)
REF_LEVEL = CONDITION_ORDER{1};      % 'E' is the reference level for dummy coding

PRIMARY_TIER_BANDS     = {'Alpha'};
SECONDARY_TIER_BANDS   = {'Delta', 'Theta', 'Beta'};
EXPLORATORY_TIER_BANDS = {'Beta_low_exploratory', 'Beta_high_exploratory', ...
    'Gamma_low_exploratory', 'Gamma_high1_exploratory', 'Gamma_high2_exploratory'};
% Beta_low/Beta_high added (a sub-division of the unchanged locked
% Beta[13,30] band); each gets its own separate BH-FDR
% family exactly like the pre-existing Gamma exploratory bands.
expectedTierBands = [PRIMARY_TIER_BANDS, SECONDARY_TIER_BANDS, EXPLORATORY_TIER_BANDS];
if ~isempty(setxor(sort(expectedTierBands), sort(ALL_BANDS')))
    error('Phase5A:BandMismatch', ...
        ['The hard-coded band tiers do not exactly match PARAMS3.BandDefinitionsHz. ' ...
         'Tiers: %s | Phase3 bands: %s'], strjoin(expectedTierBands, ','), strjoin(ALL_BANDS', ','));
end

nBands = numel(ALL_BANDS);
nMetrics = numel(METRICS);
nROI = numel(ROI_NAMES);

TASK_TYPES = {'full', 'nback', 'sart'};
TASK_DISPLAY = struct('full', 'Full', 'nback', 'Nback', 'sart', 'SART');

%% ========================================================================
%  TASK LOOP -- everything below is repeated independently per Task
% ========================================================================
for taskI = 1:numel(TASK_TYPES)
    TASK_KEY = TASK_TYPES{taskI};
    TASK_LABEL = TASK_DISPLAY.(TASK_KEY);

    STATS_DIR = fullfile(ANALYSIS_ROOT, TASK_LABEL, '08_roi_stats', RUN_MODE);
    if ~exist(STATS_DIR, 'dir'); mkdir(STATS_DIR); end

    fprintf('\n########################################################\n');
    fprintf('### TASK = %s | RUN_MODE = %s\n', TASK_LABEL, RUN_MODE_LABEL);
    fprintf('########################################################\n\n');

    %% Section 1b: filter to this Task, then to Condition in {E,N,R}
    thisTaskT = BandPowerT_ALL(strcmp(BandPowerT_ALL.TaskType, TASK_KEY), :);
    nBaselineExcludedRows = sum(strcmp(thisTaskT.Condition, 'Baseline'));
    BandPowerT = thisTaskT(ismember(thisTaskT.Condition, CONDITION_ORDER), :);
    fprintf('Task=%s: %d total channel-level rows loaded, %d Baseline rows excluded, %d E/N/R rows retained for ROI aggregation.\n', ...
        TASK_LABEL, height(thisTaskT), nBaselineExcludedRows, height(BandPowerT));

    % ---- NEW: apply the shared noisy-subject exclusion list (see
    % EXCLUSION_CONFIG_cognitive.m), only when RUN_MODE=sensitivity_excl_noisy.
    % Never applied to the full_sample run. Applied to the channel-level
    % table BEFORE ROI aggregation (Section 3), so every ROI average for
    % an excluded subject is simply never computed. ----
    nExcludedRowsSensitivity = 0;
    if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
        subjNumForExcl = local_subject_numeric(BandPowerT.Subject);
        excludeMaskSensitivity = ismember(subjNumForExcl, EXCLUDE_NOISY_SUBJECTS);
        nExcludedRowsSensitivity = sum(excludeMaskSensitivity);
        BandPowerT = BandPowerT(~excludeMaskSensitivity, :);
        fprintf('Task=%s: RUN_MODE=%s -- excluded %d channel-level row(s) belonging to subject(s) {%s}; %d rows remain for ROI aggregation.\n', ...
            TASK_LABEL, RUN_MODE_LABEL, nExcludedRowsSensitivity, ...
            strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','), height(BandPowerT));
    end

    if isempty(BandPowerT)
        warning('Task=%s: no E/N/R rows found after filtering -- skipping this Task entirely.', TASK_LABEL);
        continue;
    end

    integrityIssues = {};

    %% ====================================================================
    %  SECTION 2: PRE-AGGREGATION VALIDATION (all performed before any
    %  model is fit). Structural membership/consistency checks are hard-
    %  gated (preValFails), exactly as in the original script. Exact
    %  recording-count literals are NOT gated -- see header item 4.
    % ====================================================================
    fprintf('--------------------------------------------------------\n');
    fprintf('SECTION 2: PRE-AGGREGATION VALIDATION (Task=%s)\n', TASK_LABEL);
    fprintf('--------------------------------------------------------\n');
    preValFails = 0;

    % --- 2a. All five ROI definitions loaded from X07 ---
    expectedROINames = {'Frontal','Central','Temporal','Parietal','Occipital'};
    if numel(ROI_NAMES) ~= 5 || ~isempty(setxor(ROI_NAMES, expectedROINames))
        fprintf('  *** FAIL: ROI_NAMES is not exactly the 5 expected ROIs. Got: %s\n', strjoin(ROI_NAMES, ', '));
        preValFails = preValFails + 1;
    else
        fprintf('  PASS: all five ROI definitions loaded (%s).\n', strjoin(ROI_NAMES, ', '));
    end
    for r = 1:numel(ROI_NAMES)
        if ~isfield(ROI_DEFINITION, ROI_NAMES{r}) || isempty(ROI_DEFINITION.(ROI_NAMES{r}))
            fprintf('  *** FAIL: ROI_DEFINITION.%s is missing or empty.\n', ROI_NAMES{r});
            preValFails = preValFails + 1;
        end
    end

    % --- 2b. All expected channels present in this Task's BandPowerT ---
    bandPowerChannels = unique(BandPowerT.Channel);
    refChannels = ChannelDefT.Label;
    missingFromBandPower = refChannels(~ismember(upper(refChannels), upper(bandPowerChannels)));
    extraInBandPower = bandPowerChannels(~ismember(upper(bandPowerChannels), upper(refChannels)));
    if ~isempty(missingFromBandPower)
        fprintf('  *** FAIL: %d reference channel(s) not found in BandPowerT: %s\n', ...
            numel(missingFromBandPower), strjoin(missingFromBandPower, ', '));
        preValFails = preValFails + 1;
    else
        fprintf('  PASS: all %d reference channels are present in BandPowerT.\n', numel(refChannels));
    end
    if ~isempty(extraInBandPower)
        fprintf('  *** FAIL: %d channel(s) in BandPowerT are not in the reference channel set: %s\n', ...
            numel(extraInBandPower), strjoin(extraInBandPower, ', '));
        preValFails = preValFails + 1;
    end
    allROIChannels = [ROI_DEFINITION.Frontal, ROI_DEFINITION.Central, ROI_DEFINITION.Temporal, ...
                      ROI_DEFINITION.Parietal, ROI_DEFINITION.Occipital];
    missingROIChan = allROIChannels(~ismember(upper(allROIChannels), upper(bandPowerChannels)));
    if ~isempty(missingROIChan)
        fprintf('  *** FAIL: %d ROI-member channel(s) not found in BandPowerT: %s\n', ...
            numel(missingROIChan), strjoin(missingROIChan, ', '));
        preValFails = preValFails + 1;
    else
        fprintf('  PASS: every channel referenced by ROI_DEFINITION is present in BandPowerT.\n');
    end

    % --- 2c. P9/P10 not included in any ROI average ---
    p9p10InRoi = intersect(upper(allROIChannels), upper({'P9','P10'}));
    if ~isempty(p9p10InRoi)
        fprintf('  *** FAIL: P9/P10 found inside ROI_DEFINITION (%s) -- must be excluded from ROI averaging.\n', strjoin(p9p10InRoi, ', '));
        preValFails = preValFails + 1;
    else
        fprintf('  PASS: P9 and P10 are absent from all five ROI_DEFINITION lists (excluded from ROI averaging).\n');
    end
    if ~all(ismember(upper({'P9','P10'}), upper(UNASSIGNED_FOR_ROI)))
        fprintf('  *** FAIL: P9/P10 not both present in UNASSIGNED_FOR_ROI as expected.\n');
        preValFails = preValFails + 1;
    else
        fprintf('  PASS: P9 and P10 are both explicitly recorded as unassigned-for-ROI (per X07).\n');
    end

    % --- 2d. Recording-level design verification (DATA-DRIVEN -- see
    %          header item 4; no exact-count literal is asserted here,
    %          only structural correctness). Mirrors X05's Section 2. ---
    comboT = unique(BandPowerT(:, {'Subject','BlockNumber','Condition'}), 'rows');
    nRecordings = height(comboT);
    fprintf('\n  %d unique Subject x BlockNumber x Condition recordings found (Task=%s).\n', nRecordings, TASK_LABEL);

    subjNumeric = local_subject_numeric(comboT.Subject);
    uniqSubj = unique(subjNumeric);
    nSubjects = numel(uniqSubj);

    subjPairs = cell(nSubjects, 1);
    pairTypeCount = containers.Map('KeyType', 'char', 'ValueType', 'double');
    duplicateSubjCondFlag = false;
    for i = 1:nSubjects
        rowsI = subjNumeric == uniqSubj(i);
        condsI = cellstr(comboT.Condition(rowsI));
        if numel(condsI) ~= 2
            % Reported, not hard-gated -- see header item 4: the exact
            % per-subject recording count is a property of this
            % not-yet-fully-audited dataset, not assumed in advance.
            integrityIssues{end+1} = sprintf('Subject %g has %d recordings (expected 2) for Task=%s.', uniqSubj(i), numel(condsI), TASK_LABEL); %#ok<AGROW>
        elseif strcmp(condsI{1}, condsI{2})
            % A genuine logical-impossibility check (not a count guess) --
            % kept hard-gated exactly as the original script.
            integrityIssues{end+1} = sprintf('Subject %g has the SAME condition (%s) in both recordings for Task=%s -- duplicate.', uniqSubj(i), condsI{1}, TASK_LABEL); %#ok<AGROW>
            duplicateSubjCondFlag = true;
            preValFails = preValFails + 1;
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
    fprintf('  Condition-pair coverage across %d subjects (Task=%s):\n', nSubjects, TASK_LABEL);
    pairKeys = keys(pairTypeCount);
    for k = 1:numel(pairKeys)
        fprintf('    %s : %d subjects\n', strrep(pairKeys{k}, '_', ' & '), pairTypeCount(pairKeys{k}));
    end
    if duplicateSubjCondFlag
        fprintf('  *** FAIL: at least one subject has a duplicated condition. See integrity issues.\n');
    else
        fprintf('  PASS: no subject has a duplicated condition.\n');
    end

    condCounts = struct();
    for c = 1:numel(CONDITION_ORDER)
        condCounts.(CONDITION_ORDER{c}) = sum(strcmp(cellstr(comboT.Condition), CONDITION_ORDER{c}));
    end
    fprintf('  Recording-level condition counts (Task=%s, reported -- not gated against a hard-coded literal, see header item 4): E=%d, N=%d, R=%d\n', ...
        TASK_LABEL, condCounts.E, condCounts.N, condCounts.R);

    % Structural design-assumption check (identical logic to X05): if
    % every subject in this Task's E/N/R data has all 3 conditions, the
    % incomplete-block-design justification for using an LMM instead of
    % RM-ANOVA no longer applies.
    allSubjectsComplete3 = false;
    condsPerSubjCounts = cellfun(@numel, subjPairs);
    if all(condsPerSubjCounts >= 3)
        allSubjectsComplete3 = true;
        integrityIssues{end+1} = sprintf(['Task=%s: every subject appears to have 3+ recordings/conditions -- the incomplete-design ' ...
            'justification for using an LMM instead of RM-ANOVA no longer applies. STOP and re-examine.'], TASK_LABEL); %#ok<AGROW>
    end

    fprintf('\nPre-aggregation STRUCTURAL validation failures: %d\n', preValFails);
    if preValFails > 0
        error('Phase5A:PreValidationFailed', ...
            'Task=%s: pre-aggregation validation failed (%d structural issue(s) above). Fix before aggregating/fitting.', TASK_LABEL, preValFails);
    end
    if allSubjectsComplete3
        error('Phase5A:DesignAssumptionViolated', ...
            'Task=%s: design verification found a complete 3-condition-per-subject structure. Stopping before ROI aggregation/fitting -- see integrityIssues.', TASK_LABEL);
    end

    %% ====================================================================
    %  SECTION 3: ROI AGGREGATION (channel-level BandPowerT -> ROI-level
    %     long table). Arithmetic-mean-of-channels rule; LogAbsolutePower
    %     = log(mean AbsolutePower), matching X04's own GlobalBandPowerT
    %     construction exactly (see header). P9/P10 are never included
    %     (their rows are simply never referenced by any ROI's channel
    %     list).
    % ====================================================================
    fprintf('\n--------------------------------------------------------\n');
    fprintf('SECTION 3: ROI AGGREGATION (Task=%s)\n', TASK_LABEL);
    fprintf('--------------------------------------------------------\n');

    BandPowerT.ChannelUC = upper(BandPowerT.Channel);

    roiRows = cell(0, 7);   % Subject, BlockNumber, Condition, ROI, Band, AbsolutePower, RelativePower
    roiCols = {'Subject','BlockNumber','Condition','ROI','Band','AbsolutePower','RelativePower'};

    recCombos = unique(BandPowerT(:, {'Subject','BlockNumber','Condition'}), 'rows');
    nonPositiveROIAbs = {};

    for i = 1:height(recCombos)
        subjStr = recCombos.Subject{i};
        blockNum = recCombos.BlockNumber(i);
        condStr = recCombos.Condition{i};
        recMask = strcmp(BandPowerT.Subject, subjStr) & BandPowerT.BlockNumber == blockNum & strcmp(BandPowerT.Condition, condStr);

        for r = 1:nROI
            roiName = ROI_NAMES{r};
            roiChansUC = upper(ROI_DEFINITION.(roiName));
            chanMask = ismember(BandPowerT.ChannelUC, roiChansUC);

            for b = 1:nBands
                bandName = ALL_BANDS{b};
                rowMask = recMask & chanMask & strcmp(BandPowerT.Band, bandName);
                nChanFound = sum(rowMask);
                if nChanFound ~= numel(roiChansUC)
                    integrityIssues{end+1} = sprintf('Sub%s Block%g Cond%s ROI=%s Band=%s: found %d/%d expected ROI channels.', ...
                        subjStr, blockNum, condStr, roiName, bandName, nChanFound, numel(roiChansUC)); %#ok<AGROW>
                end
                roiAbs = mean(BandPowerT.AbsolutePower(rowMask));
                roiRel = mean(BandPowerT.RelativePower(rowMask));
                if roiAbs <= 0
                    nonPositiveROIAbs{end+1} = sprintf('Sub%s Block%g Cond%s ROI=%s Band=%s', subjStr, blockNum, condStr, roiName, bandName); %#ok<AGROW>
                end
                roiRows(end+1, :) = {subjStr, blockNum, condStr, roiName, bandName, roiAbs, roiRel}; %#ok<AGROW>
            end
        end
    end

    RoiBandPowerT = local_build_table(roiRows, [1 3 4 5], [2 6 7], roiCols);
    RoiBandPowerT.LogAbsolutePower = log(RoiBandPowerT.AbsolutePower);
    if ~isempty(nonPositiveROIAbs)
        RoiBandPowerT.LogAbsolutePower(RoiBandPowerT.AbsolutePower <= 0) = NaN;
        integrityIssues{end+1} = sprintf('%d ROI x Band x recording combination(s) had non-positive ROI AbsolutePower -- LogAbsolutePower set to NaN there (never silently patched): %s', ...
            numel(nonPositiveROIAbs), strjoin(nonPositiveROIAbs, '; ')); %#ok<AGROW>
    end

    fprintf('ROI aggregation complete (Task=%s): %d rows (%d recordings x %d ROIs x %d bands).\n', ...
        TASK_LABEL, height(RoiBandPowerT), height(recCombos), nROI, nBands);

    %% ====================================================================
    %  SECTION 4: POST-AGGREGATION VALIDATION (still before any model is
    %     fit). Duplicate-key check is hard-gated (structural, not a
    %     guessed literal). Per-ROI-per-band condition counts are checked
    %     for CONSISTENCY against the dynamically-computed recording-level
    %     condCounts from Section 2 -- never against a hard-coded literal.
    %     See header item 4.
    % ====================================================================
    fprintf('\n--------------------------------------------------------\n');
    fprintf('SECTION 4: POST-AGGREGATION VALIDATION (Task=%s)\n', TASK_LABEL);
    fprintf('--------------------------------------------------------\n');
    postValFails = 0;

    % --- 4a. Every Subject x BlockNumber x ROI x Band observation unique
    %          (Metric is a column, not a row-defining key). ---
    keyT = RoiBandPowerT(:, {'Subject','BlockNumber','ROI','Band'});
    [uKeyT, ~, ~] = unique(keyT, 'rows');
    dupCount = height(keyT) - height(uKeyT);
    if dupCount > 0
        fprintf('  *** FAIL: %d duplicate Subject x BlockNumber x ROI x Band row(s) found.\n', dupCount);
        postValFails = postValFails + 1;
    else
        fprintf('  PASS: every Subject x BlockNumber x ROI x Band x Metric observation is uniquely defined (%d rows, %d unique keys).\n', ...
            height(keyT), height(uKeyT));
    end

    % --- 4b. No unexpected missing values ---
    nMissingLog = sum(isnan(RoiBandPowerT.LogAbsolutePower));
    nMissingRel = sum(isnan(RoiBandPowerT.RelativePower));
    fprintf('  Missing (NaN) LogAbsolutePower values: %d / %d\n', nMissingLog, height(RoiBandPowerT));
    fprintf('  Missing (NaN) RelativePower values   : %d / %d\n', nMissingRel, height(RoiBandPowerT));
    if nMissingLog > 0 || nMissingRel > 0
        fprintf('  *** NOTE: missing values present (see integrity issues above for exactly which ROI x Band x recording). Not treated as fatal -- fitlme will list-wise exclude affected rows, reported explicitly rather than silently imputed.\n');
    else
        fprintf('  PASS: no missing values in either metric.\n');
    end

    % --- 4c. Per-ROI, per-band recording-level condition counts must
    %          match the dynamically-computed Section-2 condCounts (i.e.
    %          ROI aggregation lost/duplicated nothing) -- NOT a
    %          hard-coded literal. See header item 4. ---
    condCountMismatch = false;
    for r = 1:nROI
        for b = 1:nBands
            subT = RoiBandPowerT(strcmp(RoiBandPowerT.ROI, ROI_NAMES{r}) & strcmp(RoiBandPowerT.Band, ALL_BANDS{b}), :);
            cE = sum(strcmp(subT.Condition, 'E')); cN = sum(strcmp(subT.Condition, 'N')); cR = sum(strcmp(subT.Condition, 'R'));
            if cE ~= condCounts.E || cN ~= condCounts.N || cR ~= condCounts.R
                condCountMismatch = true;
                integrityIssues{end+1} = sprintf('ROI=%s Band=%s (Task=%s): condition counts E=%d/N=%d/R=%d do not match the recording-level counts E=%d/N=%d/R=%d.', ...
                    ROI_NAMES{r}, ALL_BANDS{b}, TASK_LABEL, cE, cN, cR, condCounts.E, condCounts.N, condCounts.R); %#ok<AGROW>
            end
        end
    end
    if condCountMismatch
        fprintf('  *** FAIL: at least one ROI x Band combination does not preserve the recording-level condition counts (E=%d/N=%d/R=%d). See integrity issues.\n', ...
            condCounts.E, condCounts.N, condCounts.R);
        postValFails = postValFails + 1;
    else
        fprintf('  PASS: every ROI x Band combination (%d x %d = %d checks) preserves the recording-level condition counts (E=%d, N=%d, R=%d).\n', ...
            nROI, nBands, nROI*nBands, condCounts.E, condCounts.N, condCounts.R);
    end

    fprintf('\nPost-aggregation validation failures: %d\n', postValFails);
    if postValFails > 0
        error('Phase5A:PostValidationFailed', ...
            'Task=%s: post-aggregation validation failed (%d issue(s) above). Fix before fitting models.', TASK_LABEL, postValFails);
    end

    %% ====================================================================
    %  SECTION 5: MODEL FITTING -- 70 independent LMMs (5 ROI x 7 bands x
    %     2 metrics), per Task
    % ====================================================================
    fprintf('\n--------------------------------------------------------\n');
    fprintf('SECTION 5: MODEL FITTING (Task=%s, ROI x Band x Metric, %d total fits)\n', TASK_LABEL, nROI*nBands*nMetrics);
    fprintf('--------------------------------------------------------\n');

    omnibusRows = cell(0, 7);   % ROI, Metric, Band, FStat, DF1, DF2, pValue
    lmeStore = cell(nROI, nBands, nMetrics);
    convergenceWarnings = {};

    for r = 1:nROI
        roiName = ROI_NAMES{r};
        for m = 1:nMetrics
            metricName = METRICS{m};
            for b = 1:nBands
                bandName = ALL_BANDS{b};
                subT = RoiBandPowerT(strcmp(RoiBandPowerT.ROI, roiName) & strcmp(RoiBandPowerT.Band, bandName), :);
                subT.Subject = categorical(local_subject_numeric(subT.Subject));
                subT.Condition = categorical(cellstr(subT.Condition), CONDITION_ORDER);
                formula = sprintf('%s ~ Condition + (1|Subject)', metricName);
                try
                    lastwarn('');
                    lme = fitlme(subT, formula, 'FitMethod', 'REML');
                    [warnMsg, ~] = lastwarn();
                    if ~isempty(warnMsg)
                        convergenceWarnings{end+1} = sprintf('%s / %s / %s (Task=%s): %s', roiName, metricName, bandName, TASK_LABEL, warnMsg); %#ok<AGROW>
                    end
                catch ME
                    integrityIssues{end+1} = sprintf('fitlme FAILED for ROI=%s %s / %s (Task=%s): %s', roiName, metricName, bandName, TASK_LABEL, ME.message); %#ok<AGROW>
                    lmeStore{r, b, m} = [];
                    continue;
                end
                lmeStore{r, b, m} = lme;
                aov = anova(lme, 'DFMethod', 'satterthwaite');
                rowIdx = find(strcmp(aov.Term, 'Condition'), 1);
                if isempty(rowIdx)
                    integrityIssues{end+1} = sprintf('anova(lme) has no Condition row for ROI=%s %s / %s (Task=%s).', roiName, metricName, bandName, TASK_LABEL); %#ok<AGROW>
                    continue;
                end
                omnibusRows(end+1, :) = {roiName, metricName, bandName, aov.FStat(rowIdx), aov.DF1(rowIdx), aov.DF2(rowIdx), aov.pValue(rowIdx)}; %#ok<AGROW>
            end
        end
    end

    OmnibusT = local_build_table(omnibusRows, [1 2 3], [4 5 6 7], {'ROI','Metric','Band','FStat','DF1','DF2','pValue'});

    tierOf = containers.Map();
    for i = 1:numel(PRIMARY_TIER_BANDS);     tierOf(PRIMARY_TIER_BANDS{i}) = 'Primary'; end
    for i = 1:numel(SECONDARY_TIER_BANDS);   tierOf(SECONDARY_TIER_BANDS{i}) = 'Secondary'; end
    for i = 1:numel(EXPLORATORY_TIER_BANDS); tierOf(EXPLORATORY_TIER_BANDS{i}) = 'Exploratory'; end
    OmnibusT.Tier = cellfun(@(b) tierOf(b), OmnibusT.Band, 'UniformOutput', false);

    fprintf('Model fitting complete (Task=%s): %d / %d omnibus tests fit successfully.\n', TASK_LABEL, height(OmnibusT), nROI*nBands*nMetrics);

    %% ====================================================================
    %  SECTION 6: FDR -- separately within each Band x Metric family,
    %     across the 5 ROIs, per Task
    % ====================================================================
    fprintf('\n--------------------------------------------------------\n');
    fprintf('SECTION 6: BH-FDR (Task=%s, Band x Metric family, across 5 ROIs)\n', TASK_LABEL);
    fprintf('--------------------------------------------------------\n');

    FDR_ALPHA = 0.05;
    OmnibusT.pFDR = nan(height(OmnibusT), 1);
    OmnibusT.SignificantFDR = false(height(OmnibusT), 1);
    OmnibusT.FamilySize = nan(height(OmnibusT), 1);

    nFamilies = 0;
    for m = 1:nMetrics
        metricName = METRICS{m};
        for b = 1:nBands
            bandName = ALL_BANDS{b};
            idx = strcmp(OmnibusT.Metric, metricName) & strcmp(OmnibusT.Band, bandName);
            nInFamily = sum(idx);
            if nInFamily == 0; continue; end
            nFamilies = nFamilies + 1;
            OmnibusT.FamilySize(idx) = nInFamily;
            if nInFamily ~= nROI
                integrityIssues{end+1} = sprintf('FDR family %s / %s (Task=%s) has %d members, expected %d (one per ROI) -- some ROI fit may have failed.', ...
                    metricName, bandName, TASK_LABEL, nInFamily, nROI); %#ok<AGROW>
            end
            pvals = OmnibusT.pValue(idx);
            [adjP, reject] = local_fdr_bh(pvals, FDR_ALPHA);
            OmnibusT.pFDR(idx) = adjP;
            OmnibusT.SignificantFDR(idx) = reject;
        end
    end
    fprintf('BH-FDR applied within %d Band x Metric families (expected %d = %d bands x %d metrics), each spanning up to %d ROIs.\n', ...
        nFamilies, nBands*nMetrics, nBands, nMetrics, nROI);

    fprintf('\n--- ROI-level omnibus results (Task=%s) ---\n', TASK_LABEL);
    disp(OmnibusT);

    %% ====================================================================
    %  SECTION 7: PAIRWISE CONTRASTS (E-N, E-R, N-R -- preserved EXACTLY
    %     from the original 08 script's own convention, see header item
    %     5; NOT the N-E/R-E/R-N convention used by 05/06_..._cognitive.m),
    %     gated on the Band x Metric FDR result, Holm-corrected within
    %     each surviving ROI x Band x Metric
    % ====================================================================
    fprintf('\n--------------------------------------------------------\n');
    fprintf('SECTION 7: PAIRWISE CONTRASTS (Task=%s, gated on FDR, Holm-corrected)\n', TASK_LABEL);
    fprintf('--------------------------------------------------------\n');

    pairwiseRows = cell(0, 12);   % ROI,Metric,Band,Contrast,Estimate,SE,DF,tStat,pValue,pHolm,CI_Low,CI_High
    CONTRASTS = {'E','N','E - N'; 'E','R','E - R'; 'N','R','N - R'};   % {lvl2, lvl1, label}; label = lvl2 - lvl1 -- IDENTICAL to the original 08 script

    sigMask = OmnibusT.SignificantFDR;
    for row = find(sigMask)'
        roiName = OmnibusT.ROI{row};
        metricName = OmnibusT.Metric{row};
        bandName = OmnibusT.Band{row};
        r = find(strcmp(ROI_NAMES, roiName));
        b = find(strcmp(ALL_BANDS, bandName));
        m = find(strcmp(METRICS, metricName));
        lme = lmeStore{r, b, m};
        if isempty(lme); continue; end
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
            pairwiseRows(end+1, :) = [{roiName, metricName, bandName}, contrastData(c, 1:6), {adjP(c)}, contrastData(c, 7:8)]; %#ok<AGROW>
        end
    end
    PairwiseT = local_build_table(pairwiseRows, [1 2 3 4], [5 6 7 8 9 10 11 12], ...
        {'ROI','Metric','Band','Contrast','Estimate','SE','DF','tStat','pValue','pHolm','CI_Low','CI_High'});

    fprintf('\n--- Pairwise contrasts (Task=%s, ROI x Band x Metric significant under Band x Metric FDR) ---\n', TASK_LABEL);
    disp(PairwiseT);

    %% --------------------------------------------------------------
    %  SECTION 8: MODEL DIAGNOSTICS (Primary + Secondary tier bands only,
    %     all 5 ROIs, both metrics -- Alpha/Delta/Theta/Beta x 5 ROI x 2
    %     metrics = 40 plots per Task. Exploratory/gamma bands are
    %     skipped here, matching the original script's own diagBands
    %     convention.)
    % --------------------------------------------------------------
    fprintf('\n--------------------------------------------------------\n');
    fprintf('SECTION 8: MODEL DIAGNOSTICS (Task=%s)\n', TASK_LABEL);
    fprintf('--------------------------------------------------------\n');

    diagBands = [PRIMARY_TIER_BANDS, SECONDARY_TIER_BANDS];
    nDiagSaved = 0;
    for r = 1:nROI
        roiName = ROI_NAMES{r};
        for m = 1:nMetrics
            metricName = METRICS{m};
            for bi = 1:numel(diagBands)
                bandName = diagBands{bi};
                b = find(strcmp(ALL_BANDS, bandName));
                lme = lmeStore{r, b, m};
                if isempty(lme); continue; end
                try
                    res = residuals(lme, 'ResidualType', 'Pearson');
                    fig = figure('Visible', 'off', 'Position', [100 100 900 400]);
                    subplot(1,2,1);
                    histogram(res, 15);
                    title(sprintf('Residuals: %s / %s (%s), Task=%s', roiName, bandName, metricName, TASK_LABEL), 'Interpreter', 'none');
                    xlabel('Pearson residual'); ylabel('Count');
                    subplot(1,2,2);
                    qqplot(res);
                    title('Normal Q-Q');
                    outPng = fullfile(STATS_DIR, sprintf('diagnostic_residuals_%s_%s_%s.png', roiName, bandName, metricName));
                    saveas(fig, outPng);
                    close(fig);
                    nDiagSaved = nDiagSaved + 1;
                catch ME
                    integrityIssues{end+1} = sprintf('Residual diagnostic failed for ROI=%s %s / %s (Task=%s): %s', roiName, metricName, bandName, TASK_LABEL, ME.message); %#ok<AGROW>
                end
            end
        end
    end
    fprintf('Saved %d residual diagnostic plots (Task=%s; Primary+Secondary bands x 5 ROIs x 2 metrics; expected up to %d).\n', ...
        nDiagSaved, TASK_LABEL, nROI*numel(diagBands)*nMetrics);

    %% --------------------------------------------------------------
    %  SECTION 9: PARAMETERS AND OUTPUT
    % --------------------------------------------------------------

    PARAMS5A = struct();
    PARAMS5A.ScriptName = 'X08_roi_level_statistics_cognitive.m';
    PARAMS5A.Task = TASK_LABEL;
    PARAMS5A.GeneratedOn = datestr(now); %#ok<TNOW1,DATST>
    PARAMS5A.AnalysisLevel = sprintf('ROI-averaged band power (X07 locked 5-ROI definition, applied to X04 channel-level BandPowerT), Task=%s only, Condition in {E,N,R} only (Baseline excluded)', TASK_LABEL);
    PARAMS5A.Model = 'LinearMixedEffects_RandomIntercept';
    PARAMS5A.FixedEffectsPrimary = 'Condition (3-level: E, N, R; reference = E)';
    PARAMS5A.RandomEffects = '(1|Subject)';
    PARAMS5A.FitMethod = 'REML';
    PARAMS5A.DFMethod = 'Satterthwaite';
    PARAMS5A.Metrics = METRICS;
    PARAMS5A.ROINames = ROI_NAMES;
    PARAMS5A.UnassignedForROI = {UNASSIGNED_FOR_ROI};
    PARAMS5A.ROIAggregationRule = 'Arithmetic mean of channel-level AbsolutePower/RelativePower across ROI member channels; LogAbsolutePower = log(mean AbsolutePower), matching X04 GlobalBandPowerT construction.';
    PARAMS5A.BandTiers = struct('Primary', {PRIMARY_TIER_BANDS}, 'Secondary', {SECONDARY_TIER_BANDS}, 'Exploratory', {EXPLORATORY_TIER_BANDS});
    PARAMS5A.FDRAlpha = FDR_ALPHA;
    PARAMS5A.FDRFamilyScheme = 'Separate BH-FDR family per (Band x Metric), across the 5 ROIs (family size = 5). NOT tier-based. Bands/metrics/ROIs never pooled.';
    PARAMS5A.PairwiseContrastMethod = 'Holm-Bonferroni across the 3 pairwise contrasts (E-N, E-R, N-R -- the ORIGINAL 08 script''s own convention, preserved exactly; NOTE this differs from the N-E/R-E/R-N convention used by 05/06_..._cognitive.m), computed only for ROI x Band x Metric combinations whose Band x Metric FDR family test survives BH-FDR.';
    PARAMS5A.EffectSizePolicy = 'Model-based estimated difference + 95%% CI only. No standardized effect-size denominator computed in this script.';
    PARAMS5A.BaselineExclusionNote = sprintf('%d Baseline channel-level row(s) excluded before ROI aggregation; Baseline is aggregated upstream (X04) but never enters this or any other inferential model in this pipeline.', nBaselineExcludedRows);
    PARAMS5A.DesignVerificationNote = 'Exact recording-count literals are NOT gated (see script header item 4) -- structural checks only; exact counts are reported in DesignVerification below and in the report.';
    PARAMS5A.SourcePhase3File = bandPowerFile;
    PARAMS5A.SourcePhase5DFile = roiDefFile;
    PARAMS5A.SourcePhase5DDecision = PARAMS5D.ROIDecision;
    PARAMS5A.DesignVerification = struct('nRecordings', nRecordings, 'nSubjects', nSubjects, ...
        'PairTypeCounts', {pairTypeCount.keys}, 'ConditionCounts', condCounts, 'IntegrityIssues', {integrityIssues});
    PARAMS5A.RunDateTime = datestr(now); %#ok<TNOW1,DATST>
    PARAMS5A.RunMode = RUN_MODE;
    if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
        PARAMS5A.ExcludedSubjects = EXCLUDE_NOISY_SUBJECTS;
        PARAMS5A.ExcludedSubjectsRationale = EXCLUSION_RATIONALE_cognitive;
        PARAMS5A.RowsExcludedForSensitivity = nExcludedRowsSensitivity;
    else
        PARAMS5A.ExcludedSubjects = [];
    end

    save(fullfile(STATS_DIR, 'phase5a_parameters.mat'), 'PARAMS5A');
    save(fullfile(STATS_DIR, 'phase5a_roi_bandpower.mat'), 'RoiBandPowerT', 'PARAMS5A');
    save(fullfile(STATS_DIR, 'phase5a_roi_omnibus_results.mat'), 'OmnibusT');
    save(fullfile(STATS_DIR, 'phase5a_roi_pairwise_contrasts.mat'), 'PairwiseT');

    writetable(OmnibusT, fullfile(STATS_DIR, 'phase5a_roi_omnibus_results.xlsx'));
    writetable(OmnibusT, fullfile(STATS_DIR, 'phase5a_roi_omnibus_results.csv'));
    if ~isempty(PairwiseT)
        writetable(PairwiseT, fullfile(STATS_DIR, 'phase5a_roi_pairwise_contrasts.xlsx'));
        writetable(PairwiseT, fullfile(STATS_DIR, 'phase5a_roi_pairwise_contrasts.csv'));
    end

    fprintf('\nSaved output to: %s\n', STATS_DIR);

    %% --------------------------------------------------------------
    %  SECTION 10: REPORT
    % --------------------------------------------------------------

    reportLines = {};
    reportLines{end+1} = repmat('=', 1, 60);
    reportLines{end+1} = sprintf('PHASE 5A -- ROI-LEVEL E vs N vs R STATISTICS REPORT (Task=%s)', TASK_LABEL);
    reportLines{end+1} = repmat('=', 1, 60);
    reportLines{end+1} = '';
    reportLines{end+1} = sprintf('Baseline channel-level rows excluded before ROI aggregation: %d', nBaselineExcludedRows);
    reportLines{end+1} = sprintf('Recordings verified: %d | Subjects: %d', nRecordings, nSubjects);
    reportLines{end+1} = sprintf('Condition-pair coverage: %s', strjoin(cellfun(@(k) sprintf('%s=%d', strrep(k,'_',' & '), pairTypeCount(k)), pairKeys, 'UniformOutput', false), ', '));
    reportLines{end+1} = sprintf('Recording-level condition counts (reported, not gated against a hard-coded literal -- see header item 4): E=%d, N=%d, R=%d', condCounts.E, condCounts.N, condCounts.R);
    reportLines{end+1} = '';
    reportLines{end+1} = 'Model: Linear mixed-effects, Condition (fixed, 3-level) + (1|Subject), REML, Satterthwaite DF.';
    reportLines{end+1} = sprintf('Analysis level: ROI-averaged band power (X07 locked 5-ROI definition; X04 channel-level BandPowerT), Task=%s only.', TASK_LABEL);
    reportLines{end+1} = sprintf('ROIs (%d): %s', nROI, strjoin(ROI_NAMES, ', '));
    reportLines{end+1} = sprintf('Unassigned for ROI averaging: %s', strjoin(UNASSIGNED_FOR_ROI, ', '));
    reportLines{end+1} = 'Tiers (labeling only, not FDR family boundaries): Primary=Alpha; Secondary=Delta/Theta/Beta; Exploratory=3 gamma bands.';
    reportLines{end+1} = 'Metrics: LogAbsolutePower and RelativePower, co-primary, parallel, never pooled for FDR.';
    reportLines{end+1} = 'FDR: BH-FDR within each Band x Metric family, across the 5 ROIs (14 families x 5 ROIs = 70 tests).';
    reportLines{end+1} = 'Pairwise: E-N, E-R, N-R (this script''s own original convention -- see header item 5); computed only for FDR-surviving ROI x Band x Metric; Holm-corrected within each (family size 3).';
    reportLines{end+1} = 'Effect sizes: model-based estimated difference + 95%% CI only. No standardized effect size in this script.';
    reportLines{end+1} = '';
    reportLines{end+1} = '--- ROI-level omnibus results (all 70 tests) ---';
    for row = 1:height(OmnibusT)
        reportLines{end+1} = sprintf('  [%s] %-10s %-25s Tier=%-11s F(%g,%.1f)=%.3f p=%.4f pFDR=%.4f family_n=%d %s', ...
            OmnibusT.Metric{row}, OmnibusT.ROI{row}, OmnibusT.Band{row}, OmnibusT.Tier{row}, OmnibusT.DF1(row), OmnibusT.DF2(row), ...
            OmnibusT.FStat(row), OmnibusT.pValue(row), OmnibusT.pFDR(row), OmnibusT.FamilySize(row), ...
            local_yn(OmnibusT.SignificantFDR(row))); %#ok<AGROW>
    end
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Pairwise contrasts (FDR-surviving ROI x Band x Metric only; Holm-corrected within each) ---';
    if isempty(PairwiseT)
        reportLines{end+1} = '  (none -- no ROI x Band x Metric survived its Band x Metric FDR family)';
    else
        for row = 1:height(PairwiseT)
            reportLines{end+1} = sprintf('  [%s] %-10s %-25s %-8s est=%.4f (95%% CI %.4f to %.4f) t=%.3f df=%.1f p=%.4f pHolm=%.4f', ...
                PairwiseT.Metric{row}, PairwiseT.ROI{row}, PairwiseT.Band{row}, PairwiseT.Contrast{row}, PairwiseT.Estimate(row), ...
                PairwiseT.CI_Low(row), PairwiseT.CI_High(row), PairwiseT.tStat(row), PairwiseT.DF(row), ...
                PairwiseT.pValue(row), PairwiseT.pHolm(row)); %#ok<AGROW>
        end
    end
    reportLines{end+1} = '';
    reportLines{end+1} = sprintf('Residual diagnostic plots saved: %d', nDiagSaved);
    reportLines{end+1} = '';
    if ~isempty(integrityIssues)
        reportLines{end+1} = '--- INTEGRITY ISSUES (reported, not auto-corrected) ---';
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
    reportLines{end+1} = 'This script performs NO interpretation beyond the validation/statistics';
    reportLines{end+1} = 'specified above. Electrode-level analysis and topography (if written for';
    reportLines{end+1} = 'this dataset) are separate scripts, not this one.';
    reportLines{end+1} = '';
    reportLines{end+1} = 'Output files written to:';
    reportLines{end+1} = ['  ' STATS_DIR];
    reportLines{end+1} = repmat('=', 1, 60);

    reportText = strjoin(reportLines, '\n');
    fprintf('%s\n', reportText);
    fid = fopen(fullfile(REPORT_DIR, sprintf('phase5a_roi_level_statistics_report_%s_%s.txt', TASK_LABEL, RUN_MODE_LABEL)), 'w');
    fprintf(fid, '%s\n', reportText);
    fclose(fid);

end % TASK_TYPES loop

fprintf('\n=== X08 complete for all three Tasks. STOP for review before writing any Phase 5B/electrode-level script. ===\n');

%% ======================= LOCAL FUNCTIONS (all after top-level code) =======================

function n = local_subject_numeric(subjVal)
% Robustly extract a numeric subject id regardless of whether Subject is
% stored as numeric, categorical, char, or cellstr, and regardless of
% zero-padding (e.g. '7' vs '07'). Identical implementation to the
% original 08_roi_level_statistics.m / X05 / X06.
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
% Benjamini-Hochberg step-up FDR procedure. Identical implementation to
% the original 08_roi_level_statistics.m / X05 / X06.
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
% Holm step-down procedure (controls family-wise error rate). Identical
% implementation to the original 08_roi_level_statistics.m / X05 / X06.
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
% Build a contrast row (length = numel(coeffNames)) representing
% mean(level2) - mean(level1) under reference-coding, by name-matching
% against the fitted model's actual CoefficientNames (never assumes
% column position/order). Identical implementation to the original
% 08_roi_level_statistics.m / X05 / X06.
    cRow = zeros(1, numel(coeffNames));
    if ~strcmp(level1, refLevel)
        name1 = sprintf('%s_%s', varName, level1);
        idx1 = find(strcmp(coeffNames, name1));
        if isempty(idx1)
            error('Phase5A:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name1);
        end
        cRow(idx1) = cRow(idx1) - 1;
    end
    if ~strcmp(level2, refLevel)
        name2 = sprintf('%s_%s', varName, level2);
        idx2 = find(strcmp(coeffNames, name2));
        if isempty(idx2)
            error('Phase5A:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name2);
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
% Build a table column-by-column via the table() constructor instead of
% cell2table(), so a text column is never silently reinterpreted.
% Identical implementation to the original 08_roi_level_statistics.m /
% X05 / X06.
    nCols = numel(varNames);
    if size(rows, 2) ~= nCols
        error('Phase5A:BuildTableColumnMismatch', ...
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
            error('Phase5A:BuildTableColumnMismatch', 'Column %d not classified as text or numeric.', j);
        end
    end
    T = table(cols{:}, 'VariableNames', varNames);
end
