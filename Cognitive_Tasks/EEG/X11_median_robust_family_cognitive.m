%% ========================================================================
% X11_median_robust_family_cognitive.m
%
% Purpose: median-robust derived-metric family (exploratory extension
% family 1 of 5: X11-X15). Adds derived metrics (median-based band power,
% epoch-to-epoch coefficient of variation, and broadband total power,
% each at Electrode / ROI / Global scalp level) alongside the two locked
% co-primary metrics (LogAbsolutePower, RelativePower) computed in
% X05/X06/X08/X09/X10. Runs as three separate, parallel analyses, one per
% Task; Task is never modeled as a factor and Task x Condition is never
% fit. Never modifies, re-derives, or reads the statistical outputs of
% X05/X06/X08/X09/X10; only reads X03's raw per-window PSD output, X04's
% subject-level aggregation (BandPowerT), and X07's shared ROI/electrode
% definitions, all read-only. Never pools its new metrics with
% LogAbsolutePower/RelativePower for FDR, and never pools its 4
% band-dependent metrics with each other or with its 2 broadband
% TotalPower metrics.
%
% Method: EEG band power is right-skewed; a few outlier epochs (below the
% hard QC rejection threshold applied in X03) can inflate a mean-based
% estimate and its variance, so median-based aggregation is used as a
% more robust alternative. Because mean and median do not commute, the
% order of operations differs by level: at Electrode level,
% median-across-valid-epochs per channel (reused directly from X04's
% WindowMedian_Absolute/WindowMedian_Relative, not recomputed); at ROI
% level, median-across-epochs per channel first, then arithmetic mean
% across the ROI's member channels (ROIs have few member channels, so
% spatial averaging alone would not smooth outliers much); at Global
% level, arithmetic mean-across-channels first (per valid epoch, fresh
% from X03's psd_data), then median-across-epochs of that per-epoch
% global series (the opposite order, since 62-channel spatial averaging
% already smooths spatial noise). The same per-level order-of-operations
% logic is applied to a second, broadband pair of metrics
% (TotalPowerMean/TotalPowerMedian) built fresh from X03's
% TotalPower_1to45Hz field.
%
% What this script computes (per Task, independently): band-dependent
% (per Subject x Condition[E/N/R] x Band, at Electrode/ROI/Global level)
% AbsMedian, RelMedian, CV_Absolute (= WindowSD_Absolute / AbsolutePower),
% CV_Relative; broadband (per Subject x Condition, no Band dimension)
% TotalPowerMean, TotalPowerMedian.
%
% Statistical model (identical fitting machinery to every other script in
% this pipeline, applied to the new metrics only): Metric ~ Condition +
% (1|Subject), REML, Satterthwaite DF, Condition in {E,N,R} only
% (Baseline excluded).
%
% Multiple-comparison families (own, separate families; never pooled with
% each other): band-dependent metrics (AbsMedian, RelMedian, CV_Absolute,
% CV_Relative) get BH-FDR per Band x Metric across electrodes (family
% size up to 62) at Electrode level, per Band x Metric across ROIs
% (family size up to 5) at ROI level, and per Tier x Metric across the
% bands in that tier at Global level (Primary = Alpha alone, family size
% 1, reported as-is; Secondary = Delta/Theta/Beta, family size 3;
% Exploratory = 3 gamma bands, family size 3); broadband TotalPower
% metrics get BH-FDR per Metric across electrodes (Electrode level), per
% Metric across ROIs (ROI level), or a single test (Global level, family
% size 1). Pairwise contrasts (N-E, R-E, R-N, matching X05/X06's own
% convention) are Holm-corrected (family size 3), computed only for units
% whose own FDR family test survives BH-FDR.
%
% Judgment calls (also written into PARAMS11.JudgmentCalls and the
% per-Task report): (1) ROI-level CV_Absolute/CV_Relative = arithmetic
% mean of per-channel CV across the ROI's member channels -- a reasonable
% summary, not literally the CV of the ROI-aggregated series (that would
% require a different, epoch-level ROI-mean computation, used only at
% Global level); (2) Global-level CV_Absolute/CV_Relative = SD-across-
% epochs of the per-epoch mean-across-channels series, divided by its
% mean; (3) residual diagnostics are scoped to Global level only (Primary
% + Secondary tier bands, both band-dependent metrics), mirroring X08's
% own restriction to Primary+Secondary bands; (4) a recording with 0
% valid windows is skipped for Global/TotalPower rows and logged, rather
% than fabricating a row, matching X04's own "skip aggregation, still
% record QC" pattern.
%
% Requirements: base MATLAB + Statistics and Machine Learning Toolbox
% (fitlme, anova on LinearMixedModel, coefTest).
%
% Inputs (read-only, never modified):
%   <ANALYSIS_ROOT>/03_psd/cognitive_psd.mat            (psd_data)
%   <ANALYSIS_ROOT>/04_subject_level/phase3_subject_level_bandpower.mat (BandPowerT)
%   <ANALYSIS_ROOT>/07_roi_electrode_definitions/phase5_channel_roi_definitions.mat
%
% Outputs (one complete set per Task; nothing under X03/X04/X05/X06/X07/
% X08/X09/X10's own output folders is ever written to):
%   <ANALYSIS_ROOT>/<Full|Nback|SART>/11_median_robust_family/
%       median_robust_data.mat   (MedianRobustElectrodeT, MedianRobustROIT,
%                                  MedianRobustGlobalT, TotalPowerElectrodeT,
%                                  TotalPowerROIT, TotalPowerGlobalT --
%                                  schema required by X12/X13)
%       phase11_parameters.mat
%       phase11_omnibus_band_results.mat / .csv
%       phase11_omnibus_totalpower_results.mat / .csv
%       phase11_pairwise_band_contrasts.mat / .csv
%       phase11_pairwise_totalpower_contrasts.mat / .csv
%       diagnostic_residuals_Global_<band>_<metric>.png
%   <ANALYSIS_ROOT>/09_reports/
%       phase_median_robust_family_report_<Full|Nback|SART>.txt
%
% RUN_MODE: full_sample / sensitivity_excl_noisy (see
% EXCLUSION_CONFIG_cognitive.m); run this script once per mode. Because
% X12 and X13 read this script's per-Task output file, both look under
% this script's RUN_MODE subfolder for it.
% ========================================================================

clear; clc;

%% Section 0: Paths and existence checks
config;
DATA_ROOT_FOR_SIBLING = DATA_DIR;
[DATA_ROOT_PARENT, ~] = fileparts(DATA_ROOT_FOR_SIBLING);
ANALYSIS_ROOT = OUT_DIR;

SUBJ_DIR   = fullfile(ANALYSIS_ROOT, '04_subject_level');
PSD_DIR    = fullfile(ANALYSIS_ROOT, '03_psd');
DEF_DIR    = fullfile(ANALYSIS_ROOT, '07_roi_electrode_definitions');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');

bandPowerFile = fullfile(SUBJ_DIR, 'phase3_subject_level_bandpower.mat');
psdMatFile    = fullfile(PSD_DIR, 'cognitive_psd.mat');
roiDefFile    = fullfile(DEF_DIR, 'phase5_channel_roi_definitions.mat');

if ~exist(bandPowerFile, 'file')
    error('Phase11:MissingInput', ['X04 output not found: %s\n' ...
        'Run X04_aggregate_cognitive_subject_level.m first.'], bandPowerFile);
end
if ~exist(psdMatFile, 'file')
    error('Phase11:MissingInput', ['X03 output not found: %s\n' ...
        'Run X03_compute_cognitive_psd.m first.'], psdMatFile);
end
if ~exist(roiDefFile, 'file')
    error('Phase11:MissingInput', ['X07 output not found: %s\n' ...
        'Run X07_roi_electrode_definitions_cognitive.m first. This script never ' ...
        'reconstructs or guesses ROI membership on its own.'], roiDefFile);
end
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
    error('Phase11:BadRunMode', 'RUN_MODE must be one of: %s', strjoin(RUN_MODES_cognitive, ', '));
end
RUN_MODE_LABEL = RUN_MODES_DISPLAY_cognitive.(RUN_MODE);
fprintf('RUN_MODE = %s\n', RUN_MODE_LABEL);
if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
    fprintf('Excluding subject(s) {%s} from every model in this run -- see EXCLUSION_CONFIG_cognitive.m for rationale.\n', ...
        strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','));
end
% ---- END NEW ----

if isempty(ver('stats'))
    error('Phase11:MissingToolbox', ...
        'Statistics and Machine Learning Toolbox not detected. This script requires fitlme, anova(LinearMixedModel), and coefTest.');
end

%% Section 1: Load shared, read-only inputs (X03, X04, X07 ONLY -- never
%  X05/X06/X08/X09/X10's own saved statistical outputs)
S4 = load(bandPowerFile, 'BandPowerT', 'PARAMS3');
BandPowerT_ALL = S4.BandPowerT;
PARAMS3 = S4.PARAMS3;
clear S4;

requiredBPCols = {'Subject','TaskType','BlockNumber','Condition','Channel','Band', ...
    'AbsolutePower','RelativePower','WindowMedian_Absolute','WindowSD_Absolute', ...
    'WindowMedian_Relative','WindowSD_Relative'};
missingBP = requiredBPCols(~ismember(requiredBPCols, BandPowerT_ALL.Properties.VariableNames));
if ~isempty(missingBP)
    error('Phase11:MissingColumns', 'BandPowerT is missing expected column(s): %s', strjoin(missingBP, ', '));
end

S7 = load(roiDefFile, 'ChannelDefT', 'ROI_DEFINITION', 'ROI_NAMES', 'UNASSIGNED_FOR_ROI', 'PARAMS5D');
ChannelDefT = S7.ChannelDefT; %#ok<NASGU>
ROI_DEFINITION = S7.ROI_DEFINITION;
ROI_NAMES = S7.ROI_NAMES;
UNASSIGNED_FOR_ROI = S7.UNASSIGNED_FOR_ROI;
PARAMS5D = S7.PARAMS5D;
clear S7;

S3 = load(psdMatFile, 'psd_data', 'ChannelLabelsReference', 'PARAMS');
psd_data = S3.psd_data;
ChannelLabelsReference = S3.ChannelLabelsReference(:);
PARAMS2 = S3.PARAMS;
clear S3;

ALL_BANDS = fieldnames(PARAMS3.BandDefinitionsHz);
if ~isempty(setxor(sort(ALL_BANDS), sort(fieldnames(PARAMS2.BandDefinitionsHz))))
    error('Phase11:BandMismatch', 'X04''s PARAMS3 bands and X03''s PARAMS bands do not match.');
end
if ~isfield(psd_data, 'AbsoluteBandPower') || ~all(ismember(ALL_BANDS, fieldnames(psd_data(1).AbsoluteBandPower)))
    error('Phase11:BandMismatch', 'psd_data(1).AbsoluteBandPower does not contain all expected bands.');
end

METRICS_BAND  = {'AbsMedian', 'RelMedian', 'CV_Absolute', 'CV_Relative'};
METRICS_TOTAL = {'TotalPowerMean', 'TotalPowerMedian'};
CONDITION_ORDER = {'E', 'N', 'R'};   % Baseline excluded from this analysis (see header)
REF_LEVEL = CONDITION_ORDER{1};
CONTRASTS = {'N','E','N - E'; 'R','E','R - E'; 'R','N','R - N'};   % N-E, R-E, R-N -- matches X05/X06 convention (see header)

PRIMARY_TIER_BANDS     = {'Alpha'};
SECONDARY_TIER_BANDS   = {'Delta', 'Theta', 'Beta'};
EXPLORATORY_TIER_BANDS = {'Beta_low_exploratory', 'Beta_high_exploratory', ...
    'Gamma_low_exploratory', 'Gamma_high1_exploratory', 'Gamma_high2_exploratory'};
% Beta_low/Beta_high added (a sub-division of the
% unchanged locked Beta[13,30] band); each gets its own separate BH-FDR
% family exactly like the pre-existing Gamma exploratory bands.
expectedTierBands = [PRIMARY_TIER_BANDS, SECONDARY_TIER_BANDS, EXPLORATORY_TIER_BANDS];
if ~isempty(setxor(sort(expectedTierBands), sort(ALL_BANDS')))
    error('Phase11:BandMismatch', 'Hard-coded band tiers do not match PARAMS3.BandDefinitionsHz.');
end
tierDefs = {'Primary', PRIMARY_TIER_BANDS; 'Secondary', SECONDARY_TIER_BANDS; 'Exploratory', EXPLORATORY_TIER_BANDS};
tierOf = containers.Map();
for i = 1:numel(PRIMARY_TIER_BANDS);     tierOf(PRIMARY_TIER_BANDS{i}) = 'Primary'; end     %#ok<NOCOL>
for i = 1:numel(SECONDARY_TIER_BANDS);   tierOf(SECONDARY_TIER_BANDS{i}) = 'Secondary'; end %#ok<NOCOL>
for i = 1:numel(EXPLORATORY_TIER_BANDS); tierOf(EXPLORATORY_TIER_BANDS{i}) = 'Exploratory'; end %#ok<NOCOL>

FDR_ALPHA = 0.05;
nBands = numel(ALL_BANDS);
nROI = numel(ROI_NAMES);
nChanRef = numel(ChannelLabelsReference);

TASK_TYPES = {'full', 'nback', 'sart'};
TASK_DISPLAY = struct('full', 'Full', 'nback', 'Nback', 'sart', 'SART');

fprintf('========================================================\n');
fprintf('PHASE 11 (X11) -- MEDIAN-ROBUST DERIVED-METRIC FAMILY, COGNITIVE-TASK DATA\n');
fprintf('Exploratory ADDITION alongside the locked co-primary metrics (X05/X06/X08/X09/X10).\n');
fprintf('Running as THREE INDEPENDENT, PARALLEL analyses (Full / Nback / SART).\n');
fprintf('========================================================\n\n');

%% ========================================================================
%  TASK LOOP -- everything below is repeated independently per Task
% ========================================================================
for taskI = 1:numel(TASK_TYPES)
    TASK_KEY = TASK_TYPES{taskI};
    TASK_LABEL = TASK_DISPLAY.(TASK_KEY);

    STATS_DIR = fullfile(ANALYSIS_ROOT, TASK_LABEL, '11_median_robust_family', RUN_MODE);
    if ~exist(STATS_DIR, 'dir'); mkdir(STATS_DIR); end

    fprintf('\n########################################################\n');
    fprintf('### TASK = %s | RUN_MODE = %s\n', TASK_LABEL, RUN_MODE_LABEL);
    fprintf('########################################################\n\n');

    integrityIssues = {};

    %% Section 2: filter to this Task, then to Condition in {E,N,R} (X04 side)
    thisTaskBP = BandPowerT_ALL(strcmp(BandPowerT_ALL.TaskType, TASK_KEY), :);
    nBaselineExcludedBP = sum(strcmp(thisTaskBP.Condition, 'Baseline'));
    BandPowerT = thisTaskBP(ismember(thisTaskBP.Condition, CONDITION_ORDER), :);
    fprintf('Task=%s: %d channel-level rows loaded from X04, %d Baseline rows excluded, %d E/N/R rows retained.\n', ...
        TASK_LABEL, height(thisTaskBP), nBaselineExcludedBP, height(BandPowerT));

    % ---- NEW: apply the shared noisy-subject exclusion list (see
    % EXCLUSION_CONFIG_cognitive.m), only when RUN_MODE=sensitivity_excl_noisy.
    % Never applied to the full_sample run. Applied to the X04 side here;
    % the X03 side (taskEntriesENR) is excluded identically just below, so
    % Section 4's recCombos/taskEntriesENR count cross-check stays valid. ----
    nExcludedRowsSensitivityBP = 0;
    if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
        subjNumForExclBP = local_subject_numeric(BandPowerT.Subject);
        excludeMaskBP = ismember(subjNumForExclBP, EXCLUDE_NOISY_SUBJECTS);
        nExcludedRowsSensitivityBP = sum(excludeMaskBP);
        BandPowerT = BandPowerT(~excludeMaskBP, :);
        fprintf('Task=%s: RUN_MODE=%s -- excluded %d X04 channel-level row(s) belonging to subject(s) {%s}; %d rows remain.\n', ...
            TASK_LABEL, RUN_MODE_LABEL, nExcludedRowsSensitivityBP, ...
            strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','), height(BandPowerT));
    end

    if isempty(BandPowerT)
        warning('Task=%s: no E/N/R rows found -- skipping this Task entirely.', TASK_LABEL);
        continue;
    end

    %% Section 3: filter psd_data to this Task, then to Condition in {E,N,R} (X03 side)
    taskMaskPSD = strcmp({psd_data.TaskType}, TASK_KEY);
    taskEntries = psd_data(taskMaskPSD);
    nBaselineExcludedPSD = sum(strcmp({taskEntries.Condition}, 'Baseline'));
    condMaskPSD = ismember({taskEntries.Condition}, CONDITION_ORDER);
    taskEntriesENR = taskEntries(condMaskPSD);
    fprintf('Task=%s: %d recordings loaded from X03, %d Baseline recordings excluded, %d E/N/R recordings retained.\n', ...
        TASK_LABEL, numel(taskEntries), nBaselineExcludedPSD, numel(taskEntriesENR));

    % ---- NEW: same shared noisy-subject exclusion list, applied to the
    % X03 side so it stays consistent with the X04 side above. ----
    nExcludedRowsSensitivityPSD = 0;
    if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
        subjNumForExclPSD = local_subject_numeric({taskEntriesENR.Subject});
        excludeMaskPSD = ismember(subjNumForExclPSD, EXCLUDE_NOISY_SUBJECTS);
        nExcludedRowsSensitivityPSD = sum(excludeMaskPSD);
        taskEntriesENR = taskEntriesENR(~excludeMaskPSD);
        fprintf('Task=%s: RUN_MODE=%s -- excluded %d X03 recording(s) belonging to subject(s) {%s}; %d recordings remain.\n', ...
            TASK_LABEL, RUN_MODE_LABEL, nExcludedRowsSensitivityPSD, ...
            strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','), numel(taskEntriesENR));
    end

    recKeyMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for e = 1:numel(taskEntriesENR)
        k = sprintf('%s|%g|%s', taskEntriesENR(e).Subject, taskEntriesENR(e).BlockNumber, taskEntriesENR(e).Condition);
        if isKey(recKeyMap, k)
            integrityIssues{end+1} = sprintf('Task=%s: duplicate psd_data entry for key %s.', TASK_LABEL, k); %#ok<AGROW>
        end
        recKeyMap(k) = e;
    end

    %% Section 4: recCombos (Subject x BlockNumber x Condition), from X04
    %  side, cross-checked for count consistency against X03 side
    recCombos = unique(BandPowerT(:, {'Subject','BlockNumber','Condition'}), 'rows');
    nRecordings = height(recCombos);
    if nRecordings ~= numel(taskEntriesENR)
        integrityIssues{end+1} = sprintf('Task=%s: X04 recCombos (%d) and X03 E/N/R recordings (%d) counts differ.', ...
            TASK_LABEL, nRecordings, numel(taskEntriesENR)); %#ok<AGROW>
    end

    %% ====================================================================
    %  SECTION 5: ELECTRODE-LEVEL median-robust metrics -- X04's own
    %  WindowMedian_Absolute/WindowMedian_Relative columns used DIRECTLY
    %  (median-across-valid-epochs, per channel; no new computation
    %  needed at this level -- see header rationale).
    % ====================================================================
    MedianRobustElectrodeT = BandPowerT(:, {'Subject','TaskType','BlockNumber','Condition','Channel','Band'});
    MedianRobustElectrodeT.AbsMedian   = BandPowerT.WindowMedian_Absolute;
    MedianRobustElectrodeT.RelMedian   = BandPowerT.WindowMedian_Relative;
    MedianRobustElectrodeT.CV_Absolute = BandPowerT.WindowSD_Absolute ./ BandPowerT.AbsolutePower;
    MedianRobustElectrodeT.CV_Relative = BandPowerT.WindowSD_Relative ./ BandPowerT.RelativePower;
    nBadCV = sum(~isfinite(MedianRobustElectrodeT.CV_Absolute)) + sum(~isfinite(MedianRobustElectrodeT.CV_Relative));
    if nBadCV > 0
        integrityIssues{end+1} = sprintf(['Task=%s: %d non-finite CV_Absolute/CV_Relative value(s) found ' ...
            '(non-positive AbsolutePower/RelativePower denominator) -- left as-is, not silently patched.'], TASK_LABEL, nBadCV); %#ok<AGROW>
    end
    fprintf('Electrode-level MedianRobustElectrodeT built: %d rows (%d recordings x %d channels x %d bands).\n', ...
        height(MedianRobustElectrodeT), nRecordings, nChanRef, nBands);

    %% ====================================================================
    %  SECTION 6: ROI-LEVEL median-robust metrics -- median-across-epochs
    %  FIRST (Section 5's per-channel AbsMedian/RelMedian/CV_Absolute/
    %  CV_Relative columns), THEN arithmetic mean across each ROI's
    %  member channels (see header rationale).
    % ====================================================================
    elecChannelUC = upper(MedianRobustElectrodeT.Channel);   % local var only -- never added as a table column
    roiRows9 = cell(0, 9);
    roiCols9 = {'Subject','BlockNumber','Condition','ROI','Band','AbsMedian','RelMedian','CV_Absolute','CV_Relative'};

    for i = 1:height(recCombos)
        subjStr = recCombos.Subject{i}; blockNum = recCombos.BlockNumber(i); condStr = recCombos.Condition{i};
        recMask = strcmp(MedianRobustElectrodeT.Subject, subjStr) & MedianRobustElectrodeT.BlockNumber == blockNum & strcmp(MedianRobustElectrodeT.Condition, condStr);
        for r = 1:nROI
            roiName = ROI_NAMES{r};
            roiChansUC = upper(ROI_DEFINITION.(roiName));
            chanMask = ismember(elecChannelUC, roiChansUC);
            for b = 1:nBands
                bandName = ALL_BANDS{b};
                rowMask = recMask & chanMask & strcmp(MedianRobustElectrodeT.Band, bandName);
                nFound = sum(rowMask);
                if nFound ~= numel(roiChansUC)
                    integrityIssues{end+1} = sprintf('Sub%s Block%g Cond%s ROI=%s Band=%s (Task=%s): found %d/%d expected ROI channels.', ...
                        subjStr, blockNum, condStr, roiName, bandName, TASK_LABEL, nFound, numel(roiChansUC)); %#ok<AGROW>
                end
                roiRows9(end+1, :) = {subjStr, blockNum, condStr, roiName, bandName, ...
                    mean(MedianRobustElectrodeT.AbsMedian(rowMask)), mean(MedianRobustElectrodeT.RelMedian(rowMask)), ...
                    mean(MedianRobustElectrodeT.CV_Absolute(rowMask)), mean(MedianRobustElectrodeT.CV_Relative(rowMask))}; %#ok<AGROW>
            end
        end
    end
    MedianRobustROIT = local_build_table(roiRows9, [1 3 4 5], [2 6 7 8 9], roiCols9);
    MedianRobustROIT.TaskType = repmat({TASK_KEY}, height(MedianRobustROIT), 1);
    MedianRobustROIT = MedianRobustROIT(:, {'Subject','TaskType','BlockNumber','Condition','ROI','Band','AbsMedian','RelMedian','CV_Absolute','CV_Relative'});
    fprintf('ROI-level MedianRobustROIT built: %d rows (%d recordings x %d ROIs x %d bands).\n', ...
        height(MedianRobustROIT), nRecordings, nROI, nBands);

    %% ====================================================================
    %  SECTION 7: GLOBAL-LEVEL median-robust metrics + BROADBAND
    %  TotalPower (Electrode + Global), computed FRESH from X03's
    %  psd_data (never touched by X04). Global order-of-operations is
    %  mean-across-channels PER EPOCH first, THEN median-across-epochs --
    %  the OPPOSITE order from ROI (see header rationale).
    % ====================================================================
    globalRows = cell(0, 8);
    globalCols = {'Subject','BlockNumber','Condition','Band','AbsMedian','RelMedian','CV_Absolute','CV_Relative'};
    totalElecRows = cell(0, 6);
    totalElecCols = {'Subject','BlockNumber','Condition','Channel','TotalPowerMean','TotalPowerMedian'};
    totalGlobalRows = cell(0, 5);
    totalGlobalCols = {'Subject','BlockNumber','Condition','TotalPowerMean','TotalPowerMedian'};

    nZeroValidSkipped = 0;
    for i = 1:height(recCombos)
        subjStr = recCombos.Subject{i}; blockNum = recCombos.BlockNumber(i); condStr = recCombos.Condition{i};
        k = sprintf('%s|%g|%s', subjStr, blockNum, condStr);
        if ~isKey(recKeyMap, k)
            integrityIssues{end+1} = sprintf('Sub%s Block%g Cond%s (Task=%s): no matching X03 psd_data entry found -- Global/TotalPower rows skipped.', ...
                subjStr, blockNum, condStr, TASK_LABEL); %#ok<AGROW>
            continue;
        end
        entry = taskEntriesENR(recKeyMap(k));
        valid = entry.WindowValid;
        nValidThis = sum(valid);
        if nValidThis == 0
            nZeroValidSkipped = nZeroValidSkipped + 1;
            integrityIssues{end+1} = sprintf('Sub%s Block%g Cond%s (Task=%s): 0 valid windows -- Global/TotalPower rows skipped.', ...
                subjStr, blockNum, condStr, TASK_LABEL); %#ok<AGROW>
            continue;
        end
        if ~isequal(entry.ChannelLabels(:), ChannelLabelsReference)
            integrityIssues{end+1} = sprintf('Sub%s Block%g Cond%s (Task=%s): ChannelLabels differ from the reference channel order.', ...
                subjStr, blockNum, condStr, TASK_LABEL); %#ok<AGROW>
        end

        % --- broadband TotalPower (no Band dimension) ---
        totMat = entry.TotalPower_1to45Hz(valid, :);               % nValid x nChan
        perChanMean   = mean(totMat, 1);
        perChanMedian = median(totMat, 1);
        for ch = 1:nChanRef
            totalElecRows(end+1, :) = {subjStr, blockNum, condStr, ChannelLabelsReference{ch}, perChanMean(ch), perChanMedian(ch)}; %#ok<AGROW>
        end
        globalTotEpoch = mean(totMat, 2);                            % nValid x 1, mean across channels PER epoch
        totalGlobalRows(end+1, :) = {subjStr, blockNum, condStr, mean(globalTotEpoch), median(globalTotEpoch)}; %#ok<AGROW>

        % --- band-dependent Global median-robust metrics ---
        for b = 1:nBands
            bandName = ALL_BANDS{b};
            absMat = entry.AbsoluteBandPower.(bandName)(valid, :);   % nValid x nChan
            relMat = entry.RelativeBandPower.(bandName)(valid, :);
            globalAbsEpoch = mean(absMat, 2);    % mean across channels PER epoch
            globalRelEpoch = mean(relMat, 2);
            absMedG = median(globalAbsEpoch);
            relMedG = median(globalRelEpoch);
            cvAbsG = std(globalAbsEpoch) / mean(globalAbsEpoch);
            cvRelG = std(globalRelEpoch) / mean(globalRelEpoch);
            globalRows(end+1, :) = {subjStr, blockNum, condStr, bandName, absMedG, relMedG, cvAbsG, cvRelG}; %#ok<AGROW>
        end
    end

    MedianRobustGlobalT = local_build_table(globalRows, [1 3 4], [2 5 6 7 8], globalCols);
    MedianRobustGlobalT.TaskType = repmat({TASK_KEY}, height(MedianRobustGlobalT), 1);
    MedianRobustGlobalT = MedianRobustGlobalT(:, {'Subject','TaskType','BlockNumber','Condition','Band','AbsMedian','RelMedian','CV_Absolute','CV_Relative'});

    TotalPowerElectrodeT = local_build_table(totalElecRows, [1 3 4], [2 5 6], totalElecCols);
    TotalPowerElectrodeT.TaskType = repmat({TASK_KEY}, height(TotalPowerElectrodeT), 1);
    TotalPowerElectrodeT = TotalPowerElectrodeT(:, {'Subject','TaskType','BlockNumber','Condition','Channel','TotalPowerMean','TotalPowerMedian'});

    TotalPowerGlobalT = local_build_table(totalGlobalRows, [1 3], [2 4 5], totalGlobalCols);
    TotalPowerGlobalT.TaskType = repmat({TASK_KEY}, height(TotalPowerGlobalT), 1);
    TotalPowerGlobalT = TotalPowerGlobalT(:, {'Subject','TaskType','BlockNumber','Condition','TotalPowerMean','TotalPowerMedian'});

    fprintf(['Global-level MedianRobustGlobalT: %d rows. TotalPowerElectrodeT: %d rows. TotalPowerGlobalT: %d rows. ' ...
        '(%d recordings skipped for 0 valid windows or missing X03 entry)\n'], ...
        height(MedianRobustGlobalT), height(TotalPowerElectrodeT), height(TotalPowerGlobalT), nZeroValidSkipped);

    %% ====================================================================
    %  SECTION 8: TotalPowerROIT -- per-channel TotalPowerMean/
    %  TotalPowerMedian (Section 7) averaged across each ROI's member
    %  channels (same order-of-operations rationale as Section 6:
    %  epoch-level robust/mean statistic computed per channel FIRST, then
    %  spatial mean across the ROI).
    % ====================================================================
    totElecChannelUC = upper(TotalPowerElectrodeT.Channel);   % local var only -- never added as a table column
    totalRoiRows6 = cell(0, 6);
    totalRoiCols6 = {'Subject','BlockNumber','Condition','ROI','TotalPowerMean','TotalPowerMedian'};
    for i = 1:height(recCombos)
        subjStr = recCombos.Subject{i}; blockNum = recCombos.BlockNumber(i); condStr = recCombos.Condition{i};
        recMask = strcmp(TotalPowerElectrodeT.Subject, subjStr) & TotalPowerElectrodeT.BlockNumber == blockNum & strcmp(TotalPowerElectrodeT.Condition, condStr);
        if ~any(recMask); continue; end
        for r = 1:nROI
            roiName = ROI_NAMES{r};
            roiChansUC = upper(ROI_DEFINITION.(roiName));
            chanMask = ismember(totElecChannelUC, roiChansUC);
            rowMask = recMask & chanMask;
            totalRoiRows6(end+1, :) = {subjStr, blockNum, condStr, roiName, ...
                mean(TotalPowerElectrodeT.TotalPowerMean(rowMask)), mean(TotalPowerElectrodeT.TotalPowerMedian(rowMask))}; %#ok<AGROW>
        end
    end
    TotalPowerROIT = local_build_table(totalRoiRows6, [1 3 4], [2 5 6], totalRoiCols6);
    TotalPowerROIT.TaskType = repmat({TASK_KEY}, height(TotalPowerROIT), 1);
    TotalPowerROIT = TotalPowerROIT(:, {'Subject','TaskType','BlockNumber','Condition','ROI','TotalPowerMean','TotalPowerMedian'});
    fprintf('TotalPowerROIT built: %d rows (%d recordings x %d ROIs).\n', height(TotalPowerROIT), nRecordings, nROI);

    %% ====================================================================
    %  SECTION 9: SAVE THE OUTPUT-CONTRACT FILE -- exactly these 6 tables,
    %  exactly these names/columns (X12/X13 depend on this exact schema).
    % ====================================================================
    save(fullfile(STATS_DIR, 'median_robust_data.mat'), ...
        'MedianRobustElectrodeT', 'MedianRobustROIT', 'MedianRobustGlobalT', ...
        'TotalPowerElectrodeT', 'TotalPowerROIT', 'TotalPowerGlobalT');

    %% ====================================================================
    %  SECTION 10: BAND-DEPENDENT OMNIBUS TESTS
    %  Electrode/ROI: family = per (Band x Metric), across units.
    %  Global       : family = per (Tier x Metric), across bands in tier.
    % ====================================================================
    levelConfigsBand = struct('Name', {'Electrode', 'ROI'}, ...
        'Table', {MedianRobustElectrodeT, MedianRobustROIT}, ...
        'UnitCol', {'Channel', 'ROI'}, ...
        'Units', {ChannelLabelsReference, ROI_NAMES(:)});

    omnibusBandRows = cell(0, 12);
    omnibusBandCols = {'Level','Unit','Metric','Band','Tier','FStat','DF1','DF2','pValue','pFDR','SignificantFDR','FamilySize'};
    lmeMapBand = containers.Map();   % key: 'Level|Unit|Metric|Band'

    for lc = 1:numel(levelConfigsBand)
        levelName = levelConfigsBand(lc).Name;
        Tlevel = levelConfigsBand(lc).Table;
        unitCol = levelConfigsBand(lc).UnitCol;
        units = levelConfigsBand(lc).Units;
        for mI = 1:numel(METRICS_BAND)
            metricName = METRICS_BAND{mI};
            for b = 1:nBands
                bandName = ALL_BANDS{b};
                pvalsThisFamily = nan(numel(units), 1);
                fitResultsThisFamily = cell(numel(units), 1);
                for u = 1:numel(units)
                    unitName = units{u};
                    subT = Tlevel(strcmp(Tlevel.(unitCol), unitName) & strcmp(Tlevel.Band, bandName), :);
                    fr = local_fit_condition_omnibus(subT, metricName, CONDITION_ORDER);
                    fitResultsThisFamily{u} = fr;
                    if fr.ok
                        pvalsThisFamily(u) = fr.p;
                        lmeMapBand(sprintf('%s|%s|%s|%s', levelName, unitName, metricName, bandName)) = fr.lme;
                    else
                        integrityIssues{end+1} = sprintf('%s: fitlme FAILED for %s / %s / %s (Task=%s): %s', ...
                            levelName, unitName, metricName, bandName, TASK_LABEL, fr.errMsg); %#ok<AGROW>
                    end
                end
                okMask = ~isnan(pvalsThisFamily);
                nInFamily = sum(okMask);
                if nInFamily ~= numel(units)
                    integrityIssues{end+1} = sprintf('%s FDR family %s / %s (Task=%s) has %d/%d successful fits.', ...
                        levelName, metricName, bandName, TASK_LABEL, nInFamily, numel(units)); %#ok<AGROW>
                end
                if nInFamily == 0; continue; end
                [adjP, reject] = local_bh_fdr_or_raw(pvalsThisFamily(okMask), FDR_ALPHA);
                okIdx = find(okMask);
                for jj = 1:numel(okIdx)
                    u = okIdx(jj);
                    fr = fitResultsThisFamily{u};
                    omnibusBandRows(end+1, :) = {levelName, units{u}, metricName, bandName, tierOf(bandName), ...
                        fr.F, fr.DF1, fr.DF2, fr.p, adjP(jj), reject(jj), nInFamily}; %#ok<AGROW>
                end
            end
        end
    end

    % Global level -- tiered family across bands within tier, per metric
    for mI = 1:numel(METRICS_BAND)
        metricName = METRICS_BAND{mI};
        for t = 1:size(tierDefs, 1)
            tierName = tierDefs{t, 1};
            tierBands = tierDefs{t, 2};
            pvalsTier = nan(numel(tierBands), 1);
            frTier = cell(numel(tierBands), 1);
            for bi = 1:numel(tierBands)
                bandName = tierBands{bi};
                subT = MedianRobustGlobalT(strcmp(MedianRobustGlobalT.Band, bandName), :);
                fr = local_fit_condition_omnibus(subT, metricName, CONDITION_ORDER);
                frTier{bi} = fr;
                if fr.ok
                    pvalsTier(bi) = fr.p;
                    lmeMapBand(sprintf('Global|Global|%s|%s', metricName, bandName)) = fr.lme;
                else
                    integrityIssues{end+1} = sprintf('Global: fitlme FAILED for %s / %s (Task=%s): %s', ...
                        metricName, bandName, TASK_LABEL, fr.errMsg); %#ok<AGROW>
                end
            end
            okMask = ~isnan(pvalsTier);
            nInFamily = sum(okMask);
            if nInFamily == 0; continue; end
            [adjP, reject] = local_bh_fdr_or_raw(pvalsTier(okMask), FDR_ALPHA);
            okIdx = find(okMask);
            for jj = 1:numel(okIdx)
                bi = okIdx(jj);
                fr = frTier{bi};
                omnibusBandRows(end+1, :) = {'Global', 'Global', metricName, tierBands{bi}, tierName, ...
                    fr.F, fr.DF1, fr.DF2, fr.p, adjP(jj), reject(jj), nInFamily}; %#ok<AGROW>
            end
        end
    end

    OmnibusBandT = local_build_table(omnibusBandRows, [1 2 3 4 5], [6 7 8 9 10 11 12], omnibusBandCols);
    fprintf('Band-dependent omnibus tests complete (Task=%s): %d rows (Electrode+ROI+Global).\n', TASK_LABEL, height(OmnibusBandT));

    %% ====================================================================
    %  SECTION 11: BROADBAND TotalPower OMNIBUS TESTS
    %  Electrode/ROI: family = per Metric, across units (no Band/Tier axis).
    %  Global       : single test per metric (family size 1, reported as-is).
    % ====================================================================
    levelConfigsTotal = struct('Name', {'Electrode', 'ROI'}, ...
        'Table', {TotalPowerElectrodeT, TotalPowerROIT}, ...
        'UnitCol', {'Channel', 'ROI'}, ...
        'Units', {ChannelLabelsReference, ROI_NAMES(:)});

    omnibusTotalRows = cell(0, 10);
    omnibusTotalCols = {'Level','Unit','Metric','FStat','DF1','DF2','pValue','pFDR','SignificantFDR','FamilySize'};
    lmeMapTotal = containers.Map();   % key: 'Level|Unit|Metric'

    for lc = 1:numel(levelConfigsTotal)
        levelName = levelConfigsTotal(lc).Name;
        Tlevel = levelConfigsTotal(lc).Table;
        unitCol = levelConfigsTotal(lc).UnitCol;
        units = levelConfigsTotal(lc).Units;
        for mI = 1:numel(METRICS_TOTAL)
            metricName = METRICS_TOTAL{mI};
            pvalsThisFamily = nan(numel(units), 1);
            fitResultsThisFamily = cell(numel(units), 1);
            for u = 1:numel(units)
                unitName = units{u};
                subT = Tlevel(strcmp(Tlevel.(unitCol), unitName), :);
                fr = local_fit_condition_omnibus(subT, metricName, CONDITION_ORDER);
                fitResultsThisFamily{u} = fr;
                if fr.ok
                    pvalsThisFamily(u) = fr.p;
                    lmeMapTotal(sprintf('%s|%s|%s', levelName, unitName, metricName)) = fr.lme;
                else
                    integrityIssues{end+1} = sprintf('%s: fitlme FAILED for TotalPower %s / %s (Task=%s): %s', ...
                        levelName, unitName, metricName, TASK_LABEL, fr.errMsg); %#ok<AGROW>
                end
            end
            okMask = ~isnan(pvalsThisFamily);
            nInFamily = sum(okMask);
            if nInFamily == 0; continue; end
            [adjP, reject] = local_bh_fdr_or_raw(pvalsThisFamily(okMask), FDR_ALPHA);
            okIdx = find(okMask);
            for jj = 1:numel(okIdx)
                u = okIdx(jj);
                fr = fitResultsThisFamily{u};
                omnibusTotalRows(end+1, :) = {levelName, units{u}, metricName, fr.F, fr.DF1, fr.DF2, fr.p, adjP(jj), reject(jj), nInFamily}; %#ok<AGROW>
            end
        end
    end

    for mI = 1:numel(METRICS_TOTAL)
        metricName = METRICS_TOTAL{mI};
        fr = local_fit_condition_omnibus(TotalPowerGlobalT, metricName, CONDITION_ORDER);
        if ~fr.ok
            integrityIssues{end+1} = sprintf('Global: fitlme FAILED for TotalPower %s (Task=%s): %s', metricName, TASK_LABEL, fr.errMsg); %#ok<AGROW>
            continue;
        end
        lmeMapTotal(sprintf('Global|Global|%s', metricName)) = fr.lme;
        omnibusTotalRows(end+1, :) = {'Global', 'Global', metricName, fr.F, fr.DF1, fr.DF2, fr.p, fr.p, fr.p <= FDR_ALPHA, 1}; %#ok<AGROW>
    end

    OmnibusTotalPowerT = local_build_table(omnibusTotalRows, [1 2 3], [4 5 6 7 8 9 10], omnibusTotalCols);
    fprintf('Broadband TotalPower omnibus tests complete (Task=%s): %d rows.\n', TASK_LABEL, height(OmnibusTotalPowerT));

    %% ====================================================================
    %  SECTION 12: PAIRWISE CONTRASTS (N-E, R-E, R-N; Holm-corrected,
    %  family=3; gated on each row's own FDR-family survival)
    % ====================================================================
    pairwiseBandRows = cell(0, 13);
    pairwiseBandCols = {'Level','Unit','Metric','Band','Contrast','Estimate','SE','DF','tStat','pValue','pHolm','CI_Low','CI_High'};
    sigMaskBand = OmnibusBandT.SignificantFDR;
    for row = find(sigMaskBand)'
        levelName = OmnibusBandT.Level{row}; unitName = OmnibusBandT.Unit{row};
        metricName = OmnibusBandT.Metric{row}; bandName = OmnibusBandT.Band{row};
        key = sprintf('%s|%s|%s|%s', levelName, unitName, metricName, bandName);
        if ~isKey(lmeMapBand, key); continue; end
        lme = lmeMapBand(key);
        rows3 = local_run_pairwise(lme, REF_LEVEL, CONTRASTS, FDR_ALPHA);
        for c = 1:3
            pairwiseBandRows(end+1, :) = [{levelName, unitName, metricName, bandName}, rows3(c, :)]; %#ok<AGROW>
        end
    end
    PairwiseBandT = local_build_table(pairwiseBandRows, [1 2 3 4 5], [6 7 8 9 10 11 12 13], pairwiseBandCols);

    pairwiseTotalRows = cell(0, 12);
    pairwiseTotalCols = {'Level','Unit','Metric','Contrast','Estimate','SE','DF','tStat','pValue','pHolm','CI_Low','CI_High'};
    sigMaskTotal = OmnibusTotalPowerT.SignificantFDR;
    for row = find(sigMaskTotal)'
        levelName = OmnibusTotalPowerT.Level{row}; unitName = OmnibusTotalPowerT.Unit{row};
        metricName = OmnibusTotalPowerT.Metric{row};
        key = sprintf('%s|%s|%s', levelName, unitName, metricName);
        if ~isKey(lmeMapTotal, key); continue; end
        lme = lmeMapTotal(key);
        rows3 = local_run_pairwise(lme, REF_LEVEL, CONTRASTS, FDR_ALPHA);
        for c = 1:3
            pairwiseTotalRows(end+1, :) = [{levelName, unitName, metricName}, rows3(c, :)]; %#ok<AGROW>
        end
    end
    PairwiseTotalPowerT = local_build_table(pairwiseTotalRows, [1 2 3 4], [5 6 7 8 9 10 11 12], pairwiseTotalCols);

    fprintf('Pairwise contrasts complete (Task=%s): %d band-metric rows, %d TotalPower rows.\n', ...
        TASK_LABEL, height(PairwiseBandT), height(PairwiseTotalPowerT));

    %% ====================================================================
    %  SECTION 13: RESIDUAL DIAGNOSTICS -- Global level only, Primary +
    %  Secondary tier bands, both band-dependent metrics (judgment call
    %  #3, see header).
    % ====================================================================
    diagBands = [PRIMARY_TIER_BANDS, SECONDARY_TIER_BANDS];
    nDiagSaved = 0;
    for mI = 1:numel(METRICS_BAND)
        metricName = METRICS_BAND{mI};
        for bi = 1:numel(diagBands)
            bandName = diagBands{bi};
            key = sprintf('Global|Global|%s|%s', metricName, bandName);
            if ~isKey(lmeMapBand, key); continue; end
            lme = lmeMapBand(key);
            try
                res = residuals(lme, 'ResidualType', 'Pearson');
                fig = figure('Visible', 'off', 'Position', [100 100 900 400]);
                subplot(1,2,1);
                histogram(res, 15);
                title(sprintf('Residuals: Global %s (%s), Task=%s', bandName, metricName, TASK_LABEL), 'Interpreter', 'none');
                xlabel('Pearson residual'); ylabel('Count');
                subplot(1,2,2);
                qqplot(res);
                title('Normal Q-Q');
                outPng = fullfile(STATS_DIR, sprintf('diagnostic_residuals_Global_%s_%s.png', bandName, metricName));
                saveas(fig, outPng);
                close(fig);
                nDiagSaved = nDiagSaved + 1;
            catch ME
                integrityIssues{end+1} = sprintf('Residual diagnostic failed for Global %s / %s (Task=%s): %s', ...
                    metricName, bandName, TASK_LABEL, ME.message); %#ok<AGROW>
            end
        end
    end
    fprintf('Saved %d residual diagnostic plots (Task=%s; Global level only, Primary+Secondary bands -- judgment call, see header/PARAMS11).\n', ...
        nDiagSaved, TASK_LABEL);

    %% ====================================================================
    %  SECTION 14: PARAMETERS AND ADDITIONAL OUTPUT
    % ====================================================================
    PARAMS11 = struct();
    PARAMS11.ScriptName = 'X11_median_robust_family_cognitive.m';
    PARAMS11.Task = TASK_LABEL;
    PARAMS11.GeneratedOn = datestr(now); %#ok<TNOW1,DATST>
    PARAMS11.Purpose = ['Exploratory addition, alongside (never replacing or ' ...
        'pooling with) the locked co-primary LogAbsolutePower/RelativePower metrics of X05/X06/X08/X09/X10.'];
    PARAMS11.ElectrodeOrderOfOperations = ['Electrode: X04''s own WindowMedian_Absolute/WindowMedian_Relative columns ' ...
        'used directly (median-across-valid-epochs, per channel); no new computation.'];
    PARAMS11.ROIOrderOfOperations = ['ROI: median-across-epochs FIRST (per channel, X04''s WindowMedian columns), THEN ' ...
        'arithmetic mean across the ROI''s member channels -- protects each channel''s own temporal-outlier problem ' ...
        'before spatial pooling, since ROIs here have only 6-23 member channels.'];
    PARAMS11.GlobalOrderOfOperations = ['Global: mean-across-channels PER EPOCH first (computed fresh from X03 psd_data, ' ...
        'across all 62 reference channels), THEN median-across-epochs of that per-epoch global series -- the OPPOSITE ' ...
        'order from ROI, because at 62-channel scale the spatial-averaging step already smooths spatial noise via the ' ...
        'large channel count, so the robust step is applied on the temporal axis AFTER that.'];
    PARAMS11.Metrics_BandDependent = METRICS_BAND;
    PARAMS11.Metrics_Broadband = METRICS_TOTAL;
    PARAMS11.Model = 'LinearMixedEffects_RandomIntercept';
    PARAMS11.FixedEffectsPrimary = 'Condition (3-level: E, N, R; reference = E)';
    PARAMS11.RandomEffects = '(1|Subject)';
    PARAMS11.FitMethod = 'REML';
    PARAMS11.DFMethod = 'Satterthwaite';
    PARAMS11.FDRAlpha = FDR_ALPHA;
    PARAMS11.FDRFamilyScheme = ['Electrode/ROI levels: separate BH-FDR family per (Band x Metric), across units ' ...
        '(family size = up to 62 electrodes / up to 5 ROIs). Global level: separate BH-FDR family per (Tier x Metric), ' ...
        'across the bands in that tier (family size = 1/3/3 for Primary/Secondary/Exploratory), matching X05''s tiered ' ...
        'scheme -- used as the actual FDR-family axis only at Global level, where there is no spatial axis to pool ' ...
        'across instead. Broadband TotalPower metrics (no Band dimension): Electrode/ROI = single family per metric ' ...
        'across units (up to 62/5); Global = single test per metric (family size 1, reported as-is, never treated as ' ...
        'pre-corrected). Metrics, levels and bands/tiers are NEVER pooled into one family, and NEVER pooled with the ' ...
        'locked LogAbsolutePower/RelativePower families of X05/X06/X08/X09.'];
    PARAMS11.PairwiseContrastMethod = ['Holm-Bonferroni across the 3 pairwise contrasts (N-E, R-E, R-N -- matching ' ...
        'X05/X06''s convention, standardized here since this is a new family not tied to any old script), computed ' ...
        'only for Level x Unit x Metric x Band (or x Metric for TotalPower) combinations whose own FDR family test ' ...
        'survives BH-FDR.'];
    PARAMS11.BaselineExclusionNote = sprintf(['%d channel-level (X04) / %d recording-level (X03) Baseline row(s) ' ...
        'excluded before any modeling for Task=%s.'], nBaselineExcludedBP, nBaselineExcludedPSD, TASK_LABEL);
    PARAMS11.JudgmentCalls = { ...
        ['ROI-level CV_Absolute/CV_Relative = arithmetic mean of per-channel CV across ROI member channels -- a ' ...
         'reasonable summary, NOT literally the CV of the ROI-aggregated series.'], ...
        ['Global-level CV_Absolute/CV_Relative = SD-across-epochs of the per-epoch mean-across-channels series, ' ...
         'divided by its mean.'], ...
        ['Residual diagnostics scoped to Global level only (Primary+Secondary tier bands, both band-dependent ' ...
         'metrics) to keep plot count tractable -- Electrode (up to 62 x 4 metrics x 4 bands) and ROI (5 x 4 x 4) ' ...
         'diagnostics are not generated by this script.'], ...
        ['A recording with 0 valid windows (X03 WindowValid) is skipped for Global/TotalPower rows (logged in ' ...
         'integrityIssues), exactly like X04''s own "skip aggregation, still record QC" pattern for such recordings; ' ...
         'no Global/TotalPower row is fabricated for it.']};
    PARAMS11.SourcePhase2PSDFile = psdMatFile;
    PARAMS11.SourcePhase3BandPowerFile = bandPowerFile;
    PARAMS11.SourcePhase5DROIFile = roiDefFile;
    PARAMS11.SourcePhase5DDecision = PARAMS5D.ROIDecision;
    PARAMS11.DesignVerification = struct('nRecordings', nRecordings, 'nBaselineExcludedBP', nBaselineExcludedBP, ...
        'nBaselineExcludedPSD', nBaselineExcludedPSD, 'nZeroValidSkipped', nZeroValidSkipped, 'IntegrityIssues', {integrityIssues});
    PARAMS11.RunDateTime = datestr(now); %#ok<TNOW1,DATST>
    PARAMS11.RunMode = RUN_MODE;
    if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
        PARAMS11.ExcludedSubjects = EXCLUDE_NOISY_SUBJECTS;
        PARAMS11.ExcludedSubjectsRationale = EXCLUSION_RATIONALE_cognitive;
        PARAMS11.RowsExcludedForSensitivity = struct('BandPowerRows', nExcludedRowsSensitivityBP, 'PSDRecordings', nExcludedRowsSensitivityPSD);
    else
        PARAMS11.ExcludedSubjects = [];
    end

    save(fullfile(STATS_DIR, 'phase11_parameters.mat'), 'PARAMS11');
    save(fullfile(STATS_DIR, 'phase11_omnibus_band_results.mat'), 'OmnibusBandT');
    save(fullfile(STATS_DIR, 'phase11_omnibus_totalpower_results.mat'), 'OmnibusTotalPowerT');
    save(fullfile(STATS_DIR, 'phase11_pairwise_band_contrasts.mat'), 'PairwiseBandT');
    save(fullfile(STATS_DIR, 'phase11_pairwise_totalpower_contrasts.mat'), 'PairwiseTotalPowerT');
    writetable(OmnibusBandT, fullfile(STATS_DIR, 'phase11_omnibus_band_results.csv'));
    writetable(OmnibusTotalPowerT, fullfile(STATS_DIR, 'phase11_omnibus_totalpower_results.csv'));
    if ~isempty(PairwiseBandT)
        writetable(PairwiseBandT, fullfile(STATS_DIR, 'phase11_pairwise_band_contrasts.csv'));
    end
    if ~isempty(PairwiseTotalPowerT)
        writetable(PairwiseTotalPowerT, fullfile(STATS_DIR, 'phase11_pairwise_totalpower_contrasts.csv'));
    end

    fprintf('\nSaved output to: %s\n', STATS_DIR);

    %% ====================================================================
    %  SECTION 15: REPORT
    % ====================================================================
    reportLines = {};
    reportLines{end+1} = repmat('=', 1, 60);
    reportLines{end+1} = sprintf('PHASE 11 (X11) -- MEDIAN-ROBUST DERIVED-METRIC FAMILY REPORT (Task=%s)', TASK_LABEL);
    reportLines{end+1} = repmat('=', 1, 60);
    reportLines{end+1} = '';
    reportLines{end+1} = 'Exploratory ADDITION alongside (never replacing/pooling with) the locked co-primary';
    reportLines{end+1} = 'LogAbsolutePower/RelativePower metrics of X05/X06/X08/X09/X10.';
    reportLines{end+1} = '';
    reportLines{end+1} = sprintf('Baseline rows excluded: %d channel-level (X04) / %d recording-level (X03)', nBaselineExcludedBP, nBaselineExcludedPSD);
    reportLines{end+1} = sprintf('Recordings: %d | Reference channels: %d | ROIs: %d | Bands: %d | 0-valid-window recordings skipped: %d', ...
        nRecordings, nChanRef, nROI, nBands, nZeroValidSkipped);
    reportLines{end+1} = '';
    reportLines{end+1} = 'Order-of-operations (see script header for full rationale):';
    reportLines{end+1} = '  Electrode: median-across-valid-epochs per channel (X04 WindowMedian columns, reused directly).';
    reportLines{end+1} = '  ROI      : median-across-epochs per channel FIRST, then mean across ROI member channels.';
    reportLines{end+1} = '  Global   : mean-across-channels PER EPOCH first, then median-across-epochs of that series.';
    reportLines{end+1} = '';
    reportLines{end+1} = 'FDR family sizes used:';
    reportLines{end+1} = '  Band-dependent metrics (AbsMedian, RelMedian, CV_Absolute, CV_Relative -- each its own family):';
    reportLines{end+1} = '    Electrode: family = per (Band x Metric), across electrodes (up to 62).';
    reportLines{end+1} = '    ROI      : family = per (Band x Metric), across ROIs (up to 5).';
    reportLines{end+1} = '    Global   : family = per (Tier x Metric), across bands in that tier (Primary=1, Secondary=3, Exploratory=3).';
    reportLines{end+1} = '  Broadband TotalPower metrics (TotalPowerMean, TotalPowerMedian -- each its own family):';
    reportLines{end+1} = '    Electrode: family = per Metric, across electrodes (up to 62).';
    reportLines{end+1} = '    ROI      : family = per Metric, across ROIs (up to 5).';
    reportLines{end+1} = '    Global   : single test per Metric (family size 1, reported as-is, never treated as pre-corrected).';
    reportLines{end+1} = '';
    reportLines{end+1} = sprintf('Omnibus tests fit: %d band-dependent rows, %d TotalPower rows.', height(OmnibusBandT), height(OmnibusTotalPowerT));
    reportLines{end+1} = sprintf('Band-dependent SignificantFDR: %d / %d', sum(OmnibusBandT.SignificantFDR), height(OmnibusBandT));
    reportLines{end+1} = sprintf('TotalPower SignificantFDR: %d / %d', sum(OmnibusTotalPowerT.SignificantFDR), height(OmnibusTotalPowerT));
    reportLines{end+1} = '';
    reportLines{end+1} = sprintf('Pairwise contrasts (N-E, R-E, R-N; Holm-corrected, family=3; gated on FDR survival): %d band-metric rows, %d TotalPower rows.', ...
        height(PairwiseBandT), height(PairwiseTotalPowerT));
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Judgment calls (disclosed) ---';
    for i = 1:numel(PARAMS11.JudgmentCalls)
        reportLines{end+1} = ['  - ' PARAMS11.JudgmentCalls{i}]; %#ok<AGROW>
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
    fid = fopen(fullfile(REPORT_DIR, sprintf('phase_median_robust_family_report_%s_%s.txt', TASK_LABEL, RUN_MODE_LABEL)), 'w');
    fprintf(fid, '%s\n', reportText);
    fclose(fid);

end % TASK_TYPES loop

fprintf('\n=== X11 complete for all three Tasks. This is a purely additive exploratory family -- X05/X06/X08/X09/X10 outputs are untouched. ===\n');

%% ======================= LOCAL FUNCTIONS (all after top-level code) =======================
%  local_subject_numeric, local_fdr_bh, local_holm_bonferroni,
%  local_condition_contrast, local_build_table, local_yn, local_fmt are
%  copied VERBATIM from X05_primary_cognitive_statistics.m (self-
%  contained-file convention, same as every other script in this
%  pipeline). local_fit_condition_omnibus, local_bh_fdr_or_raw, and
%  local_run_pairwise are NEW helpers written for this script, factoring
%  out logic that X05/X08 inline once but this script must repeat across
%  6 separate analysis contexts (Electrode/ROI/Global x band-dependent/
%  broadband).

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
            error('Phase11:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name1);
        end
        cRow(idx1) = cRow(idx1) - 1;
    end
    if ~strcmp(level2, refLevel)
        name2 = sprintf('%s_%s', varName, level2);
        idx2 = find(strcmp(coeffNames, name2));
        if isempty(idx2)
            error('Phase11:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name2);
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
        error('Phase11:BuildTableColumnMismatch', ...
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
            error('Phase11:BuildTableColumnMismatch', 'Column %d not classified as text or numeric.', j);
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

function fr = local_fit_condition_omnibus(subT, metricName, conditionOrder)
% Fits Metric ~ Condition + (1|Subject), REML, and extracts the omnibus
% Satterthwaite F-test for Condition. Returns a struct rather than
% erroring out, so the many calling loops in this script can log a
% failure (via fr.errMsg) and continue, exactly like X05/X08's own
% try/catch pattern around fitlme.
    fr = struct('ok', false, 'lme', [], 'F', NaN, 'DF1', NaN, 'DF2', NaN, 'p', NaN, 'warnMsg', '', 'errMsg', '');
    if isempty(subT) || ~ismember(metricName, subT.Properties.VariableNames)
        fr.errMsg = 'Empty subT or missing metric column.';
        return;
    end
    subT.Subject = categorical(local_subject_numeric(subT.Subject));
    subT.Condition = categorical(cellstr(subT.Condition), conditionOrder);
    formula = sprintf('%s ~ Condition + (1|Subject)', metricName);
    try
        lastwarn('');
        lme = fitlme(subT, formula, 'FitMethod', 'REML');
        [warnMsg, ~] = lastwarn();
        fr.warnMsg = warnMsg;
    catch ME
        fr.errMsg = ME.message;
        return;
    end
    aov = anova(lme, 'DFMethod', 'satterthwaite');
    rowIdx = find(strcmp(aov.Term, 'Condition'), 1);
    if isempty(rowIdx)
        fr.errMsg = 'anova(lme) has no Condition row.';
        return;
    end
    fr.ok = true;
    fr.lme = lme;
    fr.F = aov.FStat(rowIdx);
    fr.DF1 = aov.DF1(rowIdx);
    fr.DF2 = aov.DF2(rowIdx);
    fr.p = aov.pValue(rowIdx);
end

function [adjP, reject] = local_bh_fdr_or_raw(pvals, alphaFDR)
% Family size 1 is reported as-is, never treated as pre-corrected --
% same convention X05 already uses for its Primary tier (family size 1).
    if numel(pvals) == 1
        adjP = pvals;
        reject = pvals <= alphaFDR;
    else
        [adjP, reject] = local_fdr_bh(pvals, alphaFDR);
    end
end

function rows3 = local_run_pairwise(lme, refLevel, CONTRASTS, alphaFDR)
% Runs the 3 pairwise contrasts defined by CONTRASTS (rows: {lvl2, lvl1,
% label}, label = lvl2 - lvl1), Holm-corrects their p-values (family
% size 3), and returns a 3x9 cell array of
% {Contrast, Estimate, SE, DF, tStat, pValue, pHolm, CI_Low, CI_High} --
% identical arithmetic to X05/X08's own inline pairwise-contrast block,
% factored into a function since this script repeats it across 6
% separate analysis contexts (Electrode/ROI/Global x band-dependent/
% broadband).
    coeffNames = lme.CoefficientNames;
    beta = lme.Coefficients.Estimate;
    covB = lme.CoefficientCovariance;
    contrastP = nan(3, 1);
    contrastData = cell(3, 8);
    for c = 1:3
        lvl2 = CONTRASTS{c, 1}; lvl1 = CONTRASTS{c, 2}; label = CONTRASTS{c, 3};
        C = local_condition_contrast(coeffNames, 'Condition', refLevel, lvl1, lvl2);
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
