%% ========================================================================
% X12_band_ratio_indices_cognitive.m
%
% Purpose: exploratory EEG band-ratio indices (exploratory extension
% family 2 of 5: X11-X15). Adds derived band-ratio metrics alongside the
% two locked co-primary metrics (LogAbsolutePower, RelativePower).
% Never modifies, reads, or depends on the statistical outputs of
% X05/X06/X08/X09/X10; reads only X04's subject-level aggregation
% (data), X07's shared ROI/electrode definitions, and sibling script
% X11's median-robust data tables (data). The four ratio indices are
% tested in their own separate FDR families -- never pooled with
% LogAbsolutePower/RelativePower, with X11's own metrics, or with each
% other.
%
% Rationale: Theta/Beta ratio (TBR) is a widely used EEG index of
% attention/vigilance, with elevated TBR linked to reduced sustained
% attention and increased mind-wandering (e.g. Putman et al., 2010;
% Angelidis et al., 2016). Alpha/Theta and Alpha/Beta are companion
% indices capturing relaxation/drowsiness and cortical arousal
% respectively. The Engagement Index (Beta / (Alpha + Theta)) originates
% from Pope, Bogart & Bartolome (1995) and is a standard index of
% cognitive workload/engagement. All four are descriptive/exploratory
% only; none replace or pool with the locked co-primary metrics.
%
% The X11 contract (sibling script, trusted as given, not re-derived):
% X11_median_robust_family_cognitive.m saves, per Task, median_robust_
% data.mat containing MedianRobustElectrodeT/MedianRobustROIT/
% MedianRobustGlobalT (Subject, TaskType, BlockNumber, Condition,
% Channel/ROI, Band, AbsMedian, RelMedian, CV_Absolute, CV_Relative).
% This script uses only the AbsMedian column (the median-based analogue
% of AbsolutePower); each table is loaded read-only and defensively
% re-filtered to this script's own Task and to Condition in {E,N,R}
% (mismatches/exclusions are counted and logged).
%
% Method: four band-ratio indices per Subject x Condition (E/N/R;
% Baseline excluded): ThetaBeta = Power(Theta)/Power(Beta); AlphaTheta =
% Power(Alpha)/Power(Theta); AlphaBeta = Power(Alpha)/Power(Beta);
% EngagementIndex = Power(Beta)/(Power(Alpha)+Power(Theta)). Each index
% has two variants: MeanBased (from X04's AbsolutePower, with ROI-level
% ratios computed as ratio-of-ROI-mean-powers, not mean-of-per-channel-
% ratios) and MedianBased (from X11's AbsMedian at the matching level).
% Computed at three spatial levels: Electrode (62 electrodes, P9/P10
% included, no ROI averaging), ROI (5 ROIs, X07's definition; ROI-level
% AbsolutePower is recomputed fresh here as the arithmetic mean of
% channel-level power across each ROI's member channels, for
% Theta/Alpha/Beta only, never read from X08's own statistical output),
% and Global (whole-scalp). LMM: RatioValue ~ Condition (reference = E) +
% (1|Subject), REML, Satterthwaite DF, fit per Task independently.
%
% FDR family structure: for each of the 3 spatial levels, for each of the
% 4 RatioType x 2 Variant = 8 combinations, BH-FDR is applied within that
% family across the level's units (Electrode: family size 62; ROI: family
% size 5; Global: family size 1, reported as-is) -- 24 independent FDR
% families per Task, never pooled across RatioType, Variant, level, or
% with any other script's families. Pairwise contrasts are N-E, R-E, R-N
% (matching X05/X06's convention), computed only for units whose omnibus
% test survives its own FDR family, Holm-corrected within that unit
% (family size 3).
%
% Requirements: base MATLAB + Statistics and Machine Learning Toolbox
% (fitlme, anova on LinearMixedModel, coefTest).
%
% Inputs (read-only, never modified):
%   <ANALYSIS_ROOT>/04_subject_level/phase3_subject_level_bandpower.mat
%     (BandPowerT, GlobalBandPowerT, PARAMS3) -- X04.
%   <ANALYSIS_ROOT>/07_roi_electrode_definitions/
%     phase5_channel_roi_definitions.mat (ChannelDefT, ROI_DEFINITION,
%     ROI_NAMES, UNASSIGNED_FOR_ROI, PARAMS5D) -- X07, shared across Tasks.
%   <ANALYSIS_ROOT>/<Full|Nback|SART>/11_median_robust_family/<RUN_MODE>/
%     median_robust_data.mat (MedianRobustElectrodeT, MedianRobustROIT,
%     MedianRobustGlobalT) -- X11 (sibling script; see contract above).
%
% Outputs (one complete set per Task; nothing under X05/X06/X08/X09/X10/
% X11's own output directories is ever written to):
%   <ANALYSIS_ROOT>/<Full|Nback|SART>/12_band_ratio_indices/
%       phase12_parameters.mat
%       phase12_band_ratio_electrode_values.mat / .csv / .xlsx
%       phase12_band_ratio_electrode_omnibus.mat / .csv / .xlsx
%       phase12_band_ratio_electrode_pairwise.mat / .csv / .xlsx
%       phase12_band_ratio_roi_values.mat / .csv / .xlsx
%       phase12_band_ratio_roi_omnibus.mat / .csv / .xlsx
%       phase12_band_ratio_roi_pairwise.mat / .csv / .xlsx
%       phase12_band_ratio_global_values.mat / .csv / .xlsx
%       phase12_band_ratio_global_omnibus.mat / .csv / .xlsx
%       phase12_band_ratio_global_pairwise.mat / .csv / .xlsx
%   <ANALYSIS_ROOT>/09_reports/
%       phase_band_ratio_indices_report_<Full|Nback|SART>.txt
%
% RUN_MODE: full_sample / sensitivity_excl_noisy (see
% EXCLUSION_CONFIG_cognitive.m); run this script once per mode. This
% script's own X04-side data (Electrode/Global) is exclusion-filtered
% directly; its X11 sibling-contract data is already exclusion-filtered
% upstream, since this script reads X11's own RUN_MODE-labeled subfolder.
% ========================================================================

clear; clc;

%% Section 0: Paths and existence checks
config;
DATA_ROOT_FOR_SIBLING = DATA_DIR;
[DATA_ROOT_PARENT, ~] = fileparts(DATA_ROOT_FOR_SIBLING);
ANALYSIS_ROOT = OUT_DIR;

SUBJ_DIR   = fullfile(ANALYSIS_ROOT, '04_subject_level');
DEF_DIR    = fullfile(ANALYSIS_ROOT, '07_roi_electrode_definitions');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');

if ~exist(SUBJ_DIR, 'dir')
    error('Phase12:MissingInput', 'X04 output directory not found: %s', SUBJ_DIR);
end
bandPowerFile = fullfile(SUBJ_DIR, 'phase3_subject_level_bandpower.mat');
if ~exist(bandPowerFile, 'file')
    error('Phase12:MissingInput', 'Required X04 file not found: %s', bandPowerFile);
end
roiDefFile = fullfile(DEF_DIR, 'phase5_channel_roi_definitions.mat');
if ~exist(roiDefFile, 'file')
    error('Phase12:MissingInput', ['X07 ROI/electrode definition file not found: %s\n' ...
        'Run X07_roi_electrode_definitions_cognitive.m first. This script will NOT ' ...
        'reconstruct or guess ROI/electrode membership on its own.'], roiDefFile);
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
    error('Phase12:BadRunMode', 'RUN_MODE must be one of: %s', strjoin(RUN_MODES_cognitive, ', '));
end
RUN_MODE_LABEL = RUN_MODES_DISPLAY_cognitive.(RUN_MODE);
fprintf('RUN_MODE = %s\n', RUN_MODE_LABEL);
if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
    fprintf('Excluding subject(s) {%s} from every model in this run -- see EXCLUSION_CONFIG_cognitive.m for rationale.\n', ...
        strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','));
end
% ---- END NEW ----

if isempty(ver('stats'))
    error('Phase12:MissingToolbox', ...
        'Statistics and Machine Learning Toolbox not detected. This script requires fitlme, anova(LinearMixedModel), and coefTest.');
end

%% Section 1: Load X04 output (shared across all three Task analyses)
S4 = load(bandPowerFile, 'BandPowerT', 'GlobalBandPowerT', 'PARAMS3');
BandPowerT_ALL = S4.BandPowerT;
GlobalBandPowerT_ALL = S4.GlobalBandPowerT;
PARAMS3 = S4.PARAMS3;
clear S4;

requiredElecCols = {'Subject','TaskType','BlockNumber','Condition','Channel','Band','AbsolutePower'};
missingElecCols = requiredElecCols(~ismember(requiredElecCols, BandPowerT_ALL.Properties.VariableNames));
if ~isempty(missingElecCols)
    error('Phase12:MissingColumns', 'BandPowerT is missing expected column(s): %s', strjoin(missingElecCols, ', '));
end
requiredGlobalCols = {'Subject','TaskType','BlockNumber','Condition','Band','AbsolutePower'};
missingGlobalCols = requiredGlobalCols(~ismember(requiredGlobalCols, GlobalBandPowerT_ALL.Properties.VariableNames));
if ~isempty(missingGlobalCols)
    error('Phase12:MissingColumns', 'GlobalBandPowerT is missing expected column(s): %s', strjoin(missingGlobalCols, ', '));
end

%% Section 2: Load X07 output (shared, Task-invariant ROI/electrode definitions)
S7 = load(roiDefFile, 'ChannelDefT', 'ROI_DEFINITION', 'ROI_NAMES', 'UNASSIGNED_FOR_ROI', 'PARAMS5D');
ChannelDefT = S7.ChannelDefT;
ROI_DEFINITION = S7.ROI_DEFINITION;
ROI_NAMES = S7.ROI_NAMES;
UNASSIGNED_FOR_ROI = S7.UNASSIGNED_FOR_ROI;
PARAMS5D = S7.PARAMS5D; %#ok<NASGU>
clear S7;

ELECTRODES = ChannelDefT.Label;   % locked 62-electrode reference order (P9/P10 included)
nElec = numel(ELECTRODES);
nROI = numel(ROI_NAMES);

%% Section 3: Parameters shared across the whole script
RATIO_TYPES = {'ThetaBeta', 'AlphaTheta', 'AlphaBeta', 'EngagementIndex'};
VARIANTS    = {'MeanBased', 'MedianBased'};
CONDITION_ORDER = {'E', 'N', 'R'};   % Baseline excluded from this analysis
REF_LEVEL = CONDITION_ORDER{1};
FDR_ALPHA = 0.05;
RATIO_BANDS_NEEDED = {'Theta', 'Alpha', 'Beta'};

allBands = fieldnames(PARAMS3.BandDefinitionsHz);
missingRatioBands = RATIO_BANDS_NEEDED(~ismember(RATIO_BANDS_NEEDED, allBands));
if ~isempty(missingRatioBands)
    error('Phase12:BandMismatch', 'PARAMS3.BandDefinitionsHz is missing required band(s): %s', strjoin(missingRatioBands, ', '));
end

TASK_TYPES = {'full', 'nback', 'sart'};
TASK_DISPLAY = struct('full', 'Full', 'nback', 'Nback', 'sart', 'SART');

fprintf('========================================================\n');
fprintf('PHASE X12 -- EXPLORATORY BAND-RATIO INDICES, COGNITIVE-TASK DATA\n');
fprintf('Running as THREE INDEPENDENT, PARALLEL analyses (Full / Nback / SART).\n');
fprintf('========================================================\n\n');

%% ========================================================================
%  TASK LOOP -- everything below is repeated independently per Task
% ========================================================================
for taskI = 1:numel(TASK_TYPES)
    TASK_KEY = TASK_TYPES{taskI};
    TASK_LABEL = TASK_DISPLAY.(TASK_KEY);

    STATS_DIR = fullfile(ANALYSIS_ROOT, TASK_LABEL, '12_band_ratio_indices', RUN_MODE);
    if ~exist(STATS_DIR, 'dir'); mkdir(STATS_DIR); end

    fprintf('\n########################################################\n');
    fprintf('### TASK = %s | RUN_MODE = %s\n', TASK_LABEL, RUN_MODE_LABEL);
    fprintf('########################################################\n\n');

    integrityIssues = {};

    %% Section 4: filter X04 electrode-level and global-level data to this Task, then to E/N/R
    thisTaskElecT = BandPowerT_ALL(strcmp(BandPowerT_ALL.TaskType, TASK_KEY), :);
    nBaselineExcludedElec = sum(strcmp(thisTaskElecT.Condition, 'Baseline'));
    BandPowerT = thisTaskElecT(ismember(thisTaskElecT.Condition, CONDITION_ORDER), :);
    fprintf('Task=%s: %d electrode-level rows loaded, %d Baseline rows excluded, %d E/N/R rows retained.\n', ...
        TASK_LABEL, height(thisTaskElecT), nBaselineExcludedElec, height(BandPowerT));

    thisTaskGlobalT = GlobalBandPowerT_ALL(strcmp(GlobalBandPowerT_ALL.TaskType, TASK_KEY), :);
    nBaselineExcludedGlobal = sum(strcmp(thisTaskGlobalT.Condition, 'Baseline'));
    GlobalBandPowerT = thisTaskGlobalT(ismember(thisTaskGlobalT.Condition, CONDITION_ORDER), :);
    fprintf('Task=%s: %d global-level rows loaded, %d Baseline rows excluded, %d E/N/R rows retained.\n', ...
        TASK_LABEL, height(thisTaskGlobalT), nBaselineExcludedGlobal, height(GlobalBandPowerT));

    % ---- NEW: apply the shared noisy-subject exclusion list (see
    % EXCLUSION_CONFIG_cognitive.m), only when RUN_MODE=sensitivity_excl_noisy.
    % Never applied to the full_sample run. Applied here to the X04 sources
    % (electrode-level BandPowerT, global-level GlobalBandPowerT); X11's
    % sibling-contract tables (Section 6 below) are already exclusion-
    % filtered upstream, since this script now reads X11's own
    % RUN_MODE-labeled output subfolder. ----
    nExcludedRowsSensitivityElec = 0;
    nExcludedRowsSensitivityGlobal = 0;
    if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
        subjNumForExclElec = local_subject_numeric(BandPowerT.Subject);
        excludeMaskElec = ismember(subjNumForExclElec, EXCLUDE_NOISY_SUBJECTS);
        nExcludedRowsSensitivityElec = sum(excludeMaskElec);
        BandPowerT = BandPowerT(~excludeMaskElec, :);

        subjNumForExclGlobal = local_subject_numeric(GlobalBandPowerT.Subject);
        excludeMaskGlobal = ismember(subjNumForExclGlobal, EXCLUDE_NOISY_SUBJECTS);
        nExcludedRowsSensitivityGlobal = sum(excludeMaskGlobal);
        GlobalBandPowerT = GlobalBandPowerT(~excludeMaskGlobal, :);

        fprintf('Task=%s: RUN_MODE=%s -- excluded %d electrode-level and %d global-level row(s) belonging to subject(s) {%s}; %d/%d rows remain.\n', ...
            TASK_LABEL, RUN_MODE_LABEL, nExcludedRowsSensitivityElec, nExcludedRowsSensitivityGlobal, ...
            strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','), height(BandPowerT), height(GlobalBandPowerT));
    end

    if isempty(BandPowerT) || isempty(GlobalBandPowerT)
        warning('Task=%s: no E/N/R rows found after filtering -- skipping this Task entirely.', TASK_LABEL);
        continue;
    end

    %% Section 5: data-driven design verification (mirrors X05/X08 Section 2)
    comboT = unique(BandPowerT(:, {'Subject','BlockNumber','Condition'}), 'rows');
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
    fprintf('Task=%s: %d recordings, %d subjects, condition-pair coverage: %s\n', TASK_LABEL, nRecordings, nSubjects, ...
        strjoin(cellfun(@(k) sprintf('%s=%d', strrep(k,'_',' & '), pairTypeCount(k)), keys(pairTypeCount), 'UniformOutput', false), ', '));

    allSubjectsComplete3 = false;
    condsPerSubjCounts = cellfun(@numel, subjPairs);
    if all(condsPerSubjCounts >= 3)
        allSubjectsComplete3 = true;
        integrityIssues{end+1} = sprintf(['Task=%s: every subject appears to have 3+ recordings/conditions -- the incomplete-design ' ...
            'justification for using an LMM instead of RM-ANOVA no longer applies. STOP and re-examine.'], TASK_LABEL);
    end
    if allSubjectsComplete3
        error('Phase12:DesignAssumptionViolated', ...
            'Task=%s: design verification found a complete 3-condition-per-subject structure. Stopping before fitting -- see integrityIssues.', TASK_LABEL);
    end

    %% Section 6: load X11 median-robust tables for this Task (sibling script contract)
    % ---- X11 saves under a RUN_MODE-labeled subfolder too, so this path
    % must match it exactly --
    % this script must be run with the SAME RUN_MODE value as the X11 run
    % it is reading, so the two stay consistent (X11 excludes noisy
    % subjects at its own source, this script never re-excludes them). ----
    MEDIAN_DIR = fullfile(ANALYSIS_ROOT, TASK_LABEL, '11_median_robust_family', RUN_MODE);
    medianFile = fullfile(MEDIAN_DIR, 'median_robust_data.mat');
    if ~exist(medianFile, 'file')
        error('Phase12:MissingInput', ['X11 median-robust file not found: %s\n' ...
            'Run X11_median_robust_family_cognitive.m first (with the SAME RUN_MODE=%s). This ' ...
            'script trusts X11''s exact contract (see header) and will not reconstruct ' ...
            'median-based band power itself.'], medianFile, RUN_MODE);
    end
    S11 = load(medianFile, 'MedianRobustElectrodeT', 'MedianRobustROIT', 'MedianRobustGlobalT');

    requiredMedElecCols = {'Subject','TaskType','BlockNumber','Condition','Channel','Band','AbsMedian'};
    missingMedElecCols = requiredMedElecCols(~ismember(requiredMedElecCols, S11.MedianRobustElectrodeT.Properties.VariableNames));
    if ~isempty(missingMedElecCols)
        error('Phase12:MissingColumns', 'MedianRobustElectrodeT (X11) is missing expected column(s): %s', strjoin(missingMedElecCols, ', '));
    end
    requiredMedRoiCols = {'Subject','TaskType','BlockNumber','Condition','ROI','Band','AbsMedian'};
    missingMedRoiCols = requiredMedRoiCols(~ismember(requiredMedRoiCols, S11.MedianRobustROIT.Properties.VariableNames));
    if ~isempty(missingMedRoiCols)
        error('Phase12:MissingColumns', 'MedianRobustROIT (X11) is missing expected column(s): %s', strjoin(missingMedRoiCols, ', '));
    end
    requiredMedGlobalCols = {'Subject','TaskType','BlockNumber','Condition','Band','AbsMedian'};
    missingMedGlobalCols = requiredMedGlobalCols(~ismember(requiredMedGlobalCols, S11.MedianRobustGlobalT.Properties.VariableNames));
    if ~isempty(missingMedGlobalCols)
        error('Phase12:MissingColumns', 'MedianRobustGlobalT (X11) is missing expected column(s): %s', strjoin(missingMedGlobalCols, ', '));
    end

    [MedianRobustElectrodeT, nElecTaskMismatch, nElecBaselineExcl] = local_filter_task_condition(S11.MedianRobustElectrodeT, TASK_KEY, CONDITION_ORDER);
    [MedianRobustROIT, nRoiTaskMismatch, nRoiBaselineExcl] = local_filter_task_condition(S11.MedianRobustROIT, TASK_KEY, CONDITION_ORDER);
    [MedianRobustGlobalT, nGlobTaskMismatch, nGlobBaselineExcl] = local_filter_task_condition(S11.MedianRobustGlobalT, TASK_KEY, CONDITION_ORDER);
    clear S11;

    if nElecTaskMismatch > 0
        integrityIssues{end+1} = sprintf('X11 MedianRobustElectrodeT (Task=%s): %d row(s) had a TaskType other than "%s" and were dropped before use.', TASK_LABEL, nElecTaskMismatch, TASK_KEY); %#ok<AGROW>
    end
    if nRoiTaskMismatch > 0
        integrityIssues{end+1} = sprintf('X11 MedianRobustROIT (Task=%s): %d row(s) had a TaskType other than "%s" and were dropped before use.', TASK_LABEL, nRoiTaskMismatch, TASK_KEY); %#ok<AGROW>
    end
    if nGlobTaskMismatch > 0
        integrityIssues{end+1} = sprintf('X11 MedianRobustGlobalT (Task=%s): %d row(s) had a TaskType other than "%s" and were dropped before use.', TASK_LABEL, nGlobTaskMismatch, TASK_KEY); %#ok<AGROW>
    end
    fprintf('Task=%s: X11 tables loaded -- Electrode %d Baseline excluded, ROI %d Baseline excluded, Global %d Baseline excluded (of any found).\n', ...
        TASK_LABEL, nElecBaselineExcl, nRoiBaselineExcl, nGlobBaselineExcl);

    %% Section 7: ROI aggregation (MEAN-BASED variant only), Theta/Alpha/Beta,
    %  recomputed independently from X04+X07 -- X08's own ROI statistics are
    %  never read here.
    BandPowerT.ChannelUC = upper(BandPowerT.Channel);
    recCombos = unique(BandPowerT(:, {'Subject','Condition'}), 'rows');
    roiRows = cell(0, 5);   % Subject, Condition, ROI, Band, AbsolutePower
    roiChanMismatch = {};
    nonPositiveROIAbs = {};

    for i = 1:height(recCombos)
        subjStr = recCombos.Subject{i};
        condStr = recCombos.Condition{i};
        recMask = strcmp(BandPowerT.Subject, subjStr) & strcmp(BandPowerT.Condition, condStr);
        for r = 1:nROI
            roiName = ROI_NAMES{r};
            roiChansUC = upper(ROI_DEFINITION.(roiName));
            chanMask = ismember(BandPowerT.ChannelUC, roiChansUC);
            for b = 1:numel(RATIO_BANDS_NEEDED)
                bandName = RATIO_BANDS_NEEDED{b};
                rowMask = recMask & chanMask & strcmp(BandPowerT.Band, bandName);
                nChanFound = sum(rowMask);
                if nChanFound ~= numel(roiChansUC)
                    roiChanMismatch{end+1} = sprintf('Sub%s Cond%s ROI=%s Band=%s: found %d/%d expected ROI channels.', ...
                        subjStr, condStr, roiName, bandName, nChanFound, numel(roiChansUC)); %#ok<AGROW>
                end
                roiAbs = mean(BandPowerT.AbsolutePower(rowMask));
                if ~(roiAbs > 0)
                    nonPositiveROIAbs{end+1} = sprintf('Sub%s Cond%s ROI=%s Band=%s', subjStr, condStr, roiName, bandName); %#ok<AGROW>
                end
                roiRows(end+1, :) = {subjStr, condStr, roiName, bandName, roiAbs}; %#ok<AGROW>
            end
        end
    end
    if ~isempty(roiChanMismatch)
        integrityIssues{end+1} = sprintf('Task=%s: %d ROI x Band x recording combination(s) did not find the expected number of ROI member channels: %s', ...
            TASK_LABEL, numel(roiChanMismatch), strjoin(roiChanMismatch, '; '));
    end
    if ~isempty(nonPositiveROIAbs)
        integrityIssues{end+1} = sprintf('Task=%s: %d ROI x Band x recording combination(s) had non-positive mean ROI AbsolutePower (ratio for these will be non-finite, flagged not patched): %s', ...
            TASK_LABEL, numel(nonPositiveROIAbs), strjoin(nonPositiveROIAbs, '; '));
    end
    RoiBandPowerT_mean = local_build_table(roiRows, [1 2 3 4], 5, {'Subject','Condition','ROI','Band','AbsolutePower'});
    fprintf('Task=%s: ROI aggregation (mean-based, Theta/Alpha/Beta only) complete: %d rows (%d recordings x %d ROIs x %d bands).\n', ...
        TASK_LABEL, height(RoiBandPowerT_mean), height(recCombos), nROI, numel(RATIO_BANDS_NEEDED));

    %% Section 8: compute the four ratios, MeanBased variant, all 3 levels
    [ElecRatio_mean, nSkipElecMean, nNonFiniteElecMean] = local_compute_band_ratios(BandPowerT, {'Subject','Condition','Channel'}, 'Band', 'AbsolutePower');
    [RoiRatio_mean, nSkipRoiMean, nNonFiniteRoiMean] = local_compute_band_ratios(RoiBandPowerT_mean, {'Subject','Condition','ROI'}, 'Band', 'AbsolutePower');
    [GlobalRatio_mean, nSkipGlobMean, nNonFiniteGlobMean] = local_compute_band_ratios(GlobalBandPowerT, {'Subject','Condition'}, 'Band', 'AbsolutePower');

    %% Section 9: compute the four ratios, MedianBased variant, all 3 levels
    [ElecRatio_median, nSkipElecMed, nNonFiniteElecMed] = local_compute_band_ratios(MedianRobustElectrodeT, {'Subject','Condition','Channel'}, 'Band', 'AbsMedian');
    [RoiRatio_median, nSkipRoiMed, nNonFiniteRoiMed] = local_compute_band_ratios(MedianRobustROIT, {'Subject','Condition','ROI'}, 'Band', 'AbsMedian');
    [GlobalRatio_median, nSkipGlobMed, nNonFiniteGlobMed] = local_compute_band_ratios(MedianRobustGlobalT, {'Subject','Condition'}, 'Band', 'AbsMedian');

    nSkippedTotal = nSkipElecMean + nSkipRoiMean + nSkipGlobMean + nSkipElecMed + nSkipRoiMed + nSkipGlobMed;
    nNonFiniteTotal = nNonFiniteElecMean + nNonFiniteRoiMean + nNonFiniteGlobMean + nNonFiniteElecMed + nNonFiniteRoiMed + nNonFiniteGlobMed;
    if nSkippedTotal > 0
        integrityIssues{end+1} = sprintf(['Task=%s: %d Subject x Condition x Unit combination(s) across all levels/variants were skipped ' ...
            'from ratio computation because Theta/Alpha/Beta were not all available (Elec/ROI/Global, Mean=%d/%d/%d, Median=%d/%d/%d).'], ...
            TASK_LABEL, nSkippedTotal, nSkipElecMean, nSkipRoiMean, nSkipGlobMean, nSkipElecMed, nSkipRoiMed, nSkipGlobMed);
    end
    if nNonFiniteTotal > 0
        integrityIssues{end+1} = sprintf(['Task=%s: %d computed ratio value(s) across all levels/variants were non-finite ' ...
            '(non-positive denominator power; flagged, never silently patched).'], TASK_LABEL, nNonFiniteTotal);
    end

    ElecRatio_mean.Variant   = repmat({'MeanBased'},   height(ElecRatio_mean),   1);
    ElecRatio_median.Variant = repmat({'MedianBased'}, height(ElecRatio_median), 1);
    ElecRatioT = [ElecRatio_mean; ElecRatio_median];

    RoiRatio_mean.Variant   = repmat({'MeanBased'},   height(RoiRatio_mean),   1);
    RoiRatio_median.Variant = repmat({'MedianBased'}, height(RoiRatio_median), 1);
    RoiRatioT = [RoiRatio_mean; RoiRatio_median];

    GlobalRatio_mean.Variant   = repmat({'MeanBased'},   height(GlobalRatio_mean),   1);
    GlobalRatio_median.Variant = repmat({'MedianBased'}, height(GlobalRatio_median), 1);
    GlobalRatioT = [GlobalRatio_mean; GlobalRatio_median];

    fprintf('Task=%s: ratio value tables built -- Electrode %d rows, ROI %d rows, Global %d rows.\n', ...
        TASK_LABEL, height(ElecRatioT), height(RoiRatioT), height(GlobalRatioT));

    %% Section 10: model fitting + FDR + pairwise contrasts, per level
    [ElecOmnibusT, ElecPairwiseT] = local_fit_ratio_family(ElecRatioT, 'Channel', ELECTRODES, RATIO_TYPES, VARIANTS, CONDITION_ORDER, REF_LEVEL, FDR_ALPHA);
    [RoiOmnibusT,  RoiPairwiseT]  = local_fit_ratio_family(RoiRatioT,  'ROI',     ROI_NAMES,  RATIO_TYPES, VARIANTS, CONDITION_ORDER, REF_LEVEL, FDR_ALPHA);
    [GlobalOmnibusT, GlobalPairwiseT] = local_fit_ratio_family(GlobalRatioT, '', {'Global'}, RATIO_TYPES, VARIANTS, CONDITION_ORDER, REF_LEVEL, FDR_ALPHA);

    nExpectedElec = numel(RATIO_TYPES) * numel(VARIANTS) * nElec;
    nExpectedRoi  = numel(RATIO_TYPES) * numel(VARIANTS) * nROI;
    nExpectedGlob = numel(RATIO_TYPES) * numel(VARIANTS) * 1;
    if height(ElecOmnibusT) ~= nExpectedElec
        integrityIssues{end+1} = sprintf('Task=%s: Electrode-level omnibus table has %d rows, expected %d (%d ratios x %d variants x %d electrodes) -- some fit(s) failed, see FamilySize per family.', ...
            TASK_LABEL, height(ElecOmnibusT), nExpectedElec, numel(RATIO_TYPES), numel(VARIANTS), nElec);
    end
    if height(RoiOmnibusT) ~= nExpectedRoi
        integrityIssues{end+1} = sprintf('Task=%s: ROI-level omnibus table has %d rows, expected %d (%d ratios x %d variants x %d ROIs) -- some fit(s) failed, see FamilySize per family.', ...
            TASK_LABEL, height(RoiOmnibusT), nExpectedRoi, numel(RATIO_TYPES), numel(VARIANTS), nROI);
    end
    if height(GlobalOmnibusT) ~= nExpectedGlob
        integrityIssues{end+1} = sprintf('Task=%s: Global-level omnibus table has %d rows, expected %d (%d ratios x %d variants x 1) -- some fit(s) failed.', ...
            TASK_LABEL, height(GlobalOmnibusT), nExpectedGlob, numel(RATIO_TYPES), numel(VARIANTS));
    end

    fprintf('\n--- Global-level band-ratio omnibus results (Task=%s) ---\n', TASK_LABEL);
    disp(GlobalOmnibusT);

    %% Section 11: PARAMS12 and save outputs
    PARAMS12 = struct();
    PARAMS12.ScriptName = 'X12_band_ratio_indices_cognitive.m';
    PARAMS12.Task = TASK_LABEL;
    PARAMS12.GeneratedOn = datestr(now); %#ok<TNOW1,DATST>
    PARAMS12.AnalysisLevel = sprintf('Exploratory band-ratio indices, Task=%s only, Condition in {E,N,R} only (Baseline excluded).', TASK_LABEL);
    PARAMS12.RatioDefinitions = struct('ThetaBeta', 'Theta/Beta', 'AlphaTheta', 'Alpha/Theta', ...
        'AlphaBeta', 'Alpha/Beta', 'EngagementIndex', 'Beta/(Alpha+Theta)');
    PARAMS12.Variants = struct('MeanBased', 'X04 BandPowerT/GlobalBandPowerT AbsolutePower (or this script''s own ROI mean-of-channel-AbsolutePower aggregation)', ...
        'MedianBased', 'X11 AbsMedian at the matching level, per the X11 contract (see script header)');
    PARAMS12.Model = 'LinearMixedEffects_RandomIntercept';
    PARAMS12.FixedEffectsPrimary = 'Condition (3-level: E, N, R; reference = E)';
    PARAMS12.RandomEffects = '(1|Subject)';
    PARAMS12.FitMethod = 'REML';
    PARAMS12.DFMethod = 'Satterthwaite';
    PARAMS12.FDRAlpha = FDR_ALPHA;
    PARAMS12.FDRFamilyScheme = 'Separate BH-FDR family per (RatioType x Variant), WITHIN each spatial level (Electrode/ROI/Global); family sizes 62/5/1 respectively; never pooled across RatioType, Variant, level, or with any locked-metric or X11 family.';
    PARAMS12.PairwiseContrastMethod = 'Holm-Bonferroni across the 3 pairwise contrasts (N-E, R-E, R-N), computed only for units whose own (RatioType x Variant x level) omnibus test survives FDR.';
    PARAMS12.BaselineExclusionNote = sprintf('Electrode-level: %d Baseline row(s) excluded. Global-level: %d Baseline row(s) excluded. X11 tables: Electrode %d / ROI %d / Global %d Baseline row(s) excluded (if any).', ...
        nBaselineExcludedElec, nBaselineExcludedGlobal, nElecBaselineExcl, nRoiBaselineExcl, nGlobBaselineExcl);
    PARAMS12.SourcePhase3File = bandPowerFile;
    PARAMS12.SourcePhase5DFile = roiDefFile;
    PARAMS12.SourcePhase11File = medianFile;
    PARAMS12.DesignVerification = struct('nRecordings', nRecordings, 'nSubjects', nSubjects, ...
        'PairTypeCounts', {pairTypeCount.keys}, 'IntegrityIssues', {integrityIssues});
    PARAMS12.RunDateTime = datestr(now); %#ok<TNOW1,DATST>
    PARAMS12.RunMode = RUN_MODE;
    if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
        PARAMS12.ExcludedSubjects = EXCLUDE_NOISY_SUBJECTS;
        PARAMS12.ExcludedSubjectsRationale = EXCLUSION_RATIONALE_cognitive;
        PARAMS12.RowsExcludedForSensitivity = struct('ElectrodeRows', nExcludedRowsSensitivityElec, 'GlobalRows', nExcludedRowsSensitivityGlobal);
    else
        PARAMS12.ExcludedSubjects = [];
    end

    save(fullfile(STATS_DIR, 'phase12_parameters.mat'), 'PARAMS12');

    save(fullfile(STATS_DIR, 'phase12_band_ratio_electrode_values.mat'), 'ElecRatioT');
    save(fullfile(STATS_DIR, 'phase12_band_ratio_electrode_omnibus.mat'), 'ElecOmnibusT');
    save(fullfile(STATS_DIR, 'phase12_band_ratio_electrode_pairwise.mat'), 'ElecPairwiseT');
    writetable(ElecRatioT, fullfile(STATS_DIR, 'phase12_band_ratio_electrode_values.csv'));
    writetable(ElecRatioT, fullfile(STATS_DIR, 'phase12_band_ratio_electrode_values.xlsx'));
    writetable(ElecOmnibusT, fullfile(STATS_DIR, 'phase12_band_ratio_electrode_omnibus.csv'));
    writetable(ElecOmnibusT, fullfile(STATS_DIR, 'phase12_band_ratio_electrode_omnibus.xlsx'));
    if ~isempty(ElecPairwiseT)
        writetable(ElecPairwiseT, fullfile(STATS_DIR, 'phase12_band_ratio_electrode_pairwise.csv'));
        writetable(ElecPairwiseT, fullfile(STATS_DIR, 'phase12_band_ratio_electrode_pairwise.xlsx'));
    end

    save(fullfile(STATS_DIR, 'phase12_band_ratio_roi_values.mat'), 'RoiRatioT');
    save(fullfile(STATS_DIR, 'phase12_band_ratio_roi_omnibus.mat'), 'RoiOmnibusT');
    save(fullfile(STATS_DIR, 'phase12_band_ratio_roi_pairwise.mat'), 'RoiPairwiseT');
    writetable(RoiRatioT, fullfile(STATS_DIR, 'phase12_band_ratio_roi_values.csv'));
    writetable(RoiRatioT, fullfile(STATS_DIR, 'phase12_band_ratio_roi_values.xlsx'));
    writetable(RoiOmnibusT, fullfile(STATS_DIR, 'phase12_band_ratio_roi_omnibus.csv'));
    writetable(RoiOmnibusT, fullfile(STATS_DIR, 'phase12_band_ratio_roi_omnibus.xlsx'));
    if ~isempty(RoiPairwiseT)
        writetable(RoiPairwiseT, fullfile(STATS_DIR, 'phase12_band_ratio_roi_pairwise.csv'));
        writetable(RoiPairwiseT, fullfile(STATS_DIR, 'phase12_band_ratio_roi_pairwise.xlsx'));
    end

    save(fullfile(STATS_DIR, 'phase12_band_ratio_global_values.mat'), 'GlobalRatioT');
    save(fullfile(STATS_DIR, 'phase12_band_ratio_global_omnibus.mat'), 'GlobalOmnibusT');
    save(fullfile(STATS_DIR, 'phase12_band_ratio_global_pairwise.mat'), 'GlobalPairwiseT');
    writetable(GlobalRatioT, fullfile(STATS_DIR, 'phase12_band_ratio_global_values.csv'));
    writetable(GlobalRatioT, fullfile(STATS_DIR, 'phase12_band_ratio_global_values.xlsx'));
    writetable(GlobalOmnibusT, fullfile(STATS_DIR, 'phase12_band_ratio_global_omnibus.csv'));
    writetable(GlobalOmnibusT, fullfile(STATS_DIR, 'phase12_band_ratio_global_omnibus.xlsx'));
    if ~isempty(GlobalPairwiseT)
        writetable(GlobalPairwiseT, fullfile(STATS_DIR, 'phase12_band_ratio_global_pairwise.csv'));
        writetable(GlobalPairwiseT, fullfile(STATS_DIR, 'phase12_band_ratio_global_pairwise.xlsx'));
    end

    fprintf('\nSaved output to: %s\n', STATS_DIR);

    %% Section 12: Report
    reportLines = {};
    reportLines{end+1} = repmat('=', 1, 60);
    reportLines{end+1} = sprintf('PHASE X12 -- EXPLORATORY BAND-RATIO INDICES REPORT (Task=%s)', TASK_LABEL);
    reportLines{end+1} = repmat('=', 1, 60);
    reportLines{end+1} = '';
    reportLines{end+1} = 'Purely additive exploratory metrics, computed';
    reportLines{end+1} = 'ALONGSIDE (not replacing or pooled with) the two locked co-primary metrics';
    reportLines{end+1} = '(LogAbsolutePower, RelativePower). Never reads X05/X06/X08/X09/X10 outputs.';
    reportLines{end+1} = '';
    reportLines{end+1} = sprintf('Baseline rows excluded -- Electrode-level: %d, Global-level: %d, X11 Electrode/ROI/Global: %d/%d/%d.', ...
        nBaselineExcludedElec, nBaselineExcludedGlobal, nElecBaselineExcl, nRoiBaselineExcl, nGlobBaselineExcl);
    reportLines{end+1} = sprintf('Recordings verified: %d | Subjects: %d', nRecordings, nSubjects);
    reportLines{end+1} = sprintf('Condition-pair coverage: %s', strjoin(cellfun(@(k) sprintf('%s=%d', strrep(k,'_',' & '), pairTypeCount(k)), keys(pairTypeCount), 'UniformOutput', false), ', '));
    reportLines{end+1} = '';
    reportLines{end+1} = 'Model: Linear mixed-effects, Condition (fixed, 3-level) + (1|Subject), REML, Satterthwaite DF.';
    reportLines{end+1} = 'Ratio definitions: ThetaBeta=Theta/Beta, AlphaTheta=Alpha/Theta, AlphaBeta=Alpha/Beta, EngagementIndex=Beta/(Alpha+Theta).';
    reportLines{end+1} = 'Variants: MeanBased (X04 AbsolutePower, own independent ROI aggregation for ROI level), MedianBased (X11 AbsMedian contract).';
    reportLines{end+1} = 'FDR: separate BH-FDR family per (RatioType x Variant) WITHIN each level; family size = 62 (Electrode), 5 (ROI), 1 (Global, trivial).';
    reportLines{end+1} = 'Pairwise: N-E, R-E, R-N (matches X05/X06 convention); Holm-corrected within unit (family size 3); gated on FDR.';
    reportLines{end+1} = '';

    reportLines{end+1} = '--- GLOBAL LEVEL: all omnibus results (family size 1, trivial FDR) ---';
    for r = 1:height(GlobalOmnibusT)
        reportLines{end+1} = sprintf('  [%s] %-16s F(%g,%.1f)=%.3f p=%.4f pFDR=%.4f family_n=%d %s', ...
            GlobalOmnibusT.Variant{r}, GlobalOmnibusT.RatioType{r}, GlobalOmnibusT.DF1(r), GlobalOmnibusT.DF2(r), ...
            GlobalOmnibusT.FStat(r), GlobalOmnibusT.pValue(r), GlobalOmnibusT.pFDR(r), GlobalOmnibusT.FamilySize(r), ...
            local_yn(GlobalOmnibusT.SignificantFDR(r))); %#ok<AGROW>
    end
    reportLines{end+1} = '--- GLOBAL LEVEL: pairwise contrasts ---';
    if isempty(GlobalPairwiseT)
        reportLines{end+1} = '  (none -- no RatioType/Variant survived FDR)';
    else
        for r = 1:height(GlobalPairwiseT)
            reportLines{end+1} = sprintf('  [%s] %-16s %-8s est=%.4f (95%% CI %.4f to %.4f) t=%.3f df=%.1f p=%.4f pHolm=%.4f', ...
                GlobalPairwiseT.Variant{r}, GlobalPairwiseT.RatioType{r}, GlobalPairwiseT.Contrast{r}, GlobalPairwiseT.Estimate(r), ...
                GlobalPairwiseT.CI_Low(r), GlobalPairwiseT.CI_High(r), GlobalPairwiseT.tStat(r), GlobalPairwiseT.DF(r), ...
                GlobalPairwiseT.pValue(r), GlobalPairwiseT.pHolm(r)); %#ok<AGROW>
        end
    end
    reportLines{end+1} = '';

    reportLines{end+1} = sprintf('--- ROI LEVEL: all omnibus results (%d ROIs x %d ratios x %d variants = %d expected tests) ---', nROI, numel(RATIO_TYPES), numel(VARIANTS), nExpectedRoi);
    for r = 1:height(RoiOmnibusT)
        reportLines{end+1} = sprintf('  [%s] %-16s %-10s F(%g,%.1f)=%.3f p=%.4f pFDR=%.4f family_n=%d %s', ...
            RoiOmnibusT.Variant{r}, RoiOmnibusT.RatioType{r}, RoiOmnibusT.ROI{r}, RoiOmnibusT.DF1(r), RoiOmnibusT.DF2(r), ...
            RoiOmnibusT.FStat(r), RoiOmnibusT.pValue(r), RoiOmnibusT.pFDR(r), RoiOmnibusT.FamilySize(r), ...
            local_yn(RoiOmnibusT.SignificantFDR(r))); %#ok<AGROW>
    end
    reportLines{end+1} = '--- ROI LEVEL: pairwise contrasts (FDR-surviving ROI x RatioType x Variant only) ---';
    if isempty(RoiPairwiseT)
        reportLines{end+1} = '  (none -- no ROI x RatioType x Variant survived FDR)';
    else
        for r = 1:height(RoiPairwiseT)
            reportLines{end+1} = sprintf('  [%s] %-16s %-10s %-8s est=%.4f (95%% CI %.4f to %.4f) t=%.3f df=%.1f p=%.4f pHolm=%.4f', ...
                RoiPairwiseT.Variant{r}, RoiPairwiseT.RatioType{r}, RoiPairwiseT.ROI{r}, RoiPairwiseT.Contrast{r}, RoiPairwiseT.Estimate(r), ...
                RoiPairwiseT.CI_Low(r), RoiPairwiseT.CI_High(r), RoiPairwiseT.tStat(r), RoiPairwiseT.DF(r), ...
                RoiPairwiseT.pValue(r), RoiPairwiseT.pHolm(r)); %#ok<AGROW>
        end
    end
    reportLines{end+1} = '';

    reportLines{end+1} = sprintf('--- ELECTRODE LEVEL: per (RatioType x Variant) family summary (%d families, %d electrodes each nominal) ---', numel(RATIO_TYPES)*numel(VARIANTS), nElec);
    for vi = 1:numel(VARIANTS)
        for ri = 1:numel(RATIO_TYPES)
            idx = strcmp(ElecOmnibusT.Variant, VARIANTS{vi}) & strcmp(ElecOmnibusT.RatioType, RATIO_TYPES{ri});
            nInFamily = sum(idx);
            nSig = sum(idx & ElecOmnibusT.SignificantFDR);
            reportLines{end+1} = sprintf('  [%s] %-16s family_size=%d significant_FDR=%d', VARIANTS{vi}, RATIO_TYPES{ri}, nInFamily, nSig); %#ok<AGROW>
        end
    end
    reportLines{end+1} = '--- ELECTRODE LEVEL: FDR-surviving omnibus results only (full table in phase12_band_ratio_electrode_omnibus.csv) ---';
    sigElecRows = find(ElecOmnibusT.SignificantFDR)';
    if isempty(sigElecRows)
        reportLines{end+1} = '  (none)';
    else
        for r = sigElecRows
            reportLines{end+1} = sprintf('  [%s] %-16s %-8s F(%g,%.1f)=%.3f p=%.4f pFDR=%.4f', ...
                ElecOmnibusT.Variant{r}, ElecOmnibusT.RatioType{r}, ElecOmnibusT.Channel{r}, ElecOmnibusT.DF1(r), ElecOmnibusT.DF2(r), ...
                ElecOmnibusT.FStat(r), ElecOmnibusT.pValue(r), ElecOmnibusT.pFDR(r)); %#ok<AGROW>
        end
    end
    reportLines{end+1} = '--- ELECTRODE LEVEL: all pairwise contrasts (full table also in phase12_band_ratio_electrode_pairwise.csv) ---';
    if isempty(ElecPairwiseT)
        reportLines{end+1} = '  (none -- no electrode x RatioType x Variant survived FDR)';
    else
        for r = 1:height(ElecPairwiseT)
            reportLines{end+1} = sprintf('  [%s] %-16s %-8s %-8s est=%.4f (95%% CI %.4f to %.4f) t=%.3f df=%.1f p=%.4f pHolm=%.4f', ...
                ElecPairwiseT.Variant{r}, ElecPairwiseT.RatioType{r}, ElecPairwiseT.Channel{r}, ElecPairwiseT.Contrast{r}, ElecPairwiseT.Estimate(r), ...
                ElecPairwiseT.CI_Low(r), ElecPairwiseT.CI_High(r), ElecPairwiseT.tStat(r), ElecPairwiseT.DF(r), ...
                ElecPairwiseT.pValue(r), ElecPairwiseT.pHolm(r)); %#ok<AGROW>
        end
    end
    reportLines{end+1} = '';

    if ~isempty(integrityIssues)
        reportLines{end+1} = '--- INTEGRITY ISSUES (reported, not auto-corrected) ---';
        for i = 1:numel(integrityIssues)
            reportLines{end+1} = ['  - ' integrityIssues{i}]; %#ok<AGROW>
        end
    else
        reportLines{end+1} = 'Integrity issues: none.';
    end
    reportLines{end+1} = '';
    reportLines{end+1} = 'This script performs NO interpretation beyond the statistics specified';
    reportLines{end+1} = 'above, and does not select electrodes/ROIs/ratio types for any later use.';
    reportLines{end+1} = '';
    reportLines{end+1} = 'Output files written to:';
    reportLines{end+1} = ['  ' STATS_DIR];
    reportLines{end+1} = repmat('=', 1, 60);

    reportText = strjoin(reportLines, '\n');
    fprintf('%s\n', reportText);
    fid = fopen(fullfile(REPORT_DIR, sprintf('phase_band_ratio_indices_report_%s_%s.txt', TASK_LABEL, RUN_MODE_LABEL)), 'w');
    fprintf(fid, '%s\n', reportText);
    fclose(fid);

end % TASK_TYPES loop

fprintf('\n=== X12 complete for all three Tasks. Do not pool these results with X05/X06/X08/X09/X10/X11. ===\n');

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
            error('Phase12:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name1);
        end
        cRow(idx1) = cRow(idx1) - 1;
    end
    if ~strcmp(level2, refLevel)
        name2 = sprintf('%s_%s', varName, level2);
        idx2 = find(strcmp(coeffNames, name2));
        if isempty(idx2)
            error('Phase12:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name2);
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
        error('Phase12:BuildTableColumnMismatch', ...
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
            error('Phase12:BuildTableColumnMismatch', 'Column %d not classified as text or numeric.', j);
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

function k = local_make_key(T, keyVars)
% Builds a per-row text key by concatenating the given (text-valued)
% columns of T with a separator, so rows sharing the same combination of
% keyVars can be matched without assuming table row order.
    parts = cell(1, numel(keyVars));
    for j = 1:numel(keyVars)
        col = T.(keyVars{j});
        if iscategorical(col)
            col = cellstr(col);
        elseif ischar(col)
            col = cellstr(col);
        elseif isnumeric(col)
            col = cellstr(num2str(col));
        end
        parts{j} = col;
    end
    k = parts{1};
    for j = 2:numel(keyVars)
        k = strcat(k, '||', parts{j});
    end
end

function [RatioLongT, nSkippedMissingBand, nNonFiniteRatios] = local_compute_band_ratios(sourceT, keyVars, bandCol, valueCol)
% Computes the four band-ratio indices (ThetaBeta, AlphaTheta, AlphaBeta,
% EngagementIndex) from a long table sourceT that has, at minimum, the
% columns in keyVars plus bandCol (band name) and valueCol (numeric power,
% either AbsolutePower [MeanBased] or AbsMedian [MedianBased]). One
% Theta/Alpha/Beta triplet is required per unique combination of keyVars;
% a combination missing any of the three is skipped and counted (never
% silently dropped without a trace). Ratio = Power(numerator band(s)) /
% Power(denominator band(s)), computed directly from valueCol -- no
% further transform (no log, no z-score) is applied here.
    keyAll = local_make_key(sourceT, keyVars);
    bandVals = cellstr(sourceT.(bandCol));
    valCol = sourceT.(valueCol);

    thetaMask = strcmp(bandVals, 'Theta');
    alphaMask = strcmp(bandVals, 'Alpha');
    betaMask  = strcmp(bandVals, 'Beta');

    thetaMap = containers.Map(keyAll(thetaMask), num2cell(valCol(thetaMask)));
    alphaMap = containers.Map(keyAll(alphaMask), num2cell(valCol(alphaMask)));
    betaMap  = containers.Map(keyAll(betaMask),  num2cell(valCol(betaMask)));

    baseT = unique(sourceT(:, keyVars), 'rows');
    baseKeys = local_make_key(baseT, keyVars);
    nBase = height(baseT);

    ratioNames = {'ThetaBeta', 'AlphaTheta', 'AlphaBeta', 'EngagementIndex'};
    nRatios = numel(ratioNames);

    rows = cell(nBase * nRatios, numel(keyVars) + 2);
    rowCount = 0;
    nSkippedMissingBand = 0;
    nNonFiniteRatios = 0;

    for i = 1:nBase
        k = baseKeys{i};
        if ~(isKey(thetaMap, k) && isKey(alphaMap, k) && isKey(betaMap, k))
            nSkippedMissingBand = nSkippedMissingBand + 1;
            continue;
        end
        th = thetaMap(k); al = alphaMap(k); be = betaMap(k);
        ratioVals = [th / be, al / th, al / be, be / (al + th)];
        baseRowCells = table2cell(baseT(i, :));
        for r = 1:nRatios
            if ~isfinite(ratioVals(r))
                nNonFiniteRatios = nNonFiniteRatios + 1;
            end
            rowCount = rowCount + 1;
            rows(rowCount, :) = [baseRowCells, {ratioNames{r}}, {ratioVals(r)}];
        end
    end
    rows = rows(1:rowCount, :);

    varNames = [keyVars, {'RatioType', 'RatioValue'}];
    textColIdx = 1:(numel(keyVars) + 1);
    numColIdx = numel(keyVars) + 2;
    RatioLongT = local_build_table(rows, textColIdx, numColIdx, varNames);
end

function [Tout, nTaskMismatch, nBaselineExcluded] = local_filter_task_condition(Tin, taskKey, conditionOrder)
% Defensively filters an X11 table to TaskType==taskKey and Condition in
% conditionOrder. Any row NOT matching taskKey is dropped and counted
% (should be 0 if X11's own per-Task file contract holds); any Baseline
% (or other non-E/N/R) row is dropped and counted, mirroring X05/X08/X09's
% own Baseline-exclusion handling. See script header, JUDGMENT CALLS item 1.
    taskCol = cellstr(Tin.TaskType);
    taskMask = strcmp(taskCol, taskKey);
    nTaskMismatch = sum(~taskMask);
    Ttask = Tin(taskMask, :);
    condColTask = cellstr(Ttask.Condition);
    condMask = ismember(condColTask, conditionOrder);
    nBaselineExcluded = sum(~condMask);
    Tout = Ttask(condMask, :);
end

function [OmnibusT, PairwiseT] = local_fit_ratio_family(RatioLongT, unitColName, unitList, ratioTypes, variants, conditionOrder, refLevel, fdrAlpha)
% Fits RatioValue ~ Condition + (1|Subject) (REML, Satterthwaite) for
% every (Variant x RatioType x unit) combination, applies BH-FDR
% separately within each (Variant x RatioType) family across unitList
% (family size = numel(unitList), reduced only by any genuinely failed
% fit -- never padded with a placeholder), and computes Holm-corrected
% N-E/R-E/R-N pairwise contrasts for every unit whose omnibus test
% survives its own family's FDR. unitColName is '' for the Global level
% (single fit per Variant x RatioType, no grouping column in RatioLongT);
% the output tables still carry a 'Unit' column (constant 'Global') for a
% uniform report format at that level's trivial family size of 1.
    hasUnit = ~isempty(unitColName);
    if hasUnit
        unitLabelColName = unitColName;
    else
        unitLabelColName = 'Unit';
    end

    omnibusRows = cell(0, 7);   % Variant, RatioType, Unit, FStat, DF1, DF2, pValue
    lmeMap = containers.Map();

    for vi = 1:numel(variants)
        variantName = variants{vi};
        for ri = 1:numel(ratioTypes)
            ratioName = ratioTypes{ri};
            for ui = 1:numel(unitList)
                unitName = unitList{ui};
                if hasUnit
                    mask = strcmp(RatioLongT.Variant, variantName) & strcmp(RatioLongT.RatioType, ratioName) & strcmp(RatioLongT.(unitColName), unitName);
                else
                    mask = strcmp(RatioLongT.Variant, variantName) & strcmp(RatioLongT.RatioType, ratioName);
                end
                subT = RatioLongT(mask, :);
                if isempty(subT)
                    continue;
                end
                subT.Subject = categorical(local_subject_numeric(subT.Subject));
                subT.Condition = categorical(cellstr(subT.Condition), conditionOrder);
                try
                    lme = fitlme(subT, 'RatioValue ~ Condition + (1|Subject)', 'FitMethod', 'REML');
                catch
                    continue;   % failed fit: this unit is dropped from this family entirely (never a placeholder p-value)
                end
                aov = anova(lme, 'DFMethod', 'satterthwaite');
                rowIdx = find(strcmp(aov.Term, 'Condition'), 1);
                if isempty(rowIdx)
                    continue;
                end
                key = sprintf('%s|%s|%s', variantName, ratioName, unitName);
                lmeMap(key) = lme;
                omnibusRows(end+1, :) = {variantName, ratioName, unitName, aov.FStat(rowIdx), aov.DF1(rowIdx), aov.DF2(rowIdx), aov.pValue(rowIdx)}; %#ok<AGROW>
            end
        end
    end

    varNamesOmnibus = {'Variant', 'RatioType', unitLabelColName, 'FStat', 'DF1', 'DF2', 'pValue'};
    OmnibusT = local_build_table(omnibusRows, [1 2 3], [4 5 6 7], varNamesOmnibus);
    OmnibusT.pFDR = nan(height(OmnibusT), 1);
    OmnibusT.SignificantFDR = false(height(OmnibusT), 1);
    OmnibusT.FamilySize = nan(height(OmnibusT), 1);

    for vi = 1:numel(variants)
        variantName = variants{vi};
        for ri = 1:numel(ratioTypes)
            ratioName = ratioTypes{ri};
            idx = strcmp(OmnibusT.Variant, variantName) & strcmp(OmnibusT.RatioType, ratioName);
            nInFamily = sum(idx);
            if nInFamily == 0; continue; end
            OmnibusT.FamilySize(idx) = nInFamily;
            pvals = OmnibusT.pValue(idx);
            if nInFamily == 1
                OmnibusT.pFDR(idx) = pvals;
                OmnibusT.SignificantFDR(idx) = pvals <= fdrAlpha;
            else
                [adjP, reject] = local_fdr_bh(pvals, fdrAlpha);
                OmnibusT.pFDR(idx) = adjP;
                OmnibusT.SignificantFDR(idx) = reject;
            end
        end
    end

    pairwiseRows = cell(0, 12);   % Variant, RatioType, Unit, Contrast, Estimate, SE, DF, tStat, pValue, pHolm, CI_Low, CI_High
    CONTRASTS = {'N','E','N - E'; 'R','E','R - E'; 'R','N','R - N'};

    sigIdx = find(OmnibusT.SignificantFDR)';
    for row = sigIdx
        variantName = OmnibusT.Variant{row};
        ratioName = OmnibusT.RatioType{row};
        unitName = OmnibusT.(unitLabelColName){row};
        key = sprintf('%s|%s|%s', variantName, ratioName, unitName);
        if ~isKey(lmeMap, key); continue; end
        lme = lmeMap(key);
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
        [adjP, ~] = local_holm_bonferroni(contrastP, fdrAlpha);
        for c = 1:3
            pairwiseRows(end+1, :) = [{variantName, ratioName, unitName}, contrastData(c, 1:6), {adjP(c)}, contrastData(c, 7:8)]; %#ok<AGROW>
        end
    end
    varNamesPairwise = {'Variant', 'RatioType', unitLabelColName, 'Contrast', 'Estimate', 'SE', 'DF', 'tStat', 'pValue', 'pHolm', 'CI_Low', 'CI_High'};
    PairwiseT = local_build_table(pairwiseRows, [1 2 3 4], [5 6 7 8 9 10 11 12], varNamesPairwise);
end
