%% ========================================================================
% X13_asymmetry_indices_cognitive.m
%
% Purpose: hemispheric asymmetry indices, Alpha band only (exploratory
% extension family 3 of 5: X11-X15). Adds derived asymmetry metrics
% alongside the two locked co-primary metrics; reads X04's channel-level
% band power and X07's shared ROI/electrode definitions read-only, and
% nothing from X05/X06/X08/X09/X10's own statistical outputs.
%
% Method: standard Frontal Alpha Asymmetry (FAA) convention -- for a
% right/left electrode or region pair, Asymmetry = log(AlphaPower_right)
% - log(AlphaPower_left), in two parallel variants never pooled with each
% other: MeanBased (LogAbsolutePower = ln(mean AbsolutePower), X04's own
% convention) and MedianBased (log(AbsMedian), log-transformed here, from
% sibling script X11's MedianRobustElectrodeT output; X11's own AbsMedian
% is intentionally not log-transformed there). Alpha is the only band
% analyzed, since Frontal Alpha Asymmetry is inherently an Alpha-band
% concept. Three levels, each its own analysis: (1) Electrode-pair level
% -- 5 fixed right/left pairs in two groups never pooled for FDR: Frontal
% (F3/F4, AF3/AF4, F7/F8, the classic FAA pairs) and Posterior (O1/O2,
% PO7/PO8) -- the Posterior group is not a standard FAA pair set; it is
% an exploratory follow-up probe targeting a posterior effect noticed
% descriptively in earlier electrode-level results (X09), not itself an
% established finding. (2) ROI level -- for each of the 5 locked ROIs, if
% it has a usable bilateral split of member channels (qualifies only if
% it has at least one member channel on each hemisphere and
% min(nLeft,nRight)/max(nLeft,nRight) >= 0.5; non-qualifying ROIs are
% excluded and the reason is reported), Asymmetry = mean(log power,
% Right-hemisphere member channels) - mean(log power, Left-hemisphere
% member channels); family size = however many ROIs actually qualify
% (for this montage, all five). (3) Global level -- all 62 electrodes
% (P9/P10 included) split into Left/Right/Midline by the standard
% 10-20/10-10 trailing-digit convention (odd = Left, even = Right, 'z'/no
% digit = Midline, excluded from the split); one test per Variant (family
% size 1).
%
% Statistical model: Asymmetry ~ Condition + (1|Subject), REML,
% Satterthwaite DF, identical LMM machinery to
% X05_primary_cognitive_statistics.m. Condition in {E, N, R} only
% (Baseline rows excluded; counts logged per Task, per data source). Runs
% as three independent, parallel per-Task analyses (Full/Nback/SART).
%
% FDR family structure (own, separate families; never pooled with the
% locked X05/X06/X08/X09 metrics, X11's median/robust family, X12's ratio
% family, or each other), per Task: Pair x Frontal x MeanBased (family
% size 3), Pair x Frontal x MedianBased (3), Pair x Posterior x MeanBased
% (2), Pair x Posterior x MedianBased (2), ROI x MeanBased (# qualifying
% ROIs), ROI x MedianBased (# qualifying ROIs), Global x MeanBased
% (1, reported as-is), Global x MedianBased (1). Pairwise contrasts (N-E,
% R-E, R-N, matching X05/X06's convention) are computed only for units
% whose family-level omnibus test survives its own family's FDR,
% Holm-corrected within that unit (family size 3).
%
% Requirements: base MATLAB + Statistics and Machine Learning Toolbox
% (fitlme, anova on LinearMixedModel, coefTest, qqplot).
%
% Inputs (read-only, never modified):
%   <ANALYSIS_ROOT>/04_subject_level/phase3_subject_level_bandpower.mat
%     (BandPowerT -- channel-level, X04 Structure B) -- MeanBased source.
%   <ANALYSIS_ROOT>/07_roi_electrode_definitions/
%     phase5_channel_roi_definitions.mat (ChannelDefT, ROI_DEFINITION,
%     ROI_NAMES, UNASSIGNED_FOR_ROI) -- shared across all three Tasks.
%   <ANALYSIS_ROOT>/<Full|Nback|SART>/11_median_robust_family/<RUN_MODE>/
%     median_robust_data.mat (MedianRobustElectrodeT) -- produced by
%     sibling script X11_median_robust_family_cognitive.m, per Task; this
%     script errors immediately, per Task, if that file is missing.
%
% Outputs (one complete set per Task):
%   <ANALYSIS_ROOT>/<Full|Nback|SART>/13_asymmetry_indices/
%       phase13_asymmetry_values.mat / .csv / .xlsx   (underlying long
%           table of every Subject x Condition x Level x Group x Unit x
%           Variant asymmetry value used by every model below)
%       phase13_asymmetry_omnibus_results.mat / .xlsx / .csv
%       phase13_asymmetry_pairwise_contrasts.mat / .xlsx / .csv
%       phase13_parameters.mat
%       diagnostic_residuals_<Level>_<Group>_<Unit>_<Variant>.png
%   <ANALYSIS_ROOT>/09_reports/
%       phase_asymmetry_indices_report_<Full|Nback|SART>.txt
%
% RUN_MODE: full_sample / sensitivity_excl_noisy (see
% EXCLUSION_CONFIG_cognitive.m); run this script once per mode, with the
% same RUN_MODE value as the X11 run it reads. This script's own X04-side
% MeanBased data is exclusion-filtered directly; its X11 sibling-contract
% MedianBased data is already exclusion-filtered upstream.
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

bandPowerFile = fullfile(SUBJ_DIR, 'phase3_subject_level_bandpower.mat');
roiDefFile    = fullfile(DEF_DIR, 'phase5_channel_roi_definitions.mat');

if ~exist(bandPowerFile, 'file')
    error('Phase13:MissingInput', ['X04 band-power file not found: %s\n' ...
        'Run X04_aggregate_cognitive_subject_level.m first.'], bandPowerFile);
end
if ~exist(roiDefFile, 'file')
    error('Phase13:MissingInput', ['X07 ROI/electrode definition file not found: %s\n' ...
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
RUN_MODE = 'sensitivity_excl_noisy'; % or: 'sensitivity_excl_noisy'
if ~ismember(RUN_MODE, RUN_MODES_cognitive)
    error('Phase13:BadRunMode', 'RUN_MODE must be one of: %s', strjoin(RUN_MODES_cognitive, ', '));
end
RUN_MODE_LABEL = RUN_MODES_DISPLAY_cognitive.(RUN_MODE);
fprintf('RUN_MODE = %s\n', RUN_MODE_LABEL);
if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
    fprintf('Excluding subject(s) {%s} from every model in this run -- see EXCLUSION_CONFIG_cognitive.m for rationale.\n', ...
        strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','));
end
% ---- END NEW ----

if isempty(ver('stats'))
    error('Phase13:MissingToolbox', ...
        'Statistics and Machine Learning Toolbox not detected. This script requires fitlme, anova(LinearMixedModel), coefTest, and qqplot.');
end

fprintf('========================================================\n');
fprintf('PHASE 13 -- HEMISPHERIC ASYMMETRY INDICES (ALPHA), COGNITIVE-TASK DATA\n');
fprintf('Running as THREE INDEPENDENT, PARALLEL analyses (Full / Nback / SART).\n');
fprintf('========================================================\n\n');

%% Section 1: Load X04 (MeanBased source, channel-level) and X07 (shared
%  ROI/electrode definitions), both read-only, loaded ONCE before the Task
%  loop -- neither depends on Task.
S3 = load(bandPowerFile, 'BandPowerT', 'PARAMS3');
BandPowerT_ALL = S3.BandPowerT;
PARAMS3 = S3.PARAMS3;
clear S3;

requiredMeanCols = {'Subject','TaskType','BlockNumber','Condition','Channel','Band', ...
                     'AbsolutePower','LogAbsolutePower'};
missingMeanCols = requiredMeanCols(~ismember(requiredMeanCols, BandPowerT_ALL.Properties.VariableNames));
if ~isempty(missingMeanCols)
    error('Phase13:MissingColumns', 'BandPowerT is missing expected column(s): %s', strjoin(missingMeanCols, ', '));
end
if ~ismember('Alpha', fieldnames(PARAMS3.BandDefinitionsHz))
    error('Phase13:NoAlphaBand', 'PARAMS3.BandDefinitionsHz has no ''Alpha'' entry -- cannot compute Frontal Alpha Asymmetry.');
end

S5D = load(roiDefFile, 'ChannelDefT', 'ROI_DEFINITION', 'ROI_NAMES', 'UNASSIGNED_FOR_ROI', 'PARAMS5D');
ChannelDefT = S5D.ChannelDefT;
ROI_DEFINITION = S5D.ROI_DEFINITION;
ROI_NAMES = S5D.ROI_NAMES;
UNASSIGNED_FOR_ROI = S5D.UNASSIGNED_FOR_ROI;
PARAMS5D = S5D.PARAMS5D;
clear S5D;

CONDITION_ORDER = {'E', 'N', 'R'};   % Baseline excluded from every model in this script
REF_LEVEL = CONDITION_ORDER{1};
FDR_ALPHA = 0.05;
ROI_LATERALITY_MIN_RATIO = 0.5;      % "badly lopsided" exclusion threshold -- see header Section 4 note

TASK_TYPES = {'full', 'nback', 'sart'};
TASK_DISPLAY = struct('full', 'Full', 'nback', 'Nback', 'sart', 'SART');

%% Section 2: Fixed electrode-pair definitions (Task-invariant; verbatim
%  from the header -- classic FAA Frontal pairs + the Posterior pairs
%  described above). Columns: {Group, PairName, LeftElectrode,
%  RightElectrode}.
PAIR_DEFINITIONS = { ...
    'Frontal',   'F3_F4',    'F3',  'F4'; ...
    'Frontal',   'AF3_AF4',  'AF3', 'AF4'; ...
    'Frontal',   'F7_F8',    'F7',  'F8'; ...
    'Posterior', 'O1_O2',    'O1',  'O2'; ...
    'Posterior', 'PO7_PO8',  'PO7', 'PO8' ...
};
requiredPairChans = upper(unique([PAIR_DEFINITIONS(:,3); PAIR_DEFINITIONS(:,4)]));
missingPairChans = requiredPairChans(~ismember(requiredPairChans, upper(ChannelDefT.Label)));
if ~isempty(missingPairChans)
    error('Phase13:MissingPairElectrodes', ...
        'The following required pair electrode(s) are not present in X07''s ChannelDefT: %s', strjoin(missingPairChans, ', '));
end

%% Section 3: Electrode laterality classification (Task-invariant; applied
%  identically to every channel in the montage via local_electrode_
%  laterality -- see header for the exact convention).
ChannelDefT.Hemisphere = cellfun(@local_electrode_laterality, ChannelDefT.Label, 'UniformOutput', false);

nGlobalLeft = sum(strcmp(ChannelDefT.Hemisphere, 'Left'));
nGlobalRight = sum(strcmp(ChannelDefT.Hemisphere, 'Right'));
nGlobalMidline = sum(strcmp(ChannelDefT.Hemisphere, 'Midline'));
fprintf('--------------------------------------------------------\n');
fprintf('SECTION 3: GLOBAL (all 62-electrode) LATERALITY SPLIT (Task-invariant)\n');
fprintf('--------------------------------------------------------\n');
fprintf('  Left=%d  Right=%d  Midline(excluded)=%d  (total=%d)\n', nGlobalLeft, nGlobalRight, nGlobalMidline, height(ChannelDefT));
if nGlobalLeft == 0 || nGlobalRight == 0
    error('Phase13:GlobalLateralitySplitFailed', 'Global Left/Right split has zero channels on one hemisphere -- cannot compute a Global asymmetry index.');
end
GlobalLeftChans  = ChannelDefT.Label(strcmp(ChannelDefT.Hemisphere, 'Left'));
GlobalRightChans = ChannelDefT.Label(strcmp(ChannelDefT.Hemisphere, 'Right'));
GlobalLeftUC  = upper(GlobalLeftChans);
GlobalRightUC = upper(GlobalRightChans);

%% Section 4: ROI bilateral qualification (Task-invariant) -- for each of
%  the 5 locked ROIs, determine its member channels' Left/Right/Midline
%  split (from ROI_DEFINITION, the same locked membership X08 uses) and
%  apply the qualification rule from the header. Every ROI is evaluated;
%  none is assumed to qualify or excluded in advance.
fprintf('\n--------------------------------------------------------\n');
fprintf('SECTION 4: ROI BILATERAL QUALIFICATION (Task-invariant; threshold=%.2f)\n', ROI_LATERALITY_MIN_RATIO);
fprintf('--------------------------------------------------------\n');

QualifiedROIs = struct('Name', {}, 'LeftChannels', {}, 'RightChannels', {}, 'MidlineChannels', {}, ...
    'LeftChannelsUC', {}, 'RightChannelsUC', {}, 'NLeft', {}, 'NRight', {}, 'NMidline', {}, ...
    'Qualifies', {}, 'Reason', {});
for r = 1:numel(ROI_NAMES)
    roiName = ROI_NAMES{r};
    chans = ROI_DEFINITION.(roiName);
    hemis = cellfun(@local_electrode_laterality, chans, 'UniformOutput', false);
    leftChans = chans(strcmp(hemis, 'Left'));
    rightChans = chans(strcmp(hemis, 'Right'));
    midChans = chans(strcmp(hemis, 'Midline'));
    nL = numel(leftChans); nR = numel(rightChans); nM = numel(midChans);

    if nL == 0 || nR == 0
        qualifies = false;
        reason = sprintf('EXCLUDED: zero channels on one hemisphere (Left=%d, Right=%d).', nL, nR);
    else
        ratio = min(nL, nR) / max(nL, nR);
        if ratio < ROI_LATERALITY_MIN_RATIO
            qualifies = false;
            reason = sprintf('EXCLUDED: badly lopsided hemisphere counts (Left=%d, Right=%d, ratio=%.2f < %.2f threshold).', ...
                nL, nR, ratio, ROI_LATERALITY_MIN_RATIO);
        else
            qualifies = true;
            reason = sprintf('QUALIFIES: Left=%d, Right=%d, ratio=%.2f (Midline excluded=%d).', nL, nR, ratio, nM);
        end
    end
    QualifiedROIs(end+1) = struct('Name', roiName, 'LeftChannels', {leftChans}, 'RightChannels', {rightChans}, ...
        'MidlineChannels', {midChans}, 'LeftChannelsUC', {upper(leftChans)}, 'RightChannelsUC', {upper(rightChans)}, ...
        'NLeft', nL, 'NRight', nR, 'NMidline', nM, 'Qualifies', qualifies, 'Reason', reason); %#ok<AGROW>
    fprintf('  ROI=%-10s Left=%-3d Right=%-3d Midline=%-3d -> %s\n', roiName, nL, nR, nM, reason);
end
nQualifyingROIs = sum([QualifiedROIs.Qualifies]);
fprintf('\n%d / %d ROIs qualify for the ROI-level asymmetry split.\n\n', nQualifyingROIs, numel(QualifiedROIs));
if nQualifyingROIs == 0
    error('Phase13:NoQualifyingROIs', 'No ROI qualifies for a bilateral asymmetry split -- cannot run the ROI-level family.');
end

fprintf('\n########################################################\n');
fprintf('### Beginning per-Task asymmetry statistics (Full / Nback / SART)\n');
fprintf('########################################################\n\n');

%% ========================================================================
%  TASK LOOP -- everything below is repeated independently per Task
% ========================================================================
for taskI = 1:numel(TASK_TYPES)
    TASK_KEY = TASK_TYPES{taskI};
    TASK_LABEL = TASK_DISPLAY.(TASK_KEY);

    STATS_DIR = fullfile(ANALYSIS_ROOT, TASK_LABEL, '13_asymmetry_indices', RUN_MODE);
    if ~exist(STATS_DIR, 'dir'); mkdir(STATS_DIR); end

    fprintf('\n########################################################\n');
    fprintf('### TASK = %s | RUN_MODE = %s\n', TASK_LABEL, RUN_MODE_LABEL);
    fprintf('########################################################\n\n');

    integrityIssues = {};

    %% Section 5: Load and filter the MeanBased source (X04, Alpha band
    %  only, this Task, Condition in {E,N,R})
    thisTaskMeanT = BandPowerT_ALL(strcmp(BandPowerT_ALL.TaskType, TASK_KEY) & strcmp(BandPowerT_ALL.Band, 'Alpha'), :);
    nBaselineExcludedMean = sum(strcmp(thisTaskMeanT.Condition, 'Baseline'));
    MeanAlphaT = thisTaskMeanT(ismember(thisTaskMeanT.Condition, CONDITION_ORDER), :);
    MeanAlphaT.ChannelUC = upper(MeanAlphaT.Channel);
    fprintf('Task=%s [MeanBased/X04]: %d Alpha-band rows loaded, %d Baseline rows excluded, %d E/N/R rows retained.\n', ...
        TASK_LABEL, height(thisTaskMeanT), nBaselineExcludedMean, height(MeanAlphaT));

    % ---- NEW: apply the shared noisy-subject exclusion list (see
    % EXCLUSION_CONFIG_cognitive.m), only when RUN_MODE=sensitivity_excl_noisy.
    % Never applied to the full_sample run. Applied here to the MeanBased
    % (X04) source only; the MedianBased (X11) source below is already
    % exclusion-filtered upstream, since this script now reads X11's own
    % RUN_MODE-labeled output subfolder. ----
    nExcludedRowsSensitivityMean = 0;
    if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
        subjNumForExclMean = local_subject_numeric(MeanAlphaT.Subject);
        excludeMaskMean = ismember(subjNumForExclMean, EXCLUDE_NOISY_SUBJECTS);
        nExcludedRowsSensitivityMean = sum(excludeMaskMean);
        MeanAlphaT = MeanAlphaT(~excludeMaskMean, :);
        fprintf('Task=%s: RUN_MODE=%s -- excluded %d MeanBased/X04 row(s) belonging to subject(s) {%s}; %d rows remain.\n', ...
            TASK_LABEL, RUN_MODE_LABEL, nExcludedRowsSensitivityMean, ...
            strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','), height(MeanAlphaT));
    end

    %% Section 6: Load and filter the MedianBased source (X11 sibling
    %  contract, per-Task file, Alpha band only, Condition in {E,N,R}).
    %  Log-transform is applied HERE (X11's own AbsMedian is intentionally
    %  NOT log-transformed -- see header).
    % ---- X11 saves under a RUN_MODE-labeled subfolder too, so this path
    % must match it exactly --
    % this script must be run with the SAME RUN_MODE value as the X11 run
    % it is reading (X11 excludes noisy subjects at its own source, this
    % script never re-excludes them). ----
    medianFile = fullfile(ANALYSIS_ROOT, TASK_LABEL, '11_median_robust_family', RUN_MODE, 'median_robust_data.mat');
    if ~exist(medianFile, 'file')
        error('Phase13:MissingInput', ['X11 median/robust file not found for Task=%s: %s\n' ...
            'Run X11_median_robust_family_cognitive.m first (with the SAME RUN_MODE=%s). This ' ...
            'script trusts X11''s documented contract exactly and will NOT reconstruct or guess ' ...
            'median-based data itself.'], TASK_LABEL, medianFile, RUN_MODE);
    end
    S11 = load(medianFile, 'MedianRobustElectrodeT');
    if ~isfield(S11, 'MedianRobustElectrodeT')
        error('Phase13:MissingInput', 'MedianRobustElectrodeT not found in %s -- does not match the documented X11 contract.', medianFile);
    end
    MedianRobustElectrodeT = S11.MedianRobustElectrodeT;
    clear S11;
    requiredMedCols = {'Subject','BlockNumber','Condition','Channel','Band','AbsMedian'};
    missingMedCols = requiredMedCols(~ismember(requiredMedCols, MedianRobustElectrodeT.Properties.VariableNames));
    if ~isempty(missingMedCols)
        error('Phase13:MissingColumns', 'MedianRobustElectrodeT is missing expected column(s): %s (does not match the documented X11 contract).', strjoin(missingMedCols, ', '));
    end
    % TaskType filter is applied defensively (a judgment call, disclosed
    % here): the X11 contract file is already Task-scoped by its own save
    % path, but if it also carries a TaskType column (as the contract
    % lists), this script filters by it rather than assuming the whole
    % file is single-Task.
    if ismember('TaskType', MedianRobustElectrodeT.Properties.VariableNames)
        thisTaskMedT = MedianRobustElectrodeT(strcmp(MedianRobustElectrodeT.TaskType, TASK_KEY) & strcmp(MedianRobustElectrodeT.Band, 'Alpha'), :);
    else
        thisTaskMedT = MedianRobustElectrodeT(strcmp(MedianRobustElectrodeT.Band, 'Alpha'), :);
    end
    nBaselineExcludedMedian = sum(strcmp(thisTaskMedT.Condition, 'Baseline'));
    MedianAlphaT = thisTaskMedT(ismember(thisTaskMedT.Condition, CONDITION_ORDER), :);
    MedianAlphaT.ChannelUC = upper(MedianAlphaT.Channel);
    fprintf('Task=%s [MedianBased/X11]: %d Alpha-band rows loaded, %d Baseline rows excluded, %d E/N/R rows retained.\n', ...
        TASK_LABEL, height(thisTaskMedT), nBaselineExcludedMedian, height(MedianAlphaT));

    nNonPositiveMedian = sum(MedianAlphaT.AbsMedian <= 0);
    MedianAlphaT.LogAbsMedian = log(MedianAlphaT.AbsMedian);
    if nNonPositiveMedian > 0
        MedianAlphaT.LogAbsMedian(MedianAlphaT.AbsMedian <= 0) = NaN;
        integrityIssues{end+1} = sprintf('Task=%s: %d MedianBased AbsMedian value(s) were non-positive -- LogAbsMedian set to NaN there, not silently patched.', TASK_LABEL, nNonPositiveMedian); %#ok<AGROW>
    end

    if isempty(MeanAlphaT) || isempty(MedianAlphaT)
        warning('Task=%s: no E/N/R Alpha-band rows found in one or both sources -- skipping this Task entirely.', TASK_LABEL);
        continue;
    end

    %% Section 7: Structural presence check -- every channel required by
    %  the pair/ROI/global definitions above must be present in BOTH
    %  sources for this Task (a genuine structural requirement, hard-gated
    %  like X08/X09's channel-presence checks, not a guessed count).
    qualifyingROIChans = {};
    for r = 1:numel(QualifiedROIs)
        if QualifiedROIs(r).Qualifies
            qualifyingROIChans = [qualifyingROIChans, QualifiedROIs(r).LeftChannels, QualifiedROIs(r).RightChannels]; %#ok<AGROW>
        end
    end
    requiredChans = upper(unique([PAIR_DEFINITIONS(:,3)', PAIR_DEFINITIONS(:,4)', qualifyingROIChans, GlobalLeftChans', GlobalRightChans']));
    missingFromMean = requiredChans(~ismember(requiredChans, MeanAlphaT.ChannelUC));
    missingFromMedian = requiredChans(~ismember(requiredChans, MedianAlphaT.ChannelUC));
    if ~isempty(missingFromMean)
        error('Phase13:MissingChannelsMean', 'Task=%s: %d required channel(s) missing from the MeanBased (X04) Alpha-band data: %s', ...
            TASK_LABEL, numel(missingFromMean), strjoin(missingFromMean, ', '));
    end
    if ~isempty(missingFromMedian)
        error('Phase13:MissingChannelsMedian', 'Task=%s: %d required channel(s) missing from the MedianBased (X11) Alpha-band data: %s', ...
            TASK_LABEL, numel(missingFromMedian), strjoin(missingFromMedian, ', '));
    end
    fprintf('Task=%s: all %d required channels present in both MeanBased and MedianBased sources.\n', TASK_LABEL, numel(requiredChans));

    %% Section 8: Recording-level design verification (mirrors X05
    %  Section 2), computed from the MeanBased source (the always-present
    %  reference), with the MedianBased source's own recording count
    %  reported alongside for comparison/consistency, never assumed equal.
    comboT = unique(MeanAlphaT(:, {'Subject','BlockNumber','Condition'}), 'rows');
    nRecordings = height(comboT);
    comboMedT = unique(MedianAlphaT(:, {'Subject','BlockNumber','Condition'}), 'rows');
    nRecordingsMedian = height(comboMedT);
    fprintf('Section 8: %d unique Subject x BlockNumber x Condition combinations (MeanBased), %d (MedianBased) for Task=%s.\n', ...
        nRecordings, nRecordingsMedian, TASK_LABEL);
    if nRecordings ~= nRecordingsMedian
        integrityIssues{end+1} = sprintf('Task=%s: MeanBased (%d) and MedianBased (%d) sources have a different number of unique recordings -- the two Variants'' models will not necessarily be fit on identical N. Reported, not silently equalized.', ...
            TASK_LABEL, nRecordings, nRecordingsMedian); %#ok<AGROW>
    end

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
    pairKeys = keys(pairTypeCount);
    fprintf('Condition-pair coverage across %d subjects (Task=%s):\n', nSubjects, TASK_LABEL);
    for k = 1:numel(pairKeys)
        fprintf('  %s : %d subjects\n', strrep(pairKeys{k}, '_', ' & '), pairTypeCount(pairKeys{k}));
    end

    allSubjectsComplete3 = false;
    condsPerSubjCounts = cellfun(@numel, subjPairs);
    if all(condsPerSubjCounts >= 3)
        allSubjectsComplete3 = true;
        integrityIssues{end+1} = sprintf(['Task=%s: every subject appears to have 3+ recordings/conditions -- the incomplete-design ' ...
            'justification for using an LMM instead of RM-ANOVA no longer applies. STOP and re-examine.'], TASK_LABEL); %#ok<AGROW>
    end
    if allSubjectsComplete3
        error('Phase13:DesignAssumptionViolated', ...
            'Task=%s: design verification found a complete 3-condition-per-subject structure. Stopping before computing/fitting any asymmetry model -- see integrityIssues.', TASK_LABEL);
    end

    %% Section 9: Build the underlying asymmetry-value table (long format)
    %  for BOTH variants, at all three levels, using each variant's own
    %  set of recordings.
    rowsMean = local_build_asymmetry_rows(MeanAlphaT, 'LogAbsolutePower', comboT, PAIR_DEFINITIONS, QualifiedROIs, GlobalLeftUC, GlobalRightUC, 'MeanBased');
    rowsMedian = local_build_asymmetry_rows(MedianAlphaT, 'LogAbsMedian', comboMedT, PAIR_DEFINITIONS, QualifiedROIs, GlobalLeftUC, GlobalRightUC, 'MedianBased');

    asymCols = {'Subject','BlockNumber','Condition','Level','FamilyGroup','UnitName','Variant','Asymmetry'};
    AsymMeanT = local_build_table(rowsMean, [1 3 4 5 6 7], [2 8], asymCols);
    AsymMedianT = local_build_table(rowsMedian, [1 3 4 5 6 7], [2 8], asymCols);
    AsymmetryT = [AsymMeanT; AsymMedianT];

    nNaNAsym = sum(isnan(AsymmetryT.Asymmetry));
    if nNaNAsym > 0
        integrityIssues{end+1} = sprintf('Task=%s: %d asymmetry value(s) are NaN (propagated from a missing/non-positive underlying channel value) -- these rows will be list-wise excluded by fitlme, never imputed.', TASK_LABEL, nNaNAsym); %#ok<AGROW>
    end

    fprintf('\nAsymmetry-value table built (Task=%s): %d rows (%d Pair + %d ROI + %d Global units) x 2 Variants.\n', ...
        TASK_LABEL, height(AsymmetryT), size(PAIR_DEFINITIONS,1), nQualifyingROIs, 1);

    %% Section 10: Fit Asymmetry ~ Condition + (1|Subject) for every
    %  (Level, FamilyGroup, UnitName, Variant) unit
    unitKeys = unique(AsymmetryT(:, {'Level','FamilyGroup','UnitName','Variant'}), 'rows');
    nUnits = height(unitKeys);
    fprintf('\n--------------------------------------------------------\n');
    fprintf('SECTION 10: MODEL FITTING (Task=%s, %d units: Pair=%d, ROI=%d, Global=%d, x2 Variants)\n', ...
        TASK_LABEL, nUnits, size(PAIR_DEFINITIONS,1)*2, nQualifyingROIs*2, 1*2);
    fprintf('--------------------------------------------------------\n');

    omnibusRows = cell(0, 8);   % Level,FamilyGroup,UnitName,Variant,FStat,DF1,DF2,pValue
    lmeStore = containers.Map('KeyType', 'char', 'ValueType', 'any');
    convergenceWarnings = {};

    for u = 1:nUnits
        lvl = unitKeys.Level{u}; grp = unitKeys.FamilyGroup{u}; nm = unitKeys.UnitName{u}; variantName = unitKeys.Variant{u};
        subT = AsymmetryT(strcmp(AsymmetryT.Level, lvl) & strcmp(AsymmetryT.FamilyGroup, grp) & ...
                           strcmp(AsymmetryT.UnitName, nm) & strcmp(AsymmetryT.Variant, variantName), :);
        subT = subT(~isnan(subT.Asymmetry), :);
        subT.Subject = categorical(local_subject_numeric(subT.Subject));
        subT.Condition = categorical(cellstr(subT.Condition), CONDITION_ORDER);
        unitKey = sprintf('%s|%s|%s|%s', lvl, grp, nm, variantName);

        formula = 'Asymmetry ~ Condition + (1|Subject)';
        try
            lastwarn('');
            lme = fitlme(subT, formula, 'FitMethod', 'REML');
            [warnMsg, ~] = lastwarn();
            if ~isempty(warnMsg)
                convergenceWarnings{end+1} = sprintf('%s (Task=%s): %s', unitKey, TASK_LABEL, warnMsg); %#ok<AGROW>
            end
        catch ME
            integrityIssues{end+1} = sprintf('fitlme FAILED for %s (Task=%s): %s', unitKey, TASK_LABEL, ME.message); %#ok<AGROW>
            continue;
        end
        lmeStore(unitKey) = lme;
        aov = anova(lme, 'DFMethod', 'satterthwaite');
        rowIdx = find(strcmp(aov.Term, 'Condition'), 1);
        if isempty(rowIdx)
            integrityIssues{end+1} = sprintf('anova(lme) has no Condition row for %s (Task=%s).', unitKey, TASK_LABEL); %#ok<AGROW>
            continue;
        end
        omnibusRows(end+1, :) = {lvl, grp, nm, variantName, aov.FStat(rowIdx), aov.DF1(rowIdx), aov.DF2(rowIdx), aov.pValue(rowIdx)}; %#ok<AGROW>
    end

    OmnibusT = local_build_table(omnibusRows, [1 2 3 4], [5 6 7 8], {'Level','FamilyGroup','UnitName','Variant','FStat','DF1','DF2','pValue'});

    %% Section 11: BH-FDR within each of the 8 pre-specified families
    %  (Level x FamilyGroup x Variant), per Task
    OmnibusT.FamilyKey = strcat(OmnibusT.Level, '|', OmnibusT.FamilyGroup, '|', OmnibusT.Variant);
    OmnibusT.pFDR = nan(height(OmnibusT), 1);
    OmnibusT.SignificantFDR = false(height(OmnibusT), 1);
    OmnibusT.FamilySize = nan(height(OmnibusT), 1);

    uFamKeys = unique(OmnibusT.FamilyKey);
    for f = 1:numel(uFamKeys)
        idx = strcmp(OmnibusT.FamilyKey, uFamKeys{f});
        OmnibusT.FamilySize(idx) = sum(idx);
        pvals = OmnibusT.pValue(idx);
        if numel(pvals) == 1
            OmnibusT.pFDR(idx) = pvals;
            OmnibusT.SignificantFDR(idx) = pvals <= FDR_ALPHA;
        else
            [adjP, reject] = local_fdr_bh(pvals, FDR_ALPHA);
            OmnibusT.pFDR(idx) = adjP;
            OmnibusT.SignificantFDR(idx) = reject;
        end
    end
    fprintf('BH-FDR applied within %d families (expected 8: Pair x {Frontal,Posterior} x {MeanBased,MedianBased}, ROI x {MeanBased,MedianBased}, Global x {MeanBased,MedianBased}).\n', numel(uFamKeys));

    fprintf('\n--- Asymmetry omnibus results (Task=%s) ---\n', TASK_LABEL);
    disp(OmnibusT(:, {'Level','FamilyGroup','UnitName','Variant','FStat','DF1','DF2','pValue','pFDR','SignificantFDR','FamilySize'}));

    %% Section 12: Pairwise contrasts (N-E, R-E, R-N), gated on FDR,
    %  Holm-corrected within each surviving unit (family size 3)
    pairwiseRows = cell(0, 13);   % Level,FamilyGroup,UnitName,Variant,Contrast,Estimate,SE,DF,tStat,pValue,pHolm,CI_Low,CI_High
    CONTRASTS = {'N','E','N - E'; 'R','E','R - E'; 'R','N','R - N'};

    sigMask = OmnibusT.SignificantFDR;
    for row = find(sigMask)'
        lvl = OmnibusT.Level{row}; grp = OmnibusT.FamilyGroup{row}; nm = OmnibusT.UnitName{row}; variantName = OmnibusT.Variant{row};
        unitKey = sprintf('%s|%s|%s|%s', lvl, grp, nm, variantName);
        if ~isKey(lmeStore, unitKey); continue; end
        lme = lmeStore(unitKey);
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
            pairwiseRows(end+1, :) = [{lvl, grp, nm, variantName}, contrastData(c, 1:6), {adjP(c)}, contrastData(c, 7:8)]; %#ok<AGROW>
        end
    end
    PairwiseT = local_build_table(pairwiseRows, [1 2 3 4 5], [6 7 8 9 10 11 12 13], ...
        {'Level','FamilyGroup','UnitName','Variant','Contrast','Estimate','SE','DF','tStat','pValue','pHolm','CI_Low','CI_High'});

    fprintf('\n--- Pairwise contrasts (Task=%s; units significant under their own family''s FDR) ---\n', TASK_LABEL);
    disp(PairwiseT);

    %% Section 13: Residual diagnostics (every fitted unit x variant)
    nDiagSaved = 0;
    for u = 1:nUnits
        lvl = unitKeys.Level{u}; grp = unitKeys.FamilyGroup{u}; nm = unitKeys.UnitName{u}; variantName = unitKeys.Variant{u};
        unitKey = sprintf('%s|%s|%s|%s', lvl, grp, nm, variantName);
        if ~isKey(lmeStore, unitKey); continue; end
        lme = lmeStore(unitKey);
        try
            res = residuals(lme, 'ResidualType', 'Pearson');
            fig = figure('Visible', 'off', 'Position', [100 100 900 400]);
            subplot(1,2,1);
            histogram(res, 15);
            title(sprintf('Residuals: %s/%s/%s (%s), Task=%s', lvl, grp, nm, variantName, TASK_LABEL), 'Interpreter', 'none');
            xlabel('Pearson residual'); ylabel('Count');
            subplot(1,2,2);
            qqplot(res);
            title('Normal Q-Q');
            outPng = fullfile(STATS_DIR, sprintf('diagnostic_residuals_%s_%s_%s_%s.png', lvl, grp, nm, variantName));
            saveas(fig, outPng);
            close(fig);
            nDiagSaved = nDiagSaved + 1;
        catch ME
            integrityIssues{end+1} = sprintf('Residual diagnostic failed for %s (Task=%s): %s', unitKey, TASK_LABEL, ME.message); %#ok<AGROW>
        end
    end
    fprintf('Saved %d residual diagnostic plots (Task=%s).\n', nDiagSaved, TASK_LABEL);

    %% Section 14: PARAMS13 and outputs
    roiQualSummary = arrayfun(@(s) sprintf('%s: %s', s.Name, s.Reason), QualifiedROIs, 'UniformOutput', false);

    PARAMS13 = struct();
    PARAMS13.ScriptName = 'X13_asymmetry_indices_cognitive.m';
    PARAMS13.Task = TASK_LABEL;
    PARAMS13.GeneratedOn = datestr(now); %#ok<TNOW1,DATST>
    PARAMS13.AnalysisLevel = sprintf('Hemispheric asymmetry indices (log power, right-minus-left), Alpha band only, Task=%s only, Condition in {E,N,R} only (Baseline excluded). Purely additive exploratory addition; never pooled with X01-X10, X11, or X12.', TASK_LABEL);
    PARAMS13.Model = 'LinearMixedEffects_RandomIntercept';
    PARAMS13.FixedEffectsPrimary = 'Condition (3-level: E, N, R; reference = E)';
    PARAMS13.RandomEffects = '(1|Subject)';
    PARAMS13.FitMethod = 'REML';
    PARAMS13.DFMethod = 'Satterthwaite';
    PARAMS13.Variants = {'MeanBased (LogAbsolutePower, X04)', 'MedianBased (log(AbsMedian), X11 contract, log-transformed in this script)'};
    PARAMS13.PairDefinitions = PAIR_DEFINITIONS;
    PARAMS13.ROIQualificationRule = sprintf('A ROI qualifies only if it has >=1 channel on each hemisphere AND min(nLeft,nRight)/max(nLeft,nRight) >= %.2f. Excluded ROIs and their exact reason are disclosed below.', ROI_LATERALITY_MIN_RATIO);
    PARAMS13.ROIQualificationSummary = {roiQualSummary};
    PARAMS13.QualifyingROINames = {ROI_NAMES([QualifiedROIs.Qualifies])};
    PARAMS13.LateralityConvention = 'Trailing digit odd = Left, even = Right, ''z''/no digit = Midline (excluded from every Left/Right split).';
    PARAMS13.GlobalLateralitySplit = struct('NLeft', nGlobalLeft, 'NRight', nGlobalRight, 'NMidlineExcluded', nGlobalMidline);
    PARAMS13.FDRAlpha = FDR_ALPHA;
    PARAMS13.FDRFamilyScheme = 'Separate BH-FDR family per (Level x FamilyGroup x Variant): Pair x Frontal (n=3), Pair x Posterior (n=2), ROI (n=# qualifying ROIs), Global (n=1) -- each x 2 Variants = 8 families/Task. Never pooled with each other, with X05/X06/X08/X09/X10, with X11, or with X12.';
    PARAMS13.PairwiseContrastMethod = 'Holm-Bonferroni across the 3 pairwise contrasts (N-E, R-E, R-N -- matching X05/X06''s convention, standardized here for this new family rather than guessed), computed only for units whose own family''s omnibus test survives BH-FDR.';
    PARAMS13.PosteriorGroupNote = 'The Posterior pair group (O1/O2, PO7/PO8) targets a posterior effect noticed descriptively in earlier electrode-level (X09) results. It is NOT a standard FAA pair set and is NOT itself an established finding -- it is an exploratory follow-up probe.';
    PARAMS13.BaselineExclusionNote = sprintf('MeanBased: %d Baseline row(s) excluded; MedianBased: %d Baseline row(s) excluded, before any asymmetry value was computed.', nBaselineExcludedMean, nBaselineExcludedMedian);
    PARAMS13.MedianContractSource = medianFile;
    PARAMS13.SourcePhase3File = bandPowerFile;
    PARAMS13.SourcePhase5DFile = roiDefFile;
    PARAMS13.DesignVerification = struct('nRecordingsMeanBased', nRecordings, 'nRecordingsMedianBased', nRecordingsMedian, ...
        'nSubjects', nSubjects, 'PairTypeCounts', {pairTypeCount.keys}, 'IntegrityIssues', {integrityIssues});
    PARAMS13.RunDateTime = datestr(now); %#ok<TNOW1,DATST>
    PARAMS13.RunMode = RUN_MODE;
    if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
        PARAMS13.ExcludedSubjects = EXCLUDE_NOISY_SUBJECTS;
        PARAMS13.ExcludedSubjectsRationale = EXCLUSION_RATIONALE_cognitive;
        PARAMS13.RowsExcludedForSensitivity = nExcludedRowsSensitivityMean;
    else
        PARAMS13.ExcludedSubjects = [];
    end

    save(fullfile(STATS_DIR, 'phase13_parameters.mat'), 'PARAMS13');
    save(fullfile(STATS_DIR, 'phase13_asymmetry_values.mat'), 'AsymmetryT', 'PARAMS13');
    save(fullfile(STATS_DIR, 'phase13_asymmetry_omnibus_results.mat'), 'OmnibusT');
    save(fullfile(STATS_DIR, 'phase13_asymmetry_pairwise_contrasts.mat'), 'PairwiseT');

    writetable(AsymmetryT, fullfile(STATS_DIR, 'phase13_asymmetry_values.csv'));
    writetable(AsymmetryT, fullfile(STATS_DIR, 'phase13_asymmetry_values.xlsx'));
    writetable(OmnibusT, fullfile(STATS_DIR, 'phase13_asymmetry_omnibus_results.xlsx'));
    writetable(OmnibusT, fullfile(STATS_DIR, 'phase13_asymmetry_omnibus_results.csv'));
    if ~isempty(PairwiseT)
        writetable(PairwiseT, fullfile(STATS_DIR, 'phase13_asymmetry_pairwise_contrasts.xlsx'));
        writetable(PairwiseT, fullfile(STATS_DIR, 'phase13_asymmetry_pairwise_contrasts.csv'));
    end

    fprintf('\nSaved output (Task=%s) to: %s\n', TASK_LABEL, STATS_DIR);

    %% Section 15: Report
    reportLines = {};
    reportLines{end+1} = repmat('=', 1, 60);
    reportLines{end+1} = sprintf('PHASE 13 -- HEMISPHERIC ASYMMETRY INDICES REPORT (Task=%s)', TASK_LABEL);
    reportLines{end+1} = repmat('=', 1, 60);
    reportLines{end+1} = '';
    reportLines{end+1} = 'Purely additive exploratory addition. Never';
    reportLines{end+1} = 'modifies or depends on X01-X10''s locked statistical outputs.';
    reportLines{end+1} = 'Convention: Asymmetry = log(AlphaPower_right) - log(AlphaPower_left) (standard';
    reportLines{end+1} = 'Frontal Alpha Asymmetry convention). Alpha band only.';
    reportLines{end+1} = 'Posterior pair group (O1/O2, PO7/PO8) targets an effect noticed';
    reportLines{end+1} = 'descriptively in earlier electrode-level results --';
    reportLines{end+1} = 'it is NOT itself an established finding.';
    reportLines{end+1} = '';
    reportLines{end+1} = sprintf('Baseline rows excluded before computing asymmetry: MeanBased=%d, MedianBased=%d', nBaselineExcludedMean, nBaselineExcludedMedian);
    reportLines{end+1} = sprintf('Recordings verified: MeanBased=%d, MedianBased=%d | Subjects: %d', nRecordings, nRecordingsMedian, nSubjects);
    reportLines{end+1} = sprintf('Condition-pair coverage: %s', strjoin(cellfun(@(k) sprintf('%s=%d', strrep(k,'_',' & '), pairTypeCount(k)), pairKeys, 'UniformOutput', false), ', '));
    reportLines{end+1} = '';
    reportLines{end+1} = 'Model: Linear mixed-effects, Condition (fixed, 3-level) + (1|Subject), REML, Satterthwaite DF.';
    reportLines{end+1} = 'Variants: MeanBased (LogAbsolutePower, X04) and MedianBased (log(AbsMedian), X11), never pooled.';
    reportLines{end+1} = '';
    reportLines{end+1} = '--- ROI bilateral qualification (Task-invariant; same for every Task) ---';
    for r = 1:numel(QualifiedROIs)
        reportLines{end+1} = sprintf('  %s', roiQualSummary{r}); %#ok<AGROW>
    end
    reportLines{end+1} = sprintf('  --> %d / %d ROIs qualify for the ROI-level family.', nQualifyingROIs, numel(QualifiedROIs));
    reportLines{end+1} = '';
    reportLines{end+1} = sprintf('--- Global laterality split (all 62 electrodes, P9/P10 included) ---');
    reportLines{end+1} = sprintf('  Left=%d  Right=%d  Midline(excluded)=%d', nGlobalLeft, nGlobalRight, nGlobalMidline);
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Asymmetry omnibus results (8 FDR families: Pair x {Frontal,Posterior} x {Mean,Median}, ROI x {Mean,Median}, Global x {Mean,Median}) ---';
    for r = 1:height(OmnibusT)
        reportLines{end+1} = sprintf('  [%-6s|%-10s|%-10s|%-11s] F(%g,%.1f)=%.3f p=%.4f pFDR=%.4f family_n=%d %s', ...
            OmnibusT.Level{r}, OmnibusT.FamilyGroup{r}, OmnibusT.UnitName{r}, OmnibusT.Variant{r}, OmnibusT.DF1(r), OmnibusT.DF2(r), ...
            OmnibusT.FStat(r), OmnibusT.pValue(r), OmnibusT.pFDR(r), OmnibusT.FamilySize(r), local_yn(OmnibusT.SignificantFDR(r))); %#ok<AGROW>
    end
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Pairwise contrasts (N-E, R-E, R-N; units significant under their own family''s FDR, Holm-corrected within unit) ---';
    if isempty(PairwiseT)
        reportLines{end+1} = '  (none -- no unit survived its family''s FDR correction)';
    else
        for r = 1:height(PairwiseT)
            reportLines{end+1} = sprintf('  [%-6s|%-10s|%-10s|%-11s] %-8s est=%.4f (95%% CI %.4f to %.4f) t=%.3f df=%.1f p=%.4f pHolm=%.4f', ...
                PairwiseT.Level{r}, PairwiseT.FamilyGroup{r}, PairwiseT.UnitName{r}, PairwiseT.Variant{r}, PairwiseT.Contrast{r}, ...
                PairwiseT.Estimate(r), PairwiseT.CI_Low(r), PairwiseT.CI_High(r), PairwiseT.tStat(r), PairwiseT.DF(r), ...
                PairwiseT.pValue(r), PairwiseT.pHolm(r)); %#ok<AGROW>
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
    fid = fopen(fullfile(REPORT_DIR, sprintf('phase_asymmetry_indices_report_%s_%s.txt', TASK_LABEL, RUN_MODE_LABEL)), 'w');
    fprintf(fid, '%s\n', reportText);
    fclose(fid);

end % TASK_TYPES loop

fprintf('\n=== X13 complete for all three Tasks. Purely additive -- X01-X10''s outputs were never read or modified. ===\n');

%% ======================= LOCAL FUNCTIONS (all after top-level code) =======================

function hemi = local_electrode_laterality(label)
% Standard 10-20/10-10 laterality convention: a trailing digit that is
% odd -> Left hemisphere, even -> Right hemisphere; a 'z'/'Z' suffix or no
% trailing digit at all -> Midline. Applied identically to every channel
% label in this script (ROI-level and Global-level splits).
    lbl = char(label);
    tok = regexp(lbl, '(\d+)$', 'tokens', 'once');
    if isempty(tok)
        hemi = 'Midline';
        return;
    end
    d = str2double(tok{1});
    if mod(d, 2) == 1
        hemi = 'Left';
    else
        hemi = 'Right';
    end
end

function val = local_lookup_value(T, subjStr, blockNum, condStr, chanUCList, valueCol)
% Returns the value of T.(valueCol) for the given Subject x BlockNumber x
% Condition and channel(s) (matched via the precomputed ChannelUC column).
% A single-channel chanUCList (electrode-pair use) returns that channel's
% value; a multi-channel chanUCList (ROI/Global use) returns the
% arithmetic mean across those channels, matching the header's "mean(log
% power over ... channels)" definition for ROI/Global asymmetry. Returns
% NaN (never silently zero) if no matching row is found, and this is
% reported upstream via the "asymmetry values are NaN" integrity check.
    mask = strcmp(T.Subject, subjStr) & T.BlockNumber == blockNum & strcmp(T.Condition, condStr) & ismember(T.ChannelUC, chanUCList);
    vals = T.(valueCol)(mask);
    if isempty(vals)
        val = NaN;
    else
        val = mean(vals);
    end
end

function rows = local_build_asymmetry_rows(T, valueCol, recCombos, pairDefs, qualifiedROIs, globalLeftUC, globalRightUC, variantName)
% Builds long-format asymmetry rows (Subject,BlockNumber,Condition,Level,
% FamilyGroup,UnitName,Variant,Asymmetry) for one Variant, at all three
% levels (Pair, ROI, Global), for every recording in recCombos.
    rows = cell(0, 8);
    for i = 1:height(recCombos)
        subjStr = recCombos.Subject{i};
        blockNum = recCombos.BlockNumber(i);
        condStr = recCombos.Condition{i};

        % --- Electrode-pair level ---
        for p = 1:size(pairDefs, 1)
            grp = pairDefs{p, 1}; nm = pairDefs{p, 2}; leftE = pairDefs{p, 3}; rightE = pairDefs{p, 4};
            rVal = local_lookup_value(T, subjStr, blockNum, condStr, {upper(rightE)}, valueCol);
            lVal = local_lookup_value(T, subjStr, blockNum, condStr, {upper(leftE)}, valueCol);
            rows(end+1, :) = {subjStr, blockNum, condStr, 'Pair', grp, nm, variantName, rVal - lVal}; %#ok<AGROW>
        end

        % --- ROI level (qualifying ROIs only) ---
        for r = 1:numel(qualifiedROIs)
            if ~qualifiedROIs(r).Qualifies; continue; end
            rVal = local_lookup_value(T, subjStr, blockNum, condStr, qualifiedROIs(r).RightChannelsUC, valueCol);
            lVal = local_lookup_value(T, subjStr, blockNum, condStr, qualifiedROIs(r).LeftChannelsUC, valueCol);
            rows(end+1, :) = {subjStr, blockNum, condStr, 'ROI', 'AllROIs', qualifiedROIs(r).Name, variantName, rVal - lVal}; %#ok<AGROW>
        end

        % --- Global level ---
        rVal = local_lookup_value(T, subjStr, blockNum, condStr, globalRightUC, valueCol);
        lVal = local_lookup_value(T, subjStr, blockNum, condStr, globalLeftUC, valueCol);
        rows(end+1, :) = {subjStr, blockNum, condStr, 'Global', 'Global', 'Global', variantName, rVal - lVal}; %#ok<AGROW>
    end
end

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
            error('Phase13:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name1);
        end
        cRow(idx1) = cRow(idx1) - 1;
    end
    if ~strcmp(level2, refLevel)
        name2 = sprintf('%s_%s', varName, level2);
        idx2 = find(strcmp(coeffNames, name2));
        if isempty(idx2)
            error('Phase13:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name2);
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
        error('Phase13:BuildTableColumnMismatch', ...
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
            error('Phase13:BuildTableColumnMismatch', 'Column %d not classified as text or numeric.', j);
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
