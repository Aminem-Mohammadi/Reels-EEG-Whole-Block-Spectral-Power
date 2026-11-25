%% ========================================================================
% X16_full_reconstructed_cognitive.m
%
% Purpose: "FullReconstructed" exploratory, purely additive analysis,
% cognitive-task data. Single pseudo-Task = 'FullReconstructed'; there is
% NO Task loop in this script because, unlike X05/X06/X08/X09/X11 (each of
% which runs three independent Task analyses: Full/Nback/SART), this
% script constructs and analyzes exactly ONE new synthetic recording type.
%
% RUN_MODE: supports a toggle ('full_sample' / 'sensitivity_excl_noisy')
% near the top of Section 0, reading the single shared exclusion vector
% from EXCLUSION_CONFIG_cognitive.m. Run this script twice (once per
% RUN_MODE value) to get the full-sample results and the
% noisy-subjects-excluded sensitivity results side by side, in a
% RUN_MODE-labeled output subfolder. The exclusion filter is applied
% ONCE, right after this script builds its pooled Nback+SART band-power
% table (Section 3/4), so the SAME filtered BandPowerT/GlobalBandPowerT
% feed all three analysis sections below (Global, ROI, Electrode) --
% consistent with every other section of this script already sharing
% those two tables. Nothing about the pooling, aggregation, model, FDR
% scheme, or pairwise contrasts is changed by this toggle -- see
% EXCLUSION_CONFIG_cognitive.m for exactly which subjects are excluded in
% sensitivity mode and why.
%
% This is a new, purely additive script. It does NOT
% replace, modify, re-derive, or touch in any way the existing Full-task
% analysis already produced by X05/X06/X08/X09 for TaskType='full' (i.e.
% the raw, continuously-recorded Full session). Nothing under
% <ANALYSIS_ROOT>/Full/... is ever read or written by this script.
%
% CONCEPT (background from a separate diagnostic,
% DIAG01_full_vs_subtasks_duration_check.m): the raw "Full" task recording
% contains extra non-task time (inter-task transitions, instructions) that
% dilutes its band-power estimates relative to Nback/SART analyzed
% separately. This script builds an ALTERNATIVE, parallel definition of
% "Full" -- called "FullReconstructed" -- by taking the UNION of each
% Subject x BlockNumber x Condition's own ALREADY-EXISTING Nback and SART
% valid analysis windows (read directly from X03's psd_data, NEVER from
% the raw Full recording) and aggregating them together as if they were
% one recording. This is IN ADDITION TO the existing Full analysis, never
% a replacement.
%
% METHODOLOGY (every choice below is either an exact replication of an
% established pattern elsewhere in this codebase, or an explicitly
% disclosed deviation with a stated reason -- see "JUDGMENT CALLS" below
% and the inline comments at each point):
%   1. WINDOW POOLING (Sections 2-3): for every Subject x BlockNumber x
%      Condition that has BOTH a TaskType='nback' AND a TaskType='sart'
%      record in X03's psd_data, this script vertically concatenates the
%      Nback record's AbsoluteBandPower.(band)(WindowValid,:) with the
%      SART record's AbsoluteBandPower.(band)(WindowValid,:) (same for
%      RelativeBandPower), for every one of the 9 bands in
%      PARAMS.BandDefinitionsHz (generic loop over fieldnames -- band
%      names are never hard-coded except for the tier lists needed for
%      FDR family labeling, exactly as instructed). A combination with
%      only ONE of the two Task recordings is SKIPPED (never fabricated)
%      and logged as an integrity issue.
%   2. AGGREGATION MATH (Section 3): applied EXACTLY as
%      X04_aggregate_cognitive_subject_level.m -- arithmetic mean across
%      the pooled windows; LogAbsolutePower = log(mean(AbsolutePower))
%      (NEVER mean(log(...))); non-positive AbsolutePower -> NaN
%      LogAbsolutePower, logged as an integrity issue, never silently
%      patched; whole-scalp Global AbsolutePower/RelativePower = arithmetic
%      mean across channels of the per-channel means (X04 Structure C),
%      with the same non-positive -> NaN handling. This produces
%      BandPowerT_FullRecon (channel-level, X04 Structure B analog) and
%      GlobalBandPowerT_FullRecon (whole-scalp, X04 Structure C analog),
%      for TaskType='FullReconstructed', across ALL conditions present in
%      the pooled data (Baseline included, exactly like X04's own
%      "aggregate everything, exclude Baseline only downstream" policy).
%   3. STATISTICAL MODEL, THREE LEVELS (Sections 5-7), each an exact
%      structural copy of the corresponding existing script, applied to
%      this one new pseudo-Task:
%        - GLOBAL level  : copies X05_primary_cognitive_statistics.m's
%          model (Metric ~ Condition + (1|Subject), REML, Satterthwaite),
%          tiered FDR family scheme (BH-FDR per Tier x Metric; Primary
%          tier n=1 reported as-is, never treated as pre-corrected), and
%          Global-level residual-diagnostic scope (Primary+Secondary tier
%          bands, both metrics).
%        - ROI level     : copies X08_roi_level_statistics_cognitive.m's
%          ROI aggregation rule (arithmetic mean of channel-level
%          AbsolutePower/RelativePower across each ROI's member channels;
%          LogAbsolutePower = log(mean AbsolutePower)), its locked ROI
%          definitions (read READ-ONLY from X07's shared output, never
%          re-derived), its FDR family scheme (BH-FDR per Band x Metric,
%          across the 5 ROIs), and its residual-diagnostic scope
%          (Primary+Secondary tier bands, all 5 ROIs, both metrics).
%        - ELECTRODE level: copies
%          X09_electrode_level_statistics_cognitive.m's per-electrode
%          model (all 62 electrodes, P9/P10 included, no aggregation of
%          any kind), its FDR family scheme (BH-FDR per Band x Metric,
%          across the 62 electrodes), and its (lack of) residual-
%          diagnostic output -- X09 produces none, and neither does this
%          script at Electrode level, per the same tractability rationale.
%      PAIRWISE-CONTRAST CONVENTION -- ALL THREE
%      levels in this script use the N-E / R-E / R-N convention (as
%      defined and used by X05/X06 and by every exploratory-family script
%      X11-X15), Holm-corrected (family size 3), gated on that row's own
%      FDR-family survival. This is DELIBERATELY NOT the E-N/E-R/N-R
%      convention that X08's and X09's OWN original scripts use for
%      themselves -- this script inherits the model/FDR/ROI-aggregation/
%      electrode-set MACHINERY from X08/X09 but the PAIRWISE LABELING
%      CONVENTION from X05/X11-X15, by design. Both conventions coexist
%      elsewhere in this codebase; this script never conflates them.
%   4. STRUCTURAL TEMPLATE (overall script organization: header, Section 0
%      paths, loading inputs, computing, FDR, pairwise, saving, report):
%      copied from X11_median_robust_family_cognitive.m, the codebase's
%      own worked example of how to add a brand-new, self-contained,
%      purely-additive script. In particular, the fit/FDR/pairwise helper
%      functions (local_fit_condition_omnibus, local_bh_fdr_or_raw,
%      local_run_pairwise) are adapted directly from X11's own equivalents
%      -- X11 factored these out for the same reason this script needs
%      them: the same fit-omnibus-then-FDR-then-pairwise logic must be
%      repeated across multiple separate analysis contexts (here: Global /
%      ROI / Electrode) within one file.
%
% RELATIONSHIP TO OTHER SCRIPTS -- this is a new, stand-alone script, not
% a modification of any existing one. It has no single direct ancestor:
% it is a combination of X04's aggregation math + X05's Global model/FDR
% + X08's ROI model/FDR/aggregation + X09's Electrode model/FDR + X11's
% overall file structure. Every one of those source behaviors is copied
% EXACTLY as documented in items 1-4 above; nothing is invented beyond
% what those four scripts already establish, except the union/pooling
% step itself (item 1), which is the one genuinely new piece of
% methodology in this family.
%
% JUDGMENT CALLS (disclosed here and in PARAMS16.JudgmentCalls / the
% report, per this project's explicit "disclose, never silently decide"
% culture):
%   1. Channel-order safety check before pooling: X04's own channel-label
%      consistency check (Section 3 of X04) only LOGS a mismatch, never
%      skips, because X04 never concatenates two different recordings'
%      window matrices column-for-column. This script DOES (Section 3
%      vertically stacks Nback's and SART's window x channel matrices), so
%      a channel-order mismatch between the two recordings would silently
%      misalign channels if not caught. This script therefore SKIPS (and
%      logs) any Subject x BlockNumber x Condition combination where
%      either the Nback or the SART record's ChannelLabels does not
%      exactly equal X03's own ChannelLabelsReference -- a stricter,
%      safety-motivated deviation from X04's log-only policy, disclosed
%      here rather than silently applied.
%   2. Zero-pooled-window combinations: if the union of the Nback and SART
%      valid-window counts for a combination is 0 (both recordings had 0
%      valid windows, or the only available one did), that combination is
%      skipped for aggregation and logged -- exactly X04's own "0 valid
%      windows -> skip aggregation, still record the attempt" policy
%      (Section 3 of X04), extended to the pooled case.
%   3. A combination where ONE of the two Task recordings has 0 valid
%      windows but the OTHER has >0 is NOT skipped -- the union still
%      proceeds using whichever recording actually has valid windows. This
%      is consistent with "use ALL available valid windows, never
%      fabricate, never require more than actually exists" (X03's own
%      DurationHandlingPolicy), but is disclosed explicitly here since it
%      was not spelled out in the task description.
%   4. Electrode-level Item B (X09's own hard-coded 62-entry canonical
%      electrode-order literal, cross-checked independently against X07's
%      output) is NOT re-implemented here. This script validates the
%      electrode set it uses (ChannelDefT.Label, read READ-ONLY from X07)
%      only for count (62) and P9/P10 presence, relying on X07 itself as
%      the single locked source of truth, rather than re-hard-coding a
%      second, independent 62-item literal purely to cross-check against
%      X07 a second time. This is a disclosed scope reduction relative to
%      X09's belt-and-suspenders double-check, chosen to avoid introducing
%      a second hand-typed 62-electrode list (and its own transcription-
%      error risk) into a brand-new script.
%   5. Residual-diagnostic scope: Global level mirrors X05 exactly
%      (Primary+Secondary tier bands, both metrics); ROI level mirrors X08
%      exactly (Primary+Secondary tier bands, all 5 ROIs, both metrics);
%      Electrode level produces NO residual diagnostics, mirroring X09
%      (which also produces none). This deliberately does NOT follow
%      X11's own precedent of restricting ALL diagnostics to Global level
%      only -- the task instructions ask this script to mirror
%      05/08/09_..._stats' own OUTPUT CONTENT specifically, which takes
%      precedence over X11's tractability-driven narrowing for this one
%      aspect.
%   6. An extra output file not named in the literal task description,
%      full_reconstructed_bandpower.mat (BandPowerT_FullRecon,
%      GlobalBandPowerT_FullRecon, both BEFORE Baseline exclusion), is
%      saved alongside the requested global/roi/electrode result files.
%      This is the direct analog of X04's own phase3_subject_level_
%      bandpower.mat, saved for the same reason X04 saves it: so the
%      exact pooled/aggregated input to every model in this script can be
%      independently inspected and reproduced without re-running the
%      whole script. Disclosed here as an addition beyond the literal
%      file list, not a silent one.
%
% REQUIREMENTS: base MATLAB + Statistics and Machine Learning Toolbox
% (fitlme, anova on LinearMixedModel, coefTest).
%
% INPUT (read-only; nothing here is ever written back to X03/X04/X07's own
% output folders, and the existing Full-task X05/X06/X08/X09 output for
% TaskType='full' is never read or written):
%   <ANALYSIS_ROOT>/03_psd/cognitive_psd.mat (psd_data, FrequencyVectorHz,
%     ChannelLabelsReference, PARAMS) -- produced by X03.
%   <ANALYSIS_ROOT>/07_roi_electrode_definitions/
%     phase5_channel_roi_definitions.mat (ChannelDefT, ROI_DEFINITION,
%     ROI_NAMES, UNASSIGNED_FOR_ROI, PARAMS5D) -- produced by X07, SHARED
%     across all Tasks.
%
% OUTPUT (created fresh):
%   <ANALYSIS_ROOT>/FullReconstructed/16_full_reconstructed/
%       full_reconstructed_bandpower.mat   (see judgment call 6)
%       global_omnibus_results.mat / .xlsx
%       global_pairwise_contrasts.mat / .xlsx
%       roi_omnibus_results.mat / .xlsx / .csv
%       roi_pairwise_contrasts.mat / .xlsx / .csv
%       electrode_omnibus_results.mat / .xlsx / .csv
%       electrode_pairwise_contrasts.mat / .xlsx / .csv
%       parameters.mat   (PARAMS16)
%       diagnostic_residuals_Global_<band>_<metric>.png
%       diagnostic_residuals_ROI_<ROI>_<band>_<metric>.png
%   <ANALYSIS_ROOT>/09_reports/
%       X16_full_reconstructed_report.txt
%
% This script performs NO interpretation beyond the validation/statistics
% specified above, and does not compare FullReconstructed to the existing
% Full/Nback/SART analyses -- that is left for a separate, not-yet-written
% comparison step.
% ========================================================================

clear; clc;

%% ------------------------------------------------------------------
%  SECTION 0: PATHS AND EXISTENCE CHECKS
% ------------------------------------------------------------------

config;
DATA_ROOT_FOR_SIBLING = DATA_DIR;
[DATA_ROOT_PARENT, ~] = fileparts(DATA_ROOT_FOR_SIBLING);
ANALYSIS_ROOT = OUT_DIR;

PSD_DIR    = fullfile(ANALYSIS_ROOT, '03_psd');
DEF_DIR    = fullfile(ANALYSIS_ROOT, '07_roi_electrode_definitions');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');
STATS_DIR_BASE = fullfile(ANALYSIS_ROOT, 'FullReconstructed', '16_full_reconstructed');

psdMatFile = fullfile(PSD_DIR, 'cognitive_psd.mat');
roiDefFile = fullfile(DEF_DIR, 'phase5_channel_roi_definitions.mat');

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
    error('Phase16:BadRunMode', 'RUN_MODE must be one of: %s', strjoin(RUN_MODES_cognitive, ', '));
end
RUN_MODE_LABEL = RUN_MODES_DISPLAY_cognitive.(RUN_MODE);
fprintf('RUN_MODE = %s\n', RUN_MODE_LABEL);
if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
    fprintf('Excluding subject(s) {%s} from every model in this run -- see EXCLUSION_CONFIG_cognitive.m for rationale.\n', ...
        strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','));
end
% ---- END NEW ----

% NOTE: unlike X05 (which has one output dir PER Task, inside a Task
% loop), X16 has a SINGLE pseudo-Task and its Global/ROI/Electrode
% sections all already share ONE output directory (STATS_DIR below) --
% there is no separate per-section subfolder to begin with. RUN_MODE is
% therefore appended once here, and every save/writetable call in every
% section (Global/ROI/Electrode) below lands in the same RUN_MODE-labeled
% folder automatically, since they all reference this same STATS_DIR.
STATS_DIR = fullfile(STATS_DIR_BASE, RUN_MODE);
if ~exist(STATS_DIR, 'dir'); mkdir(STATS_DIR); end

fprintf('========================================================\n');
fprintf('PHASE 16 (X16) -- "FULL RECONSTRUCTED" EXPLORATORY ADDITIVE ANALYSIS\n');
fprintf('Single pseudo-Task = FullReconstructed (union of Nback+SART valid windows).\n');
fprintf('This is PURELY ADDITIVE. The existing Full-task analysis (X05/X06/X08/X09,\n');
fprintf('TaskType=''full'') is never read or written by this script.\n');
fprintf('========================================================\n');
fprintf('ANALYSIS_ROOT : %s\n', ANALYSIS_ROOT);
fprintf('OUTPUT        : %s\n', STATS_DIR);
fprintf('--------------------------------------------------------\n\n');

if ~exist(psdMatFile, 'file')
    error('Phase16:MissingInput', ['X03 output not found: %s\n' ...
        'Run X03_compute_cognitive_psd.m first.'], psdMatFile);
end
if ~exist(roiDefFile, 'file')
    error('Phase16:MissingInput', ['X07 output not found: %s\n' ...
        'Run X07_roi_electrode_definitions_cognitive.m first. This script never ' ...
        'reconstructs or guesses ROI/electrode membership on its own.'], roiDefFile);
end
if isempty(ver('stats'))
    error('Phase16:MissingToolbox', ...
        'Statistics and Machine Learning Toolbox not detected. This script requires fitlme, anova(LinearMixedModel), and coefTest.');
end

%% ------------------------------------------------------------------
%  SECTION 1: LOAD X03 (psd_data) AND X07 (locked ROI/electrode defs),
%  BOTH READ-ONLY.
% ------------------------------------------------------------------

S3 = load(psdMatFile, 'psd_data', 'FrequencyVectorHz', 'ChannelLabelsReference', 'PARAMS');
psd_data = S3.psd_data;
freqVecRef = S3.FrequencyVectorHz(:);
ChannelLabelsReference = S3.ChannelLabelsReference(:);
PARAMS = S3.PARAMS;
clear S3;

S7 = load(roiDefFile, 'ChannelDefT', 'ROI_DEFINITION', 'ROI_NAMES', 'UNASSIGNED_FOR_ROI', 'PARAMS5D');
ChannelDefT = S7.ChannelDefT;
ROI_DEFINITION = S7.ROI_DEFINITION;
ROI_NAMES = S7.ROI_NAMES;
UNASSIGNED_FOR_ROI = S7.UNASSIGNED_FOR_ROI;
PARAMS5D = S7.PARAMS5D;
clear S7;

ALL_BANDS = fieldnames(PARAMS.BandDefinitionsHz);
nBands = numel(ALL_BANDS);
nChanRef = numel(ChannelLabelsReference);

PRIMARY_TIER_BANDS     = {'Alpha'};
SECONDARY_TIER_BANDS   = {'Delta', 'Theta', 'Beta'};
EXPLORATORY_TIER_BANDS = {'Beta_low_exploratory', 'Beta_high_exploratory', ...
    'Gamma_low_exploratory', 'Gamma_high1_exploratory', 'Gamma_high2_exploratory'};
expectedTierBands = [PRIMARY_TIER_BANDS, SECONDARY_TIER_BANDS, EXPLORATORY_TIER_BANDS];
if ~isempty(setxor(sort(expectedTierBands), sort(ALL_BANDS')))
    error('Phase16:BandMismatch', ...
        ['The hard-coded band tiers do not exactly match PARAMS.BandDefinitionsHz. ' ...
         'Tiers: %s | X03 bands: %s'], strjoin(expectedTierBands, ','), strjoin(ALL_BANDS', ','));
end
tierDefs = {'Primary', PRIMARY_TIER_BANDS; 'Secondary', SECONDARY_TIER_BANDS; 'Exploratory', EXPLORATORY_TIER_BANDS};
tierOf = containers.Map();
for i = 1:numel(PRIMARY_TIER_BANDS);     tierOf(PRIMARY_TIER_BANDS{i}) = 'Primary'; end     %#ok<NOCOL>
for i = 1:numel(SECONDARY_TIER_BANDS);   tierOf(SECONDARY_TIER_BANDS{i}) = 'Secondary'; end %#ok<NOCOL>
for i = 1:numel(EXPLORATORY_TIER_BANDS); tierOf(EXPLORATORY_TIER_BANDS{i}) = 'Exploratory'; end %#ok<NOCOL>

METRICS = {'LogAbsolutePower', 'RelativePower'};
CONDITION_ORDER = {'E', 'N', 'R'};   % Baseline excluded from every model (see Section 4)
REF_LEVEL = CONDITION_ORDER{1};
CONTRASTS = {'N','E','N - E'; 'R','E','R - E'; 'R','N','R - N'};   % N-E, R-E, R-N -- X05/X11-X15 convention (see header)
FDR_ALPHA = 0.05;

ROI_NAMES = ROI_NAMES(:)';
nROI = numel(ROI_NAMES);

ELECTRODES = ChannelDefT.Label;
nElec = numel(ELECTRODES);

fprintf('Loaded X03: %d psd_data recordings, %d bands, %d reference channels.\n', numel(psd_data), nBands, nChanRef);
fprintf('Loaded X07: %d ROIs (%s), %d electrodes (locked reference order).\n\n', nROI, strjoin(ROI_NAMES, ', '), nElec);

% --- Electrode-set / ROI-set validation (see JUDGMENT CALLS item 4 for
%     why Item B's independent 62-entry canonical-order literal is not
%     re-implemented here) ---
if nElec ~= 62
    error('Phase16:ElectrodeCountMismatch', 'ChannelDefT has %d electrodes, expected 62.', nElec);
end
if ~all(ismember(upper({'P9','P10'}), upper(ELECTRODES)))
    error('Phase16:P9P10Missing', 'P9 and/or P10 missing from the electrode set loaded from X07.');
end
if numel(ROI_NAMES) ~= 5 || ~isempty(setxor(ROI_NAMES, {'Frontal','Central','Temporal','Parietal','Occipital'}))
    error('Phase16:ROINameMismatch', 'ROI_NAMES loaded from X07 is not exactly the 5 expected ROIs.');
end
allROIChannels = [ROI_DEFINITION.Frontal, ROI_DEFINITION.Central, ROI_DEFINITION.Temporal, ...
                  ROI_DEFINITION.Parietal, ROI_DEFINITION.Occipital];
if ~isempty(intersect(upper(allROIChannels), upper({'P9','P10'})))
    error('Phase16:P9P10InROI', 'P9/P10 found inside ROI_DEFINITION -- must be excluded from ROI averaging.');
end

%% ------------------------------------------------------------------
%  SECTION 2: IDENTIFY Subject x BlockNumber x Condition COMBINATIONS
%  THAT HAVE BOTH A NBACK AND A SART RECORD IN X03's psd_data (never the
%  raw Full recording -- see header CONCEPT). See JUDGMENT CALLS 1-3.
% ------------------------------------------------------------------

fprintf('--------------------------------------------------------\n');
fprintf('SECTION 2: IDENTIFYING NBACK x SART COMBINATIONS\n');
fprintf('--------------------------------------------------------\n');

integrityIssues = {};

nbackIdxAll = find(strcmp({psd_data.TaskType}, 'nback'));
sartIdxAll  = find(strcmp({psd_data.TaskType}, 'sart'));
fprintf('psd_data: %d Nback recordings, %d SART recordings (all Conditions, including Baseline).\n', ...
    numel(nbackIdxAll), numel(sartIdxAll));

nbackKeyMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
for ii = nbackIdxAll
    k = local_make_key(psd_data(ii).Subject, psd_data(ii).BlockNumber, psd_data(ii).Condition);
    if isKey(nbackKeyMap, k)
        integrityIssues{end+1} = sprintf(['Duplicate Nback psd_data entry for %s -- keeping the first, later ' ...
            'duplicate(s) ignored (never fabricated/merged).'], k); %#ok<AGROW>
    else
        nbackKeyMap(k) = ii; %#ok<NASGU>
    end
end
sartKeyMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
for jj = sartIdxAll
    k = local_make_key(psd_data(jj).Subject, psd_data(jj).BlockNumber, psd_data(jj).Condition);
    if isKey(sartKeyMap, k)
        integrityIssues{end+1} = sprintf(['Duplicate SART psd_data entry for %s -- keeping the first, later ' ...
            'duplicate(s) ignored (never fabricated/merged).'], k); %#ok<AGROW>
    else
        sartKeyMap(k) = jj; %#ok<NASGU>
    end
end

comboList = struct('Subject', {}, 'BlockNumber', {}, 'Condition', {}, 'NbackIdx', {}, 'SartIdx', {});
nbKeys = keys(nbackKeyMap);
for ki = 1:numel(nbKeys)
    k = nbKeys{ki};
    if isKey(sartKeyMap, k)
        ii = nbackKeyMap(k); jj = sartKeyMap(k);
        comboList(end+1) = struct('Subject', psd_data(ii).Subject, 'BlockNumber', psd_data(ii).BlockNumber, ...
            'Condition', psd_data(ii).Condition, 'NbackIdx', ii, 'SartIdx', jj); %#ok<AGROW>
    else
        integrityIssues{end+1} = sprintf('%s: Nback recording found, no matching SART recording -- combination SKIPPED (never fabricated).', k); %#ok<AGROW>
    end
end
saKeys = keys(sartKeyMap);
for ki = 1:numel(saKeys)
    k = saKeys{ki};
    if ~isKey(nbackKeyMap, k)
        integrityIssues{end+1} = sprintf('%s: SART recording found, no matching Nback recording -- combination SKIPPED (never fabricated).', k); %#ok<AGROW>
    end
end

fprintf('%d Subject x BlockNumber x Condition combinations have BOTH Nback and SART records (all Conditions).\n', numel(comboList));
fprintf('%d one-sided (Nback-only or SART-only) combination(s) skipped -- see integrity issues.\n\n', ...
    numel(nbKeys) + numel(saKeys) - 2*numel(comboList));

%% ------------------------------------------------------------------
%  SECTION 3: WINDOW POOLING + X04-IDENTICAL AGGREGATION MATH
% ------------------------------------------------------------------

fprintf('--------------------------------------------------------\n');
fprintf('SECTION 3: POOLING NBACK+SART VALID WINDOWS AND AGGREGATING (X04-identical math)\n');
fprintf('--------------------------------------------------------\n');

bandPowerCols = {'Subject','TaskType','BlockNumber','Condition','Channel','Band', ...
    'AbsolutePower','LogAbsolutePower','RelativePower'};
bandPowerRows = cell(0, numel(bandPowerCols));

globalCols = {'Subject','TaskType','BlockNumber','Condition','Band', ...
    'AbsolutePower','LogAbsolutePower','RelativePower'};
globalRows = cell(0, numel(globalCols));

nSkippedChannelMismatch = 0;
nSkippedZeroPooled = 0;
nAggregated = 0;

for c = 1:numel(comboList)
    subjStr = comboList(c).Subject;
    blockNum = comboList(c).BlockNumber;
    condStr = comboList(c).Condition;
    nbEntry = psd_data(comboList(c).NbackIdx);
    saEntry = psd_data(comboList(c).SartIdx);

    % --- JUDGMENT CALL 1: channel-order safety check before pooling ---
    if numel(nbEntry.ChannelLabels) ~= nChanRef || ~isequal(nbEntry.ChannelLabels(:), ChannelLabelsReference) || ...
       numel(saEntry.ChannelLabels) ~= nChanRef || ~isequal(saEntry.ChannelLabels(:), ChannelLabelsReference)
        integrityIssues{end+1} = sprintf(['Sub%s Block%g Cond%s: Nback and/or SART ChannelLabels do not exactly ' ...
            'match X03''s ChannelLabelsReference -- combination SKIPPED to avoid silent channel misalignment when ' ...
            'concatenating window matrices (see header JUDGMENT CALLS item 1).'], subjStr, blockNum, condStr); %#ok<AGROW>
        nSkippedChannelMismatch = nSkippedChannelMismatch + 1;
        continue;
    end
    if numel(nbEntry.FrequenciesHz) ~= numel(freqVecRef) || max(abs(nbEntry.FrequenciesHz(:) - freqVecRef)) > 1e-9 || ...
       numel(saEntry.FrequenciesHz) ~= numel(freqVecRef) || max(abs(saEntry.FrequenciesHz(:) - freqVecRef)) > 1e-9
        integrityIssues{end+1} = sprintf('Sub%s Block%g Cond%s: Nback and/or SART FrequenciesHz differ from the X03 reference vector.', ...
            subjStr, blockNum, condStr); %#ok<AGROW>
    end

    nbValid = nbEntry.WindowValid;
    saValid = saEntry.WindowValid;
    nPooledWindows = sum(nbValid) + sum(saValid);

    % --- JUDGMENT CALL 2: zero-pooled-window combinations are skipped;
    %     JUDGMENT CALL 3: a combination is NOT skipped merely because ONE
    %     of the two recordings alone has 0 valid windows ---
    if nPooledWindows == 0
        integrityIssues{end+1} = sprintf('Sub%s Block%g Cond%s: 0 pooled valid windows (Nback+SART) -- SKIPPED for aggregation.', ...
            subjStr, blockNum, condStr); %#ok<AGROW>
        nSkippedZeroPooled = nSkippedZeroPooled + 1;
        continue;
    end

    for b = 1:nBands
        bandName = ALL_BANDS{b};
        pooledAbs = [nbEntry.AbsoluteBandPower.(bandName)(nbValid, :); saEntry.AbsoluteBandPower.(bandName)(saValid, :)];
        pooledRel = [nbEntry.RelativeBandPower.(bandName)(nbValid, :); saEntry.RelativeBandPower.(bandName)(saValid, :)];

        absMean = mean(pooledAbs, 1);   % arithmetic mean across pooled windows -- X04-identical
        relMean = mean(pooledRel, 1);
        absLog  = log(absMean);         % log(mean(...)), NEVER mean(log(...)) -- X04-identical
        if any(absMean <= 0)
            integrityIssues{end+1} = sprintf(['Sub%s Block%g Cond%s Band %s: non-positive pooled AbsolutePower found ' ...
                'in %d channel(s) -- LogAbsolutePower set to NaN there, not silently patched.'], ...
                subjStr, blockNum, condStr, bandName, sum(absMean <= 0)); %#ok<AGROW>
            absLog(absMean <= 0) = NaN;
        end

        for ch = 1:nChanRef
            bandPowerRows(end+1, :) = {subjStr, 'FullReconstructed', blockNum, condStr, ChannelLabelsReference{ch}, bandName, ...
                absMean(ch), absLog(ch), relMean(ch)}; %#ok<AGROW>
        end

        % whole-scalp: arithmetic mean across channels of the per-channel
        % means -- X04 Structure C-identical
        globalAbs = mean(absMean);
        globalRel = mean(relMean);
        globalLog = log(globalAbs);
        if globalAbs <= 0
            integrityIssues{end+1} = sprintf('Sub%s Block%g Cond%s Band %s: non-positive whole-scalp pooled AbsolutePower -- LogAbsolutePower set to NaN.', ...
                subjStr, blockNum, condStr, bandName); %#ok<AGROW>
            globalLog = NaN;
        end
        globalRows(end+1, :) = {subjStr, 'FullReconstructed', blockNum, condStr, bandName, globalAbs, globalLog, globalRel}; %#ok<AGROW>
    end

    nAggregated = nAggregated + 1;
end

BandPowerT_FullRecon = local_build_table(bandPowerRows, [1 2 4 5 6], [3 7 8 9], bandPowerCols);
GlobalBandPowerT_FullRecon = local_build_table(globalRows, [1 2 4 5], [3 6 7 8], globalCols);

fprintf('Aggregation complete: %d / %d combinations aggregated (%d skipped: channel mismatch, %d skipped: 0 pooled windows).\n', ...
    nAggregated, numel(comboList), nSkippedChannelMismatch, nSkippedZeroPooled);
fprintf('BandPowerT_FullRecon: %d rows. GlobalBandPowerT_FullRecon: %d rows.\n\n', ...
    height(BandPowerT_FullRecon), height(GlobalBandPowerT_FullRecon));

save(fullfile(STATS_DIR, 'full_reconstructed_bandpower.mat'), 'BandPowerT_FullRecon', 'GlobalBandPowerT_FullRecon');

%% ------------------------------------------------------------------
%  SECTION 4: FILTER TO Condition IN {E,N,R} (Baseline excluded before
%  any modeling, exactly like X05/X08/X09/X11's own policy), THEN SHARED
%  DESIGN VERIFICATION (the same recordings feed all three levels below).
% ------------------------------------------------------------------

fprintf('--------------------------------------------------------\n');
fprintf('SECTION 4: BASELINE EXCLUSION + DESIGN VERIFICATION\n');
fprintf('--------------------------------------------------------\n');

nBaselineExcludedChan = sum(strcmp(BandPowerT_FullRecon.Condition, 'Baseline'));
BandPowerT = BandPowerT_FullRecon(ismember(BandPowerT_FullRecon.Condition, CONDITION_ORDER), :);
nBaselineExcludedGlobal = sum(strcmp(GlobalBandPowerT_FullRecon.Condition, 'Baseline'));
GlobalBandPowerT = GlobalBandPowerT_FullRecon(ismember(GlobalBandPowerT_FullRecon.Condition, CONDITION_ORDER), :);

fprintf('Channel-level: %d Baseline rows excluded, %d E/N/R rows retained.\n', nBaselineExcludedChan, height(BandPowerT));
fprintf('Global level : %d Baseline rows excluded, %d E/N/R rows retained.\n', nBaselineExcludedGlobal, height(GlobalBandPowerT));

% ---- NEW: apply the shared noisy-subject exclusion list (see
% EXCLUSION_CONFIG_cognitive.m), only when RUN_MODE=sensitivity_excl_noisy.
% Never applied to the full_sample run. Applied ONCE here, to BOTH
% BandPowerT (channel-level) and GlobalBandPowerT (whole-scalp) -- the
% two tables that feed ALL THREE downstream sections (Global uses
% GlobalBandPowerT directly; ROI aggregation and Electrode both use
% BandPowerT) -- so Global/ROI/Electrode all see the SAME excluded-subject
% set, applied before any of the three sections' modeling begins. ----
nExcludedRowsSensitivityChan = 0;
nExcludedRowsSensitivityGlobal = 0;
if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
    subjNumForExclChan = local_subject_numeric(BandPowerT.Subject);
    excludeMaskChan = ismember(subjNumForExclChan, EXCLUDE_NOISY_SUBJECTS);
    nExcludedRowsSensitivityChan = sum(excludeMaskChan);
    BandPowerT = BandPowerT(~excludeMaskChan, :);

    subjNumForExclGlobal = local_subject_numeric(GlobalBandPowerT.Subject);
    excludeMaskGlobal = ismember(subjNumForExclGlobal, EXCLUDE_NOISY_SUBJECTS);
    nExcludedRowsSensitivityGlobal = sum(excludeMaskGlobal);
    GlobalBandPowerT = GlobalBandPowerT(~excludeMaskGlobal, :);

    fprintf('RUN_MODE=%s -- excluded %d channel-level row(s) and %d Global-level row(s) belonging to subject(s) {%s}; %d Global-level rows remain.\n', ...
        RUN_MODE_LABEL, nExcludedRowsSensitivityChan, nExcludedRowsSensitivityGlobal, ...
        strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','), height(GlobalBandPowerT));
end
% ---- END NEW ----

if isempty(GlobalBandPowerT)
    error('Phase16:NoData', 'No E/N/R rows remain after Baseline exclusion -- cannot proceed.');
end

BandPowerT.ChannelUC = upper(BandPowerT.Channel);

comboT = unique(GlobalBandPowerT(:, {'Subject','BlockNumber','Condition'}), 'rows');
nRecordings = height(comboT);
fprintf('%d unique Subject x BlockNumber x Condition FullReconstructed recordings found (E/N/R only).\n', nRecordings);

subjNumeric = local_subject_numeric(comboT.Subject);
uniqSubj = unique(subjNumeric);
nSubjects = numel(uniqSubj);

subjPairs = cell(nSubjects, 1);
pairTypeCount = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:nSubjects
    rowsI = subjNumeric == uniqSubj(i);
    condsI = cellstr(comboT.Condition(rowsI));
    if numel(condsI) ~= 2
        integrityIssues{end+1} = sprintf('Subject %g has %d FullReconstructed recordings (expected 2).', uniqSubj(i), numel(condsI)); %#ok<AGROW>
    elseif strcmp(condsI{1}, condsI{2})
        integrityIssues{end+1} = sprintf('Subject %g has the SAME condition (%s) in both FullReconstructed recordings -- not a valid pair.', uniqSubj(i), condsI{1}); %#ok<AGROW>
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
fprintf('Condition-pair coverage across %d subjects:\n', nSubjects);
for k = 1:numel(pairKeys)
    fprintf('  %s : %d subjects\n', strrep(pairKeys{k}, '_', ' & '), pairTypeCount(pairKeys{k}));
end

allSubjectsComplete3 = false;
condsPerSubjCounts = cellfun(@numel, subjPairs);
if all(condsPerSubjCounts >= 3)
    allSubjectsComplete3 = true;
    integrityIssues{end+1} = ['Every subject appears to have 3+ FullReconstructed recordings/conditions -- the ' ...
        'incomplete-design justification for using an LMM instead of RM-ANOVA no longer applies. STOP and re-examine.'];
end

condList = cellstr(comboT.Condition);
blockList = comboT.BlockNumber;
confoundFlag = false;
for c = 1:numel(CONDITION_ORDER)
    thisCond = CONDITION_ORDER{c};
    blocksForCond = blockList(strcmp(condList, thisCond));
    uBlocks = unique(blocksForCond);
    fprintf('  Condition %s: BlockNumber counts = %s\n', thisCond, mat2str(arrayfun(@(u) sum(blocksForCond == u), uBlocks)'));
    if numel(uBlocks) < 2
        confoundFlag = true;
        integrityIssues{end+1} = sprintf(['Condition %s occurs at only one BlockNumber value -- Condition and block ' ...
            'order are FULLY confounded for this condition.'], thisCond); %#ok<AGROW>
    end
end

if allSubjectsComplete3
    error('Phase16:DesignAssumptionViolated', ...
        'Design verification found a complete 3-condition-per-subject structure. Stopping before fitting -- see integrityIssues.');
end

if ~isempty(integrityIssues)
    fprintf(2, '\n*** DESIGN VERIFICATION / AGGREGATION FLAGGED %d ISSUE(S) SO FAR -- see report. ***\n\n', numel(integrityIssues));
end
fprintf('\n');

%% ====================================================================
%  SECTION 5: GLOBAL-LEVEL ANALYSIS (X05-identical model/FDR/pairwise)
% ====================================================================

fprintf('--------------------------------------------------------\n');
fprintf('SECTION 5: GLOBAL-LEVEL LMM/FDR/PAIRWISE (%d bands x %d metrics = %d fits)\n', nBands, numel(METRICS), nBands*numel(METRICS));
fprintf('--------------------------------------------------------\n');

globalOmnibusRows = cell(0, 6);
lmeStoreGlobal = containers.Map();   % key: 'Metric|Band'
convergenceWarningsGlobal = {};

for mI = 1:numel(METRICS)
    metricName = METRICS{mI};
    for b = 1:nBands
        bandName = ALL_BANDS{b};
        subT = GlobalBandPowerT(strcmp(GlobalBandPowerT.Band, bandName), :);
        fr = local_fit_condition_omnibus(subT, metricName, CONDITION_ORDER);
        if ~fr.ok
            integrityIssues{end+1} = sprintf('Global: fitlme FAILED for %s / %s: %s', metricName, bandName, fr.errMsg); %#ok<AGROW>
            continue;
        end
        if ~isempty(fr.warnMsg)
            convergenceWarningsGlobal{end+1} = sprintf('%s / %s: %s', metricName, bandName, fr.warnMsg); %#ok<AGROW>
        end
        lmeStoreGlobal(sprintf('%s|%s', metricName, bandName)) = fr.lme;
        globalOmnibusRows(end+1, :) = {metricName, bandName, fr.F, fr.DF1, fr.DF2, fr.p}; %#ok<AGROW>
    end
end
GlobalOmnibusT = local_build_table(globalOmnibusRows, [1 2], [3 4 5 6], {'Metric','Band','FStat','DF1','DF2','pValue'});

GlobalOmnibusT.Tier = cellfun(@(b) tierOf(b), GlobalOmnibusT.Band, 'UniformOutput', false);
GlobalOmnibusT.pFDR = nan(height(GlobalOmnibusT), 1);
GlobalOmnibusT.SignificantFDR = false(height(GlobalOmnibusT), 1);
GlobalOmnibusT.FamilySize = nan(height(GlobalOmnibusT), 1);

for t = 1:size(tierDefs, 1)
    tierBands = tierDefs{t, 2};
    for mI = 1:numel(METRICS)
        metricName = METRICS{mI};
        idx = strcmp(GlobalOmnibusT.Metric, metricName) & ismember(GlobalOmnibusT.Band, tierBands);
        if sum(idx) == 0; continue; end
        GlobalOmnibusT.FamilySize(idx) = sum(idx);
        pvals = GlobalOmnibusT.pValue(idx);
        [adjP, reject] = local_bh_fdr_or_raw(pvals, FDR_ALPHA);
        GlobalOmnibusT.pFDR(idx) = adjP;
        GlobalOmnibusT.SignificantFDR(idx) = reject;
    end
end

fprintf('Global-level omnibus tests fit: %d / %d expected.\n', height(GlobalOmnibusT), nBands*numel(METRICS));

GlobalPairwiseRows = cell(0, 11);
sigMaskGlobal = GlobalOmnibusT.SignificantFDR;
for row = find(sigMaskGlobal)'
    metricName = GlobalOmnibusT.Metric{row};
    bandName = GlobalOmnibusT.Band{row};
    key = sprintf('%s|%s', metricName, bandName);
    if ~isKey(lmeStoreGlobal, key); continue; end
    lme = lmeStoreGlobal(key);
    rows3 = local_run_pairwise(lme, REF_LEVEL, CONTRASTS, FDR_ALPHA);
    for c = 1:3
        GlobalPairwiseRows(end+1, :) = [{metricName, bandName}, rows3(c, :)]; %#ok<AGROW>
    end
end
GlobalPairwiseT = local_build_table(GlobalPairwiseRows, [1 2 3], [4 5 6 7 8 9 10 11], ...
    {'Metric','Band','Contrast','Estimate','SE','DF','tStat','pValue','pHolm','CI_Low','CI_High'});

fprintf('Global-level pairwise contrasts computed: %d rows (%d FDR-surviving Band x Metric x 3 contrasts).\n\n', ...
    height(GlobalPairwiseT), sum(sigMaskGlobal));

% --- Residual diagnostics: Primary+Secondary tier bands, both metrics --
%     mirrors X05's own diagnostic scope exactly (see JUDGMENT CALLS 5) ---
diagBandsGlobal = unique([PRIMARY_TIER_BANDS, SECONDARY_TIER_BANDS]);
nDiagGlobal = 0;
for mI = 1:numel(METRICS)
    metricName = METRICS{mI};
    for bi = 1:numel(diagBandsGlobal)
        bandName = diagBandsGlobal{bi};
        key = sprintf('%s|%s', metricName, bandName);
        if ~isKey(lmeStoreGlobal, key); continue; end
        lme = lmeStoreGlobal(key);
        try
            res = residuals(lme, 'ResidualType', 'Pearson');
            fig = figure('Visible', 'off', 'Position', [100 100 900 400]);
            subplot(1,2,1);
            histogram(res, 15);
            title(sprintf('Residuals: Global %s (%s)', bandName, metricName), 'Interpreter', 'none');
            xlabel('Pearson residual'); ylabel('Count');
            subplot(1,2,2);
            qqplot(res);
            title('Normal Q-Q');
            outPng = fullfile(STATS_DIR, sprintf('diagnostic_residuals_Global_%s_%s.png', bandName, metricName));
            saveas(fig, outPng);
            close(fig);
            nDiagGlobal = nDiagGlobal + 1;
        catch ME
            integrityIssues{end+1} = sprintf('Global residual diagnostic failed for %s / %s: %s', metricName, bandName, ME.message); %#ok<AGROW>
        end
    end
end
fprintf('Global-level residual diagnostic plots saved: %d.\n\n', nDiagGlobal);

save(fullfile(STATS_DIR, 'global_omnibus_results.mat'), 'GlobalOmnibusT');
save(fullfile(STATS_DIR, 'global_pairwise_contrasts.mat'), 'GlobalPairwiseT');
writetable(GlobalOmnibusT, fullfile(STATS_DIR, 'global_omnibus_results.xlsx'));
if ~isempty(GlobalPairwiseT)
    writetable(GlobalPairwiseT, fullfile(STATS_DIR, 'global_pairwise_contrasts.xlsx'));
end

%% ------------------------------------------------------------------
%  SECTION 6a: ROI AGGREGATION (X08-identical rule: arithmetic mean of
%  channel-level AbsolutePower/RelativePower across ROI member channels;
%  LogAbsolutePower = log(mean AbsolutePower). ROI membership loaded
%  READ-ONLY from X07, never redefined.)
% ------------------------------------------------------------------

fprintf('--------------------------------------------------------\n');
fprintf('SECTION 6a: ROI AGGREGATION\n');
fprintf('--------------------------------------------------------\n');

bandPowerChannels = unique(BandPowerT.Channel);
missingFromBandPower = ChannelLabelsReference(~ismember(upper(ChannelLabelsReference), upper(bandPowerChannels)));
if ~isempty(missingFromBandPower)
    error('Phase16:MissingChannels', 'Reference channel(s) not found in BandPowerT: %s', strjoin(missingFromBandPower, ', '));
end
missingROIChan = allROIChannels(~ismember(upper(allROIChannels), upper(bandPowerChannels)));
if ~isempty(missingROIChan)
    error('Phase16:MissingROIChannels', 'ROI-member channel(s) not found in BandPowerT: %s', strjoin(missingROIChan, ', '));
end

roiRows = cell(0, 7);
roiCols = {'Subject','BlockNumber','Condition','ROI','Band','AbsolutePower','RelativePower'};
nonPositiveROIAbs = {};

for i = 1:height(comboT)
    subjStr = comboT.Subject{i}; blockNum = comboT.BlockNumber(i); condStr = comboT.Condition{i};
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
    integrityIssues{end+1} = sprintf('%d ROI x Band x recording combination(s) had non-positive ROI AbsolutePower -- LogAbsolutePower set to NaN there.', numel(nonPositiveROIAbs)); %#ok<AGROW>
end

fprintf('ROI aggregation complete: %d rows (%d recordings x %d ROIs x %d bands).\n\n', ...
    height(RoiBandPowerT), height(comboT), nROI, nBands);

%% ------------------------------------------------------------------
%  SECTION 6b: ROI-LEVEL LMM/FDR/PAIRWISE (X08-identical FDR family
%  scheme: BH-FDR per Band x Metric family, across the 5 ROIs)
% ------------------------------------------------------------------

fprintf('--------------------------------------------------------\n');
fprintf('SECTION 6b: ROI-LEVEL LMM/FDR/PAIRWISE (%d ROIs x %d bands x %d metrics = %d fits)\n', ...
    nROI, nBands, numel(METRICS), nROI*nBands*numel(METRICS));
fprintf('--------------------------------------------------------\n');

roiOmnibusRows = cell(0, 7);
lmeMapROI = containers.Map();   % key: 'ROI|Metric|Band'

for r = 1:nROI
    roiName = ROI_NAMES{r};
    for mI = 1:numel(METRICS)
        metricName = METRICS{mI};
        for b = 1:nBands
            bandName = ALL_BANDS{b};
            subT = RoiBandPowerT(strcmp(RoiBandPowerT.ROI, roiName) & strcmp(RoiBandPowerT.Band, bandName), :);
            fr = local_fit_condition_omnibus(subT, metricName, CONDITION_ORDER);
            if ~fr.ok
                integrityIssues{end+1} = sprintf('ROI=%s: fitlme FAILED for %s / %s: %s', roiName, metricName, bandName, fr.errMsg); %#ok<AGROW>
                continue;
            end
            lmeMapROI(sprintf('%s|%s|%s', roiName, metricName, bandName)) = fr.lme;
            roiOmnibusRows(end+1, :) = {roiName, metricName, bandName, fr.F, fr.DF1, fr.DF2, fr.p}; %#ok<AGROW>
        end
    end
end
RoiOmnibusT = local_build_table(roiOmnibusRows, [1 2 3], [4 5 6 7], {'ROI','Metric','Band','FStat','DF1','DF2','pValue'});

RoiOmnibusT.Tier = cellfun(@(b) tierOf(b), RoiOmnibusT.Band, 'UniformOutput', false);
RoiOmnibusT.pFDR = nan(height(RoiOmnibusT), 1);
RoiOmnibusT.SignificantFDR = false(height(RoiOmnibusT), 1);
RoiOmnibusT.FamilySize = nan(height(RoiOmnibusT), 1);

for mI = 1:numel(METRICS)
    metricName = METRICS{mI};
    for b = 1:nBands
        bandName = ALL_BANDS{b};
        idx = strcmp(RoiOmnibusT.Metric, metricName) & strcmp(RoiOmnibusT.Band, bandName);
        nInFamily = sum(idx);
        if nInFamily == 0; continue; end
        RoiOmnibusT.FamilySize(idx) = nInFamily;
        if nInFamily ~= nROI
            integrityIssues{end+1} = sprintf('ROI FDR family %s / %s has %d members, expected %d.', metricName, bandName, nInFamily, nROI); %#ok<AGROW>
        end
        pvals = RoiOmnibusT.pValue(idx);
        [adjP, reject] = local_fdr_bh(pvals, FDR_ALPHA);
        RoiOmnibusT.pFDR(idx) = adjP;
        RoiOmnibusT.SignificantFDR(idx) = reject;
    end
end

fprintf('ROI-level omnibus tests fit: %d / %d expected.\n', height(RoiOmnibusT), nROI*nBands*numel(METRICS));

RoiPairwiseRows = cell(0, 12);
sigMaskROI = RoiOmnibusT.SignificantFDR;
for row = find(sigMaskROI)'
    roiName = RoiOmnibusT.ROI{row}; metricName = RoiOmnibusT.Metric{row}; bandName = RoiOmnibusT.Band{row};
    key = sprintf('%s|%s|%s', roiName, metricName, bandName);
    if ~isKey(lmeMapROI, key); continue; end
    lme = lmeMapROI(key);
    rows3 = local_run_pairwise(lme, REF_LEVEL, CONTRASTS, FDR_ALPHA);
    for c = 1:3
        RoiPairwiseRows(end+1, :) = [{roiName, metricName, bandName}, rows3(c, :)]; %#ok<AGROW>
    end
end
RoiPairwiseT = local_build_table(RoiPairwiseRows, [1 2 3 4], [5 6 7 8 9 10 11 12], ...
    {'ROI','Metric','Band','Contrast','Estimate','SE','DF','tStat','pValue','pHolm','CI_Low','CI_High'});

fprintf('ROI-level pairwise contrasts computed: %d rows (%d FDR-surviving ROI x Band x Metric x 3 contrasts).\n\n', ...
    height(RoiPairwiseT), sum(sigMaskROI));

% --- Residual diagnostics: Primary+Secondary bands, all 5 ROIs, both
%     metrics -- mirrors X08's own diagnostic scope exactly ---
diagBandsROI = [PRIMARY_TIER_BANDS, SECONDARY_TIER_BANDS];
nDiagROI = 0;
for r = 1:nROI
    roiName = ROI_NAMES{r};
    for mI = 1:numel(METRICS)
        metricName = METRICS{mI};
        for bi = 1:numel(diagBandsROI)
            bandName = diagBandsROI{bi};
            key = sprintf('%s|%s|%s', roiName, metricName, bandName);
            if ~isKey(lmeMapROI, key); continue; end
            lme = lmeMapROI(key);
            try
                res = residuals(lme, 'ResidualType', 'Pearson');
                fig = figure('Visible', 'off', 'Position', [100 100 900 400]);
                subplot(1,2,1);
                histogram(res, 15);
                title(sprintf('Residuals: %s / %s (%s)', roiName, bandName, metricName), 'Interpreter', 'none');
                xlabel('Pearson residual'); ylabel('Count');
                subplot(1,2,2);
                qqplot(res);
                title('Normal Q-Q');
                outPng = fullfile(STATS_DIR, sprintf('diagnostic_residuals_ROI_%s_%s_%s.png', roiName, bandName, metricName));
                saveas(fig, outPng);
                close(fig);
                nDiagROI = nDiagROI + 1;
            catch ME
                integrityIssues{end+1} = sprintf('ROI residual diagnostic failed for %s / %s / %s: %s', roiName, metricName, bandName, ME.message); %#ok<AGROW>
            end
        end
    end
end
fprintf('ROI-level residual diagnostic plots saved: %d.\n\n', nDiagROI);

save(fullfile(STATS_DIR, 'roi_omnibus_results.mat'), 'RoiOmnibusT');
save(fullfile(STATS_DIR, 'roi_pairwise_contrasts.mat'), 'RoiPairwiseT');
writetable(RoiOmnibusT, fullfile(STATS_DIR, 'roi_omnibus_results.xlsx'));
writetable(RoiOmnibusT, fullfile(STATS_DIR, 'roi_omnibus_results.csv'));
if ~isempty(RoiPairwiseT)
    writetable(RoiPairwiseT, fullfile(STATS_DIR, 'roi_pairwise_contrasts.xlsx'));
    writetable(RoiPairwiseT, fullfile(STATS_DIR, 'roi_pairwise_contrasts.csv'));
end

%% ------------------------------------------------------------------
%  SECTION 7: ELECTRODE-LEVEL LMM/FDR/PAIRWISE (X09-identical FDR family
%  scheme: BH-FDR per Band x Metric family, across the 62 electrodes; NO
%  aggregation of any kind; P9/P10 included; NO residual diagnostics,
%  mirroring X09 which also produces none.)
% ------------------------------------------------------------------

fprintf('--------------------------------------------------------\n');
fprintf('SECTION 7: ELECTRODE-LEVEL LMM/FDR/PAIRWISE (%d electrodes x %d bands x %d metrics = %d fits)\n', ...
    nElec, nBands, numel(METRICS), nElec*nBands*numel(METRICS));
fprintf('--------------------------------------------------------\n');

elecOmnibusRows = cell(0, 7);
lmeMapElec = containers.Map();
nOKElec = 0; nFailedElec = 0;

for e = 1:nElec
    elecName = ELECTRODES{e};
    for mI = 1:numel(METRICS)
        metricName = METRICS{mI};
        for b = 1:nBands
            bandName = ALL_BANDS{b};
            subT = BandPowerT(strcmp(BandPowerT.ChannelUC, upper(elecName)) & strcmp(BandPowerT.Band, bandName), :);
            fr = local_fit_condition_omnibus(subT, metricName, CONDITION_ORDER);
            if ~fr.ok
                nFailedElec = nFailedElec + 1;
                integrityIssues{end+1} = sprintf('Electrode=%s: fitlme FAILED for %s / %s: %s', elecName, metricName, bandName, fr.errMsg); %#ok<AGROW>
                continue;
            end
            nOKElec = nOKElec + 1;
            lmeMapElec(sprintf('%s|%s|%s', elecName, metricName, bandName)) = fr.lme;
            elecOmnibusRows(end+1, :) = {elecName, metricName, bandName, fr.F, fr.DF1, fr.DF2, fr.p}; %#ok<AGROW>
        end
    end
    if mod(e, 10) == 0 || e == nElec
        fprintf('  ... %d/%d electrodes processed (%d OK, %d failed so far)\n', e, nElec, nOKElec, nFailedElec);
    end
end
ElectrodeOmnibusT = local_build_table(elecOmnibusRows, [1 2 3], [4 5 6 7], {'Electrode','Metric','Band','FStat','DF1','DF2','pValue'});
ElectrodeOmnibusT.Tier = cellfun(@(b) tierOf(b), ElectrodeOmnibusT.Band, 'UniformOutput', false);
ElectrodeOmnibusT.pFDR = nan(height(ElectrodeOmnibusT), 1);
ElectrodeOmnibusT.SignificantFDR = false(height(ElectrodeOmnibusT), 1);
ElectrodeOmnibusT.FamilySize = nan(height(ElectrodeOmnibusT), 1);

for mI = 1:numel(METRICS)
    metricName = METRICS{mI};
    for b = 1:nBands
        bandName = ALL_BANDS{b};
        idx = strcmp(ElectrodeOmnibusT.Metric, metricName) & strcmp(ElectrodeOmnibusT.Band, bandName);
        nInFamily = sum(idx);
        if nInFamily == 0; continue; end
        ElectrodeOmnibusT.FamilySize(idx) = nInFamily;
        if nInFamily ~= nElec
            integrityIssues{end+1} = sprintf('Electrode FDR family %s / %s has %d members, expected %d.', metricName, bandName, nInFamily, nElec); %#ok<AGROW>
        end
        pvals = ElectrodeOmnibusT.pValue(idx);
        [adjP, reject] = local_fdr_bh(pvals, FDR_ALPHA);
        ElectrodeOmnibusT.pFDR(idx) = adjP;
        ElectrodeOmnibusT.SignificantFDR(idx) = reject;
    end
end

fprintf('Electrode-level omnibus tests fit: %d OK, %d failed, out of %d expected.\n', nOKElec, nFailedElec, nElec*nBands*numel(METRICS));

ElectrodePairwiseRows = cell(0, 12);
sigMaskElec = ElectrodeOmnibusT.SignificantFDR;
for row = find(sigMaskElec)'
    elecName = ElectrodeOmnibusT.Electrode{row}; metricName = ElectrodeOmnibusT.Metric{row}; bandName = ElectrodeOmnibusT.Band{row};
    key = sprintf('%s|%s|%s', elecName, metricName, bandName);
    if ~isKey(lmeMapElec, key); continue; end
    lme = lmeMapElec(key);
    rows3 = local_run_pairwise(lme, REF_LEVEL, CONTRASTS, FDR_ALPHA);
    for c = 1:3
        ElectrodePairwiseRows(end+1, :) = [{elecName, metricName, bandName}, rows3(c, :)]; %#ok<AGROW>
    end
end
ElectrodePairwiseT = local_build_table(ElectrodePairwiseRows, [1 2 3 4], [5 6 7 8 9 10 11 12], ...
    {'Electrode','Metric','Band','Contrast','Estimate','SE','DF','tStat','pValue','pHolm','CI_Low','CI_High'});

fprintf('Electrode-level pairwise contrasts computed: %d rows (%d FDR-surviving Electrode x Band x Metric x 3 contrasts).\n\n', ...
    height(ElectrodePairwiseT), sum(sigMaskElec));

save(fullfile(STATS_DIR, 'electrode_omnibus_results.mat'), 'ElectrodeOmnibusT');
save(fullfile(STATS_DIR, 'electrode_pairwise_contrasts.mat'), 'ElectrodePairwiseT');
writetable(ElectrodeOmnibusT, fullfile(STATS_DIR, 'electrode_omnibus_results.xlsx'));
writetable(ElectrodeOmnibusT, fullfile(STATS_DIR, 'electrode_omnibus_results.csv'));
if ~isempty(ElectrodePairwiseT)
    writetable(ElectrodePairwiseT, fullfile(STATS_DIR, 'electrode_pairwise_contrasts.xlsx'));
    writetable(ElectrodePairwiseT, fullfile(STATS_DIR, 'electrode_pairwise_contrasts.csv'));
end

%% ====================================================================
%  SECTION 8: PARAMETERS AND OUTPUT
% ====================================================================

PARAMS16 = struct();
PARAMS16.ScriptName = 'X16_full_reconstructed_cognitive.m';
PARAMS16.GeneratedOn = datestr(now); %#ok<TNOW1,DATST>
PARAMS16.Purpose = ['Purely additive, exploratory analysis. Builds a synthetic ' ...
    '''FullReconstructed'' pseudo-Task from the UNION of each Subject x BlockNumber x Condition''s own ' ...
    'already-existing Nback and SART valid analysis windows (X03 psd_data), aggregated with X04''s exact ' ...
    'arithmetic-mean / log(mean(.)) formula. Does NOT replace or modify the existing Full-task analysis ' ...
    '(X05/X06/X08/X09, TaskType=''full'').'];
PARAMS16.WindowPoolingRule = ['Vertical concatenation of Nback''s and SART''s own WindowValid-filtered ' ...
    'AbsoluteBandPower/RelativeBandPower window x channel matrices, per band, for every Subject x BlockNumber x ' ...
    'Condition with BOTH a Nback and a SART record.'];
PARAMS16.AggregationMath = ['Identical to X04_aggregate_cognitive_subject_level.m: arithmetic mean across pooled ' ...
    'windows; LogAbsolutePower = log(mean(AbsolutePower)), never mean(log(.)); non-positive AbsolutePower -> NaN, ' ...
    'logged, never silently patched; whole-scalp Global = arithmetic mean across channels of the per-channel means.'];
PARAMS16.Model = 'LinearMixedEffects_RandomIntercept';
PARAMS16.FixedEffectsPrimary = 'Condition (3-level: E, N, R; reference = E)';
PARAMS16.RandomEffects = '(1|Subject)';
PARAMS16.FitMethod = 'REML';
PARAMS16.DFMethod = 'Satterthwaite';
PARAMS16.Metrics = METRICS;
PARAMS16.BandTiers = struct('Primary', {PRIMARY_TIER_BANDS}, 'Secondary', {SECONDARY_TIER_BANDS}, 'Exploratory', {EXPLORATORY_TIER_BANDS});
PARAMS16.FDRAlpha = FDR_ALPHA;
PARAMS16.FDRFamilyScheme = struct( ...
    'Global',    'BH-FDR per (Tier x Metric), across bands in that tier -- X05-identical (Primary n=1 reported as-is, never treated as pre-corrected).', ...
    'ROI',       'BH-FDR per (Band x Metric), across the 5 ROIs -- X08-identical.', ...
    'Electrode', 'BH-FDR per (Band x Metric), across the 62 electrodes -- X09-identical.');
PARAMS16.PairwiseContrastMethod = ['Holm-Bonferroni across N-E, R-E, R-N (X05/X11-X15 convention, used at ALL ' ...
    'THREE levels here per explicit instruction -- NOT X08/X09''s own E-N/E-R/N-R convention), computed only for ' ...
    'Level x Unit x Band x Metric combinations whose own FDR family test survives BH-FDR.'];
PARAMS16.ROIAggregationRule = ['Arithmetic mean of channel-level AbsolutePower/RelativePower across ROI member ' ...
    'channels; LogAbsolutePower = log(mean AbsolutePower) -- X08-identical.'];
PARAMS16.ROINames = ROI_NAMES;
PARAMS16.UnassignedForROI = {UNASSIGNED_FOR_ROI};
PARAMS16.Electrodes = ELECTRODES;
PARAMS16.BaselineExclusionNote = sprintf(['%d channel-level / %d Global-level Baseline row(s) excluded before ' ...
    'any modeling; Baseline is aggregated (Section 3) but never enters any model in this script.'], ...
    nBaselineExcludedChan, nBaselineExcludedGlobal);
PARAMS16.RunMode = RUN_MODE;
if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
    PARAMS16.ExcludedSubjects = EXCLUDE_NOISY_SUBJECTS;
    PARAMS16.ExcludedSubjectsRationale = EXCLUSION_RATIONALE_cognitive;
    PARAMS16.RowsExcludedForSensitivity = struct( ...
        'ChannelLevel', nExcludedRowsSensitivityChan, ...
        'GlobalLevel', nExcludedRowsSensitivityGlobal);
else
    PARAMS16.ExcludedSubjects = [];
end
PARAMS16.JudgmentCalls = { ...
    ['Channel-order safety check before pooling: a combination is skipped (not just logged) if either the Nback ' ...
     'or SART ChannelLabels does not exactly equal X03''s ChannelLabelsReference -- stricter than X04''s own ' ...
     'log-only policy, because this script (unlike X04) concatenates two recordings'' window matrices ' ...
     'column-for-column. See header item 1.'], ...
    ['Zero-pooled-window combinations (union of Nback+SART valid windows = 0) are skipped for aggregation and ' ...
     'logged, mirroring X04''s own "skip aggregation, still record the attempt" policy. See header item 2.'], ...
    ['A combination where only ONE of the two Task recordings has 0 valid windows is NOT skipped -- the union ' ...
     'proceeds using whichever recording has valid windows. See header item 3.'], ...
    ['Electrode-level Item B (X09''s own independent 62-entry canonical-order literal, cross-checked a second ' ...
     'time against X07) is not re-implemented here; this script validates only electrode count (62) and P9/P10 ' ...
     'presence, relying on X07 as the single locked source. See header item 4.'], ...
    ['Residual-diagnostic scope mirrors X05 (Global) and X08 (ROI) exactly (Primary+Secondary tier bands); ' ...
     'Electrode level produces none, mirroring X09. This deliberately does not follow X11''s own Global-only ' ...
     'diagnostic scope. See header item 5.'], ...
    ['full_reconstructed_bandpower.mat (the X04 Structure B/C analog for this pseudo-Task) is saved as an extra ' ...
     'output beyond the literal file list in the task description, for independent inspection/reproducibility -- ' ...
     'see header item 6.']};
PARAMS16.SourcePsdFile = psdMatFile;
PARAMS16.SourceRoiDefFile = roiDefFile;
PARAMS16.SourceRoiDefDecision = PARAMS5D.ROIDecision;
PARAMS16.DesignVerification = struct('nRecordings', nRecordings, 'nSubjects', nSubjects, ...
    'PairTypeCounts', {pairTypeCount.keys}, 'ConfoundFlag', confoundFlag, ...
    'nCombosBothTasks', numel(comboList), 'nAggregated', nAggregated, ...
    'nSkippedChannelMismatch', nSkippedChannelMismatch, 'nSkippedZeroPooled', nSkippedZeroPooled, ...
    'IntegrityIssues', {integrityIssues});
PARAMS16.RunDateTime = datestr(now); %#ok<TNOW1,DATST>

save(fullfile(STATS_DIR, 'parameters.mat'), 'PARAMS16');

%% ====================================================================
%  SECTION 9: REPORT
% ====================================================================

reportLines = {};
reportLines{end+1} = repmat('=', 1, 60);
reportLines{end+1} = 'PHASE 16 (X16) -- FULL RECONSTRUCTED EXPLORATORY ADDITIVE ANALYSIS REPORT';
reportLines{end+1} = repmat('=', 1, 60);
reportLines{end+1} = '';
reportLines{end+1} = 'PURELY ADDITIVE: does not replace or modify the existing Full-task analysis';
reportLines{end+1} = '(X05/X06/X08/X09, TaskType=''full''). See script header for full methodology.';
reportLines{end+1} = '';
reportLines{end+1} = sprintf('Nback x SART combinations with BOTH records: %d', numel(comboList));
reportLines{end+1} = sprintf('  Aggregated successfully: %d', nAggregated);
reportLines{end+1} = sprintf('  Skipped (channel-order mismatch): %d', nSkippedChannelMismatch);
reportLines{end+1} = sprintf('  Skipped (0 pooled valid windows): %d', nSkippedZeroPooled);
reportLines{end+1} = sprintf('Baseline rows excluded before modeling: %d channel-level, %d Global-level', nBaselineExcludedChan, nBaselineExcludedGlobal);
reportLines{end+1} = sprintf('FullReconstructed recordings (E/N/R only) verified: %d | Subjects: %d', nRecordings, nSubjects);
reportLines{end+1} = sprintf('Condition-pair coverage: %s', strjoin(cellfun(@(k) sprintf('%s=%d', strrep(k,'_',' & '), pairTypeCount(k)), pairKeys, 'UniformOutput', false), ', '));
reportLines{end+1} = sprintf('Condition x BlockNumber confound detected: %d', confoundFlag);
reportLines{end+1} = '';
reportLines{end+1} = 'Model: Linear mixed-effects, Condition (fixed, 3-level, ref=E) + (1|Subject), REML, Satterthwaite DF.';
reportLines{end+1} = 'Metrics: LogAbsolutePower and RelativePower, co-primary, parallel, never pooled for FDR.';
reportLines{end+1} = 'Pairwise convention (all 3 levels): N-E, R-E, R-N (X05/X11-X15), Holm-corrected, gated on FDR.';
reportLines{end+1} = '';
reportLines{end+1} = '--- GLOBAL LEVEL ---';
reportLines{end+1} = sprintf('Omnibus tests: %d / %d expected. Significant after tiered FDR: %d.', ...
    height(GlobalOmnibusT), nBands*numel(METRICS), sum(GlobalOmnibusT.SignificantFDR));
for r = 1:height(GlobalOmnibusT)
    reportLines{end+1} = sprintf('  [%s] %-25s Tier=%-11s F(%g,%.1f)=%.3f p=%.4f pFDR=%.4f family_n=%d %s', ...
        GlobalOmnibusT.Metric{r}, GlobalOmnibusT.Band{r}, GlobalOmnibusT.Tier{r}, GlobalOmnibusT.DF1(r), GlobalOmnibusT.DF2(r), ...
        GlobalOmnibusT.FStat(r), GlobalOmnibusT.pValue(r), GlobalOmnibusT.pFDR(r), GlobalOmnibusT.FamilySize(r), ...
        local_yn(GlobalOmnibusT.SignificantFDR(r))); %#ok<AGROW>
end
reportLines{end+1} = '';
if isempty(GlobalPairwiseT)
    reportLines{end+1} = '  Pairwise: (none -- no Band/Metric survived tiered FDR)';
else
    reportLines{end+1} = '  Pairwise contrasts:';
    for r = 1:height(GlobalPairwiseT)
        reportLines{end+1} = sprintf('    [%s] %-25s %-8s est=%.4f (95%% CI %.4f to %.4f) t=%.3f df=%.1f p=%.4f pHolm=%.4f', ...
            GlobalPairwiseT.Metric{r}, GlobalPairwiseT.Band{r}, GlobalPairwiseT.Contrast{r}, GlobalPairwiseT.Estimate(r), ...
            GlobalPairwiseT.CI_Low(r), GlobalPairwiseT.CI_High(r), GlobalPairwiseT.tStat(r), GlobalPairwiseT.DF(r), ...
            GlobalPairwiseT.pValue(r), GlobalPairwiseT.pHolm(r)); %#ok<AGROW>
    end
end
reportLines{end+1} = '';
reportLines{end+1} = '--- ROI LEVEL ---';
reportLines{end+1} = sprintf('Omnibus tests: %d / %d expected. Significant after Band x Metric FDR (across %d ROIs): %d.', ...
    height(RoiOmnibusT), nROI*nBands*numel(METRICS), nROI, sum(RoiOmnibusT.SignificantFDR));
sigRowsROI = find(RoiOmnibusT.SignificantFDR);
if isempty(sigRowsROI)
    reportLines{end+1} = '  (no ROI x Band x Metric survived FDR)';
else
    for ri = sigRowsROI'
        reportLines{end+1} = sprintf('  [%s] %-10s %-25s F(%g,%.1f)=%.3f p=%.4f pFDR=%.4f', ...
            RoiOmnibusT.Metric{ri}, RoiOmnibusT.ROI{ri}, RoiOmnibusT.Band{ri}, RoiOmnibusT.DF1(ri), RoiOmnibusT.DF2(ri), ...
            RoiOmnibusT.FStat(ri), RoiOmnibusT.pValue(ri), RoiOmnibusT.pFDR(ri)); %#ok<AGROW>
    end
end
reportLines{end+1} = sprintf('  Pairwise contrast rows: %d', height(RoiPairwiseT));
reportLines{end+1} = '';
reportLines{end+1} = '--- ELECTRODE LEVEL ---';
reportLines{end+1} = sprintf('Omnibus tests: %d OK, %d failed, out of %d expected. Significant after Band x Metric FDR (across %d electrodes): %d.', ...
    nOKElec, nFailedElec, nElec*nBands*numel(METRICS), nElec, sum(ElectrodeOmnibusT.SignificantFDR));
sigRowsElec = find(ElectrodeOmnibusT.SignificantFDR);
if isempty(sigRowsElec)
    reportLines{end+1} = '  (no Electrode x Band x Metric survived FDR)';
else
    for ri = sigRowsElec'
        reportLines{end+1} = sprintf('  [%s] %-6s %-25s F(%g,%.1f)=%.3f p=%.4f pFDR=%.4f', ...
            ElectrodeOmnibusT.Metric{ri}, ElectrodeOmnibusT.Electrode{ri}, ElectrodeOmnibusT.Band{ri}, ...
            ElectrodeOmnibusT.DF1(ri), ElectrodeOmnibusT.DF2(ri), ElectrodeOmnibusT.FStat(ri), ElectrodeOmnibusT.pValue(ri), ElectrodeOmnibusT.pFDR(ri)); %#ok<AGROW>
    end
end
reportLines{end+1} = sprintf('  Pairwise contrast rows: %d', height(ElectrodePairwiseT));
reportLines{end+1} = '';
reportLines{end+1} = '--- JUDGMENT CALLS (disclosed) ---';
for i = 1:numel(PARAMS16.JudgmentCalls)
    reportLines{end+1} = ['  - ' PARAMS16.JudgmentCalls{i}]; %#ok<AGROW>
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
if ~isempty(convergenceWarningsGlobal)
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Global-level model convergence warnings ---';
    for i = 1:numel(convergenceWarningsGlobal)
        reportLines{end+1} = ['  - ' convergenceWarningsGlobal{i}]; %#ok<AGROW>
    end
end
reportLines{end+1} = '';
reportLines{end+1} = 'This script does not compare FullReconstructed to the existing Full/Nback/SART';
reportLines{end+1} = 'analyses -- that is a separate, not-yet-written step.';
reportLines{end+1} = '';
reportLines{end+1} = 'Output files written to:';
reportLines{end+1} = ['  ' STATS_DIR];
reportLines{end+1} = repmat('=', 1, 60);

reportText = strjoin(reportLines, '\n');
fprintf('%s\n', reportText);
fid = fopen(fullfile(REPORT_DIR, sprintf('X16_full_reconstructed_report_%s.txt', RUN_MODE_LABEL)), 'w');
fprintf(fid, '%s\n', reportText);
fclose(fid);

fprintf('\n=== X16 complete. Purely additive -- X03/X04/X05/X06/X07/X08/X09/X10/X11-X15 outputs are all untouched. ===\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%  (must follow ALL top-level executable statements per MATLAB's rule for
%  local functions in a script file. Self-contained: nothing here is
%  imported/shared from another file, exactly as every other script in
%  this codebase. local_subject_numeric, local_fdr_bh,
%  local_holm_bonferroni, local_condition_contrast, local_build_table,
%  local_yn are copied VERBATIM from X05_primary_cognitive_statistics.m.
%  local_fit_condition_omnibus, local_bh_fdr_or_raw, and local_run_pairwise
%  are adapted directly from X11_median_robust_family_cognitive.m's own
%  equivalents (see header item 4). local_make_key is a small new helper
%  needed only for Section 2's Nback/SART key matching.)
% ========================================================================

function k = local_make_key(subjStr, blockNum, condStr)
    k = sprintf('%s|%g|%s', subjStr, blockNum, condStr);
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
            error('Phase16:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name1);
        end
        cRow(idx1) = cRow(idx1) - 1;
    end
    if ~strcmp(level2, refLevel)
        name2 = sprintf('%s_%s', varName, level2);
        idx2 = find(strcmp(coeffNames, name2));
        if isempty(idx2)
            error('Phase16:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name2);
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
        error('Phase16:BuildTableColumnMismatch', ...
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
            error('Phase16:BuildTableColumnMismatch', 'Column %d not classified as text or numeric.', j);
        end
    end
    T = table(cols{:}, 'VariableNames', varNames);
end

function fr = local_fit_condition_omnibus(subT, metricName, conditionOrder)
% Fits Metric ~ Condition + (1|Subject), REML, and extracts the omnibus
% Satterthwaite F-test for Condition. Returns a struct rather than
% erroring out, so the many calling loops in this script can log a
% failure (via fr.errMsg) and continue -- adapted directly from
% X11_median_robust_family_cognitive.m's own equivalent.
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
% Family size 1 is reported as-is, never treated as pre-corrected -- same
% convention X05 already uses for its Primary tier (family size 1).
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
% identical arithmetic to X05/X08/X09's own inline pairwise-contrast
% block, adapted from X11's own factored-out version since this script
% repeats it across 3 separate analysis levels (Global/ROI/Electrode).
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
