%% ========================================================================
% X08_roi_level_statistics.m
%
% Regional band-power statistics for the five regions defined in X07 (same model, FDR
% across regions per band and measure, Holm across the three contrasts).
% Inputs : 04_subject_level, 07_roi_electrode_definitions.
% Outputs: 08_roi_stats/<RUN_MODE>/, 09_reports/.
% Run twice, with RUN_MODE = 'full_sample' (N = 30) and RUN_MODE = 'sensitivity_excl_noisy'
% (cleaned sample, N = 28; excluded subjects listed in EXCLUSION_CONFIG_intervention.m).
% Outputs are written to <OUT_DIR>/<folder>/<RUN_MODE>/ and reports to <OUT_DIR>/09_reports/.
% ========================================================================

clear; clc;

%% ------------------------------------------------------------------
%  0. SETTINGS, PATHS, EXISTENCE CHECKS
% ------------------------------------------------------------------

% >>> Must match the ANALYSIS_ROOT used in Phases 1-5D <<<
config;
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

SUBJ_DIR   = fullfile(ANALYSIS_ROOT, '04_subject_level');
DEF_DIR    = fullfile(ANALYSIS_ROOT, '07_roi_electrode_definitions');
STATS_DIR  = fullfile(ANALYSIS_ROOT, '08_roi_stats', RUN_MODE);
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');

bandPowerFile = fullfile(SUBJ_DIR, 'phase3_subject_level_bandpower.mat');
roiDefFile    = fullfile(DEF_DIR, 'phase5_channel_roi_definitions.mat');

if ~exist(STATS_DIR, 'dir');  mkdir(STATS_DIR);  end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

fprintf('========================================================\n');
fprintf('PHASE 5A -- PREDEFINED 5-ROI STATISTICAL ANALYSIS\n');
fprintf('========================================================\n');
fprintf('ANALYSIS_ROOT   : %s\n', ANALYSIS_ROOT);
fprintf('bandPowerFile    : %s\n', bandPowerFile);
fprintf('roiDefFile       : %s\n', roiDefFile);
fprintf('OUTPUT           : %s\n', STATS_DIR);
fprintf('--------------------------------------------------------\n\n');

if ~exist(bandPowerFile, 'file')
    error('Phase5A:MissingInput', ['Phase 3 band-power file not found: %s\n' ...
        'Run 04_aggregate_subject_level.m first.'], bandPowerFile);
end
if ~exist(roiDefFile, 'file')
    error('Phase5A:MissingInput', ['Phase 5D ROI definition file not found: %s\n' ...
        'Run 07_roi_electrode_definitions.m first. This script will NOT ' ...
        'reconstruct or guess ROI membership on its own.'], roiDefFile);
end
if isempty(ver('stats'))
    error('Phase5A:MissingToolbox', ...
        'Statistics and Machine Learning Toolbox not detected. This script requires fitlme, anova(LinearMixedModel), coefTest, and qqplot.');
end

%% ------------------------------------------------------------------
%  1. LOAD PHASE 3 (Structure B) AND PHASE 5D (locked ROI definitions)
% ------------------------------------------------------------------

S3 = load(bandPowerFile, 'BandPowerT', 'PARAMS3');
BandPowerT = S3.BandPowerT;
BandPowerT_beforeRunMode = BandPowerT;
[BandPowerT, RUN_MODE_N_EXCL_ROWS, RUN_MODE_EXCL_SUBJ] = apply_run_mode_filter(BandPowerT, RUN_MODE, EXCLUDE_NOISY_SUBJECTS);
EXPC = expected_condition_counts(BandPowerT_beforeRunMode, RUN_MODE_EXCL_SUBJ, EXPECTED_N_SUBJECTS_TOTAL);
clear BandPowerT_beforeRunMode;
fprintf('Expected recordings per condition (RUN_MODE=%s): E=%d, N=%d, R=%d\n', RUN_MODE_LABEL, EXPC.E, EXPC.N, EXPC.R);
if RUN_MODE_N_EXCL_ROWS > 0; fprintf('RUN_MODE=%s: removed %d row(s) of subject(s) {%s} from BandPowerT.\n', RUN_MODE_LABEL, RUN_MODE_N_EXCL_ROWS, num2str(RUN_MODE_EXCL_SUBJ)); end

PARAMS3 = S3.PARAMS3;
clear S3;

requiredCols = {'Subject','InterventionNumber','Condition','Channel','Band', ...
                'AbsolutePower','LogAbsolutePower','RelativePower'};
missingCols = requiredCols(~ismember(requiredCols, BandPowerT.Properties.VariableNames));
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
CONDITION_ORDER = {'E', 'N', 'R'};
REF_LEVEL = CONDITION_ORDER{1};   % 'E' is the reference level for dummy coding -- identical to Phase 4

PRIMARY_TIER_BANDS     = {'Alpha'};
SECONDARY_TIER_BANDS   = {'Delta', 'Theta', 'Beta'};
EXPLORATORY_TIER_BANDS = {'Gamma_low_exploratory', 'Gamma_high1_exploratory', 'Gamma_high2_exploratory'};
EXPLORATORY_BETA_TIER_BANDS = {'Beta_low_exploratory', 'Beta_high_exploratory'};
expectedTierBands = [PRIMARY_TIER_BANDS, SECONDARY_TIER_BANDS, EXPLORATORY_TIER_BANDS, EXPLORATORY_BETA_TIER_BANDS];
if ~isempty(setxor(sort(expectedTierBands), sort(ALL_BANDS')))
    error('Phase5A:BandMismatch', ...
        ['The hard-coded band tiers do not exactly match PARAMS3.BandDefinitionsHz. ' ...
         'Tiers: %s | Phase3 bands: %s'], strjoin(expectedTierBands, ','), strjoin(ALL_BANDS', ','));
end

nBands = numel(ALL_BANDS);
nMetrics = numel(METRICS);
nROI = numel(ROI_NAMES);

integrityIssues = {};

%% ========================================================================
%  2. PRE-AGGREGATION VALIDATION (all performed before any model is fit)
% ========================================================================
fprintf('--------------------------------------------------------\n');
fprintf('SECTION 2: PRE-AGGREGATION VALIDATION\n');
fprintf('--------------------------------------------------------\n');
preValFails = 0;

% --- 2a. All five ROI definitions loaded from Phase 5D ---
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

% --- 2b. All expected channels present in BandPowerT ---
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
    fprintf('  *** FAIL: %d channel(s) in BandPowerT are not in the 62-channel reference: %s\n', ...
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
    fprintf('  PASS: P9 and P10 are absent from all five ROI_DEFINITION lists (confirmed excluded from ROI averaging).\n');
end
if ~all(ismember(upper({'P9','P10'}), upper(UNASSIGNED_FOR_ROI)))
    fprintf('  *** FAIL: P9/P10 not both present in UNASSIGNED_FOR_ROI as expected.\n');
    preValFails = preValFails + 1;
else
    fprintf('  PASS: P9 and P10 are both explicitly recorded as unassigned-for-ROI (per Phase 5D).\n');
end

comboT = unique(BandPowerT(:, {'Subject','InterventionNumber','Condition'}), 'rows');
nRecordings = height(comboT);
fprintf('\n  %d unique Subject x InterventionNumber x Condition recordings found.\n', nRecordings);

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
        integrityIssues{end+1} = sprintf('Subject %g has %d recordings (expected 2).', uniqSubj(i), numel(condsI)); %#ok<AGROW>
        preValFails = preValFails + 1;
    elseif strcmp(condsI{1}, condsI{2})
        integrityIssues{end+1} = sprintf('Subject %g has the SAME condition (%s) in both recordings -- duplicate.', uniqSubj(i), condsI{1}); %#ok<AGROW>
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
fprintf('  Condition-pair coverage across %d subjects:\n', nSubjects);
pairKeys = keys(pairTypeCount);
for k = 1:numel(pairKeys)
    fprintf('    %s : %d subjects\n', strrep(pairKeys{k}, '_', ' & '), pairTypeCount(pairKeys{k}));
end
if duplicateSubjCondFlag
    fprintf('  *** FAIL: at least one subject has a duplicated condition. See integrity issues.\n');
else
    fprintf('  PASS: no subject has a duplicated condition (each subject contributes exactly 2 DIFFERENT conditions).\n');
end

condCounts = struct();
for c = 1:numel(CONDITION_ORDER)
    condCounts.(CONDITION_ORDER{c}) = sum(strcmp(cellstr(comboT.Condition), CONDITION_ORDER{c}));
end
fprintf('  Recording-level condition counts: E=%d, N=%d, R=%d\n', condCounts.E, condCounts.N, condCounts.R);
if condCounts.E ~= EXPC.E || condCounts.N ~= EXPC.N || condCounts.R ~= EXPC.R
    fprintf('  *** FAIL: expected E=%d, N=%d, R=%d. Got E=%d, N=%d, R=%d.\n', EXPC.E, EXPC.N, EXPC.R, condCounts.E, condCounts.N, condCounts.R);
    integrityIssues{end+1} = sprintf('Recording-level condition counts do not match expected E=%d/N=%d/R=%d (got E=%d/N=%d/R=%d).', ...
        EXPC.E, EXPC.N, EXPC.R, condCounts.E, condCounts.N, condCounts.R); %#ok<AGROW>
    preValFails = preValFails + 1;
else
    fprintf('  PASS: recording-level condition counts match the expected E=%d, N=%d, R=%d.\n', EXPC.E, EXPC.N, EXPC.R);
end

fprintf('\nPre-aggregation validation failures: %d\n', preValFails);
if preValFails > 0
    error('Phase5A:PreValidationFailed', ...
        'Pre-aggregation validation failed (%d issue(s) above). Fix before aggregating/fitting.', preValFails);
end

%% ========================================================================
%  3. ROI AGGREGATION (Structure B -> ROI-level long table)
%     Arithmetic-mean-of-channels rule; LogAbsolutePower = log(mean
%     AbsolutePower), matching Phase 3's own Structure C construction
%     exactly. P9/P10 are never included (their rows are
%     simply never referenced by any ROI's channel list).
% ========================================================================
fprintf('\n--------------------------------------------------------\n');
fprintf('SECTION 3: ROI AGGREGATION\n');
fprintf('--------------------------------------------------------\n');

BandPowerT.ChannelUC = upper(BandPowerT.Channel);

roiRows = cell(0, 7);   % Subject, InterventionNumber, Condition, ROI, Band, AbsolutePower, RelativePower
roiCols = {'Subject','InterventionNumber','Condition','ROI','Band','AbsolutePower','RelativePower'};

recCombos = unique(BandPowerT(:, {'Subject','InterventionNumber','Condition'}), 'rows');
nonPositiveROIAbs = {};

for i = 1:height(recCombos)
    subjStr = recCombos.Subject{i};
    intNum  = recCombos.InterventionNumber(i);
    condStr = recCombos.Condition{i};
    recMask = strcmp(BandPowerT.Subject, subjStr) & BandPowerT.InterventionNumber == intNum & strcmp(BandPowerT.Condition, condStr);

    for r = 1:nROI
        roiName = ROI_NAMES{r};
        roiChansUC = upper(ROI_DEFINITION.(roiName));
        chanMask = ismember(BandPowerT.ChannelUC, roiChansUC);

        for b = 1:nBands
            bandName = ALL_BANDS{b};
            rowMask = recMask & chanMask & strcmp(BandPowerT.Band, bandName);
            nChanFound = sum(rowMask);
            if nChanFound ~= numel(roiChansUC)
                integrityIssues{end+1} = sprintf('Sub%s Int%g Cond%s ROI=%s Band=%s: found %d/%d expected ROI channels.', ...
                    subjStr, intNum, condStr, roiName, bandName, nChanFound, numel(roiChansUC)); %#ok<AGROW>
            end
            roiAbs = mean(BandPowerT.AbsolutePower(rowMask));
            roiRel = mean(BandPowerT.RelativePower(rowMask));
            if roiAbs <= 0
                nonPositiveROIAbs{end+1} = sprintf('Sub%s Int%g Cond%s ROI=%s Band=%s', subjStr, intNum, condStr, roiName, bandName); %#ok<AGROW>
            end
            roiRows(end+1, :) = {subjStr, intNum, condStr, roiName, bandName, roiAbs, roiRel}; %#ok<AGROW>
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

fprintf('ROI aggregation complete: %d rows (%d recordings x %d ROIs x %d bands).\n', ...
    height(RoiBandPowerT), height(recCombos), nROI, nBands);

%% ========================================================================
%  4. POST-AGGREGATION VALIDATION (still before any model is fit)
% ========================================================================
fprintf('\n--------------------------------------------------------\n');
fprintf('SECTION 4: POST-AGGREGATION VALIDATION\n');
fprintf('--------------------------------------------------------\n');
postValFails = 0;

% --- 4a. Every Subject x Intervention x ROI x Band observation unique
%          (Metric is a column, not a row-defining key, so uniqueness of
%          this key set is exactly "every ... x Metric observation is
%          uniquely defined": both metric columns exist once per key). ---
keyT = RoiBandPowerT(:, {'Subject','InterventionNumber','ROI','Band'});
[uKeyT, ~, ic] = unique(keyT, 'rows');
dupCount = height(keyT) - height(uKeyT);
if dupCount > 0
    fprintf('  *** FAIL: %d duplicate Subject x InterventionNumber x ROI x Band row(s) found.\n', dupCount);
    postValFails = postValFails + 1;
else
    fprintf('  PASS: every Subject x InterventionNumber x ROI x Band x Metric observation is uniquely defined (%d rows, %d unique keys).\n', ...
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

% --- 4c. Per-ROI, per-band recording-level condition counts still E=20/N=20/R=20 ---
condCountMismatch = false;
for r = 1:nROI
    for b = 1:nBands
        subT = RoiBandPowerT(strcmp(RoiBandPowerT.ROI, ROI_NAMES{r}) & strcmp(RoiBandPowerT.Band, ALL_BANDS{b}), :);
        cE = sum(strcmp(subT.Condition, 'E')); cN = sum(strcmp(subT.Condition, 'N')); cR = sum(strcmp(subT.Condition, 'R'));
        if cE ~= EXPC.E || cN ~= EXPC.N || cR ~= EXPC.R
            condCountMismatch = true;
            integrityIssues{end+1} = sprintf('ROI=%s Band=%s: condition counts E=%d/N=%d/R=%d (expected %d/%d/%d).', ...
                ROI_NAMES{r}, ALL_BANDS{b}, cE, cN, cR, EXPC.E, EXPC.N, EXPC.R); %#ok<AGROW>
        end
    end
end
if condCountMismatch
    fprintf('  *** FAIL: at least one ROI x Band combination does not preserve the expected E/N/R counts. See integrity issues.\n');
    postValFails = postValFails + 1;
else
    fprintf('  PASS: every ROI x Band combination (%d x %d = %d checks) preserves E=%d, N=%d, R=%d recordings.\n', nROI, nBands, nROI*nBands, EXPC.E, EXPC.N, EXPC.R);
end

fprintf('\nPost-aggregation validation failures: %d\n', postValFails);
if postValFails > 0
    error('Phase5A:PostValidationFailed', ...
        'Post-aggregation validation failed (%d issue(s) above). Fix before fitting models.', postValFails);
end

%% ========================================================================
%  5. MODEL FITTING -- 70 independent LMMs (5 ROI x 7 bands x 2 metrics)
% ========================================================================
fprintf('\n--------------------------------------------------------\n');
fprintf('SECTION 5: MODEL FITTING (ROI x Band x Metric, %d total fits)\n', nROI*nBands*nMetrics);
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
                    convergenceWarnings{end+1} = sprintf('%s / %s / %s: %s', roiName, metricName, bandName, warnMsg); %#ok<AGROW>
                end
            catch ME
                integrityIssues{end+1} = sprintf('fitlme FAILED for ROI=%s %s / %s: %s', roiName, metricName, bandName, ME.message); %#ok<AGROW>
                lmeStore{r, b, m} = [];
                continue;
            end
            lmeStore{r, b, m} = lme;
            aov = anova(lme, 'DFMethod', 'satterthwaite');
            rowIdx = find(strcmp(aov.Term, 'Condition'), 1);
            if isempty(rowIdx)
                integrityIssues{end+1} = sprintf('anova(lme) has no Condition row for ROI=%s %s / %s.', roiName, metricName, bandName); %#ok<AGROW>
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
for i = 1:numel(EXPLORATORY_BETA_TIER_BANDS); tierOf(EXPLORATORY_BETA_TIER_BANDS{i}) = 'ExploratoryBeta'; end
OmnibusT.Tier = cellfun(@(b) tierOf(b), OmnibusT.Band, 'UniformOutput', false);

fprintf('Model fitting complete: %d / %d omnibus tests fit successfully.\n', height(OmnibusT), nROI*nBands*nMetrics);

%% ========================================================================
%  6. FDR -- separately within each Band x Metric family, across the 5 ROIs
% ========================================================================
fprintf('\n--------------------------------------------------------\n');
fprintf('SECTION 6: BH-FDR (Band x Metric family, across 5 ROIs)\n');
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
            integrityIssues{end+1} = sprintf('FDR family %s / %s has %d members, expected %d (one per ROI) -- some ROI fit may have failed.', ...
                metricName, bandName, nInFamily, nROI); %#ok<AGROW>
        end
        pvals = OmnibusT.pValue(idx);
        [adjP, reject] = local_fdr_bh(pvals, FDR_ALPHA);
        OmnibusT.pFDR(idx) = adjP;
        OmnibusT.SignificantFDR(idx) = reject;
    end
end
fprintf('BH-FDR applied within %d Band x Metric families (expected %d = %d bands x %d metrics), each spanning up to %d ROIs.\n', ...
    nFamilies, nBands*nMetrics, nBands, nMetrics, nROI);

fprintf('\n--- ROI-level omnibus results ---\n');
disp(OmnibusT);

%% ========================================================================
%  7. PAIRWISE CONTRASTS (E-N, E-R, N-R), gated on the Band x Metric FDR
%     result, Holm-corrected within each surviving ROI x Band x Metric
% ========================================================================
fprintf('\n--------------------------------------------------------\n');
fprintf('SECTION 7: PAIRWISE CONTRASTS (gated on FDR, Holm-corrected)\n');
fprintf('--------------------------------------------------------\n');

pairwiseRows = cell(0, 12);   % ROI,Metric,Band,Contrast,Estimate,SE,DF,tStat,pValue,pHolm,CI_Low,CI_High
CONTRASTS = {'E','N','E - N'; 'E','R','E - R'; 'N','R','N - R'};   % {lvl2, lvl1, label}; label = lvl2 - lvl1

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

fprintf('\n--- Pairwise contrasts (ROI x Band x Metric significant under Band x Metric FDR) ---\n');
disp(PairwiseT);

allRows = cell(0, 13);
for r = 1:nROI
    for b = 1:nBands
        for m = 1:nMetrics
            lme = lmeStore{r, b, m};
            if isempty(lme); continue; end
            cp = nan(3, 1); cdat = cell(3, 8);
            for c = 1:3
                C = local_condition_contrast(lme.CoefficientNames, 'Condition', REF_LEVEL, CONTRASTS{c, 2}, CONTRASTS{c, 1});
                est = C * lme.Coefficients.Estimate; se = sqrt(C * lme.CoefficientCovariance * C');
                [p, ~, ~, df2c] = coefTest(lme, C, 0, 'DFMethod', 'satterthwaite');
                tc = tinv(0.975, df2c);
                cp(c) = p; cdat(c, :) = {CONTRASTS{c, 3}, est, se, df2c, est / se, p, est - tc*se, est + tc*se};
            end
            adj = local_holm_bonferroni(cp, FDR_ALPHA);
            omniRow = strcmp(OmnibusT.ROI, ROI_NAMES{r}) & strcmp(OmnibusT.Band, ALL_BANDS{b}) & strcmp(OmnibusT.Metric, METRICS{m});
            gatePassed = any(OmnibusT.SignificantFDR(omniRow));
            for c = 1:3
                allRows(end+1, :) = [{ROI_NAMES{r}, METRICS{m}, ALL_BANDS{b}}, cdat(c, 1:6), {adj(c)}, cdat(c, 7:8), {double(gatePassed)}]; %#ok<AGROW>
            end
        end
    end
end
AllContrastsT = local_build_table(allRows, [1 2 3 4], [5 6 7 8 9 10 11 12 13], ...
    {'ROI','Metric','Band','Contrast','Estimate','SE','DF','tStat','pValue','pHolm_withinModel','CI_Low','CI_High','OmnibusPassedFDR'});

%% ------------------------------------------------------------------
%  8. MODEL DIAGNOSTICS (Primary + Secondary tier bands only, all 5 ROIs,
%     both metrics -- Alpha/Delta/Theta/Beta x 5 ROI x 2 metrics = 40 plots.
%     Exploratory/gamma bands are skipped here to keep diagnostic output
%     scoped to the bands carrying interpretive weight, matching Phase 4's
%     own diagBands convention.)
% ------------------------------------------------------------------
fprintf('\n--------------------------------------------------------\n');
fprintf('SECTION 8: MODEL DIAGNOSTICS\n');
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
                title(sprintf('Residuals: %s / %s (%s)', roiName, bandName, metricName));
                xlabel('Pearson residual'); ylabel('Count');
                subplot(1,2,2);
                qqplot(res);
                title('Normal Q-Q');
                outPng = fullfile(STATS_DIR, sprintf('diagnostic_residuals_%s_%s_%s.png', roiName, bandName, metricName));
                saveas(fig, outPng);
                close(fig);
                nDiagSaved = nDiagSaved + 1;
            catch ME
                integrityIssues{end+1} = sprintf('Residual diagnostic failed for ROI=%s %s / %s: %s', roiName, metricName, bandName, ME.message); %#ok<AGROW>
            end
        end
    end
end
fprintf('Saved %d residual diagnostic plots (Primary+Secondary bands x 5 ROIs x 2 metrics; expected up to %d).\n', ...
    nDiagSaved, nROI*numel(diagBands)*nMetrics);

%% ------------------------------------------------------------------
%  9. PARAMETERS AND OUTPUT
% ------------------------------------------------------------------

PARAMS5A = struct();
PARAMS5A.ScriptName = '08_roi_level_statistics.m';
PARAMS5A.GeneratedOn = datestr(now);
PARAMS5A.AnalysisLevel = 'ROI-averaged band power (Phase 5D locked 5-ROI definition, applied to Phase 3 Structure B / BandPowerT)';
PARAMS5A.Model = 'LinearMixedEffects_RandomIntercept';
PARAMS5A.FixedEffectsPrimary = 'Condition (3-level: E, N, R; reference = E)';
PARAMS5A.RandomEffects = '(1|Subject)';
PARAMS5A.FitMethod = 'REML';
PARAMS5A.DFMethod = 'Satterthwaite';
PARAMS5A.Metrics = METRICS;
PARAMS5A.ROINames = ROI_NAMES;
PARAMS5A.UnassignedForROI = {UNASSIGNED_FOR_ROI};
PARAMS5A.ROIAggregationRule = 'Arithmetic mean of Structure B per-channel AbsolutePower/RelativePower across ROI member channels; LogAbsolutePower = log(mean AbsolutePower), matching Phase 3 Structure C construction.';
PARAMS5A.BandTiers = struct('Primary', {PRIMARY_TIER_BANDS}, 'Secondary', {SECONDARY_TIER_BANDS}, 'Exploratory', {EXPLORATORY_TIER_BANDS}, 'ExploratoryBeta', {EXPLORATORY_BETA_TIER_BANDS});
PARAMS5A.FDRAlpha = FDR_ALPHA;
PARAMS5A.FDRFamilyScheme = 'Separate BH-FDR family per (Band x Metric), across the 5 ROIs (family size = 5). NOT tier-based. Bands/metrics/ROIs never pooled.';
PARAMS5A.PairwiseContrastMethod = 'Holm-Bonferroni across the 3 pairwise contrasts (E-N, E-R, N-R), computed only for ROI x Band x Metric combinations whose Band x Metric FDR family test survives BH-FDR.';
PARAMS5A.EffectSizePolicy = 'Model-based estimated difference + 95% CI only. No standardized effect-size denominator computed in this script (decision deferred).';
PARAMS5A.SourcePhase3File = bandPowerFile;
PARAMS5A.SourcePhase5DFile = roiDefFile;
PARAMS5A.SourcePhase5DDecision = PARAMS5D.ROIDecision;
PARAMS5A.DesignVerification = struct('nRecordings', nRecordings, 'nSubjects', nSubjects, ...
    'PairTypeCounts', {pairTypeCount.keys}, 'ConditionCounts', condCounts, 'IntegrityIssues', {integrityIssues});
PARAMS5A.RunDateTime = datestr(now);

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
save(fullfile(STATS_DIR, 'phase5a_roi_all_contrasts_ungated.mat'), 'AllContrastsT');
writetable(AllContrastsT, fullfile(STATS_DIR, 'phase5a_roi_all_contrasts_ungated.csv'));

fprintf('\nSaved output to: %s\n', STATS_DIR);

%% ------------------------------------------------------------------
%  10. REPORT
% ------------------------------------------------------------------

reportLines = {};
reportLines{end+1} = sprintf('RUN_MODE = %s | subjects excluded in this run: %s', RUN_MODE_LABEL, local_v3_list(RUN_MODE_EXCL_SUBJ));
reportLines{end+1} = repmat('=', 1, 60);
reportLines{end+1} = 'PHASE 5A -- ROI-LEVEL E vs N vs R STATISTICS REPORT';
reportLines{end+1} = repmat('=', 1, 60);
reportLines{end+1} = '';
reportLines{end+1} = sprintf('Recordings verified: %d | Subjects: %d', nRecordings, nSubjects);
reportLines{end+1} = sprintf('Condition-pair coverage: %s', strjoin(cellfun(@(k) sprintf('%s=%d', strrep(k,'_',' & '), pairTypeCount(k)), pairKeys, 'UniformOutput', false), ', '));
reportLines{end+1} = sprintf('Recording-level condition counts: E=%d, N=%d, R=%d', condCounts.E, condCounts.N, condCounts.R);
reportLines{end+1} = '';
reportLines{end+1} = 'Model: Linear mixed-effects, Condition (fixed, 3-level) + (1|Subject), REML, Satterthwaite DF.';
reportLines{end+1} = 'Analysis level: ROI-averaged band power (Phase 5D locked 5-ROI definition; Phase 3 Structure B).';
reportLines{end+1} = sprintf('ROIs (%d): %s', nROI, strjoin(ROI_NAMES, ', '));
reportLines{end+1} = sprintf('Unassigned for ROI averaging: %s (still analyzed individually at the electrode level, Phase 5B, not in this script)', strjoin(UNASSIGNED_FOR_ROI, ', '));
reportLines{end+1} = 'Tiers (labeling only, not FDR family boundaries): Primary=Alpha; Secondary=Delta/Theta/Beta; Exploratory=3 gamma bands; ExploratoryBeta=Beta_low 13-20 + Beta_high 20-30 (own family).';
reportLines{end+1} = 'Metrics: LogAbsolutePower and RelativePower, co-primary, parallel, never pooled for FDR.';
reportLines{end+1} = 'FDR: BH-FDR within each Band x Metric family, across the 5 ROIs (14 families x 5 ROIs = 70 tests).';
reportLines{end+1} = 'Pairwise: E-N, E-R, N-R; computed only for FDR-surviving ROI x Band x Metric; Holm-corrected within each (family size 3).';
reportLines{end+1} = 'Effect sizes: model-based estimated difference + 95% CI only. No standardized effect size in this script.';
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
reportLines{end+1} = '--- [v2] ALL contrasts for models whose omnibus did NOT survive FDR (UNGATED, descriptive only) ---';
reportLines{end+1} = '    (models that passed FDR are already listed above; these rows carry no inferential claim)';
for row = find(AllContrastsT.OmnibusPassedFDR == 0)'
    reportLines{end+1} = sprintf('  [%s] %-10s %-25s %-6s est=%.4f (95%% CI %.4f to %.4f) p=%.4f pHolm(within model)=%.4f', ...
        AllContrastsT.Metric{row}, AllContrastsT.ROI{row}, AllContrastsT.Band{row}, AllContrastsT.Contrast{row}, ...
        AllContrastsT.Estimate(row), AllContrastsT.CI_Low(row), AllContrastsT.CI_High(row), ...
        AllContrastsT.pValue(row), AllContrastsT.pHolm_withinModel(row)); %#ok<AGROW>
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
reportLines{end+1} = 'specified above. Electrode-level analysis is 09_electrode_level_statistics.m;';
reportLines{end+1} = 'topography is 10_descriptive_topography.m; literature-based ROIs are 11/12.';
reportLines{end+1} = '';
reportLines{end+1} = 'Output files written to:';
reportLines{end+1} = ['  ' STATS_DIR];
reportLines{end+1} = repmat('=', 1, 60);

reportText = strjoin(reportLines, '\n');
fprintf('%s\n', reportText);
fid = fopen(fullfile(REPORT_DIR, sprintf('phase5a_roi_level_statistics_report_%s.txt', RUN_MODE_LABEL)), 'w');
fprintf(fid, '%s\n', reportText);
fclose(fid);

fprintf('\n=== 08_roi_level_statistics.m complete. STOP for review before writing 09_electrode_level_statistics.m. ===\n');

%% ======================= LOCAL FUNCTIONS (all after top-level code) =======================

function n = local_subject_numeric(subjVal)
% Robustly extract a numeric subject id regardless of whether Subject is
% stored as numeric, categorical, char, or cellstr, and regardless of
% zero-padding (e.g. '7' vs '07'). Identical implementation to Phase 4.
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
% Phase 4 / Phase 4 gamma robustness check.
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
% implementation to Phase 4.
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
% column position/order). Identical implementation to Phase 4.
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
% Identical implementation to Phase 4 / Phase 4 gamma robustness check.
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

function s = local_v3_list(v)
    if isempty(v); s = 'none'; else; s = strjoin(cellstr(num2str(v(:))), ', '); end
end
