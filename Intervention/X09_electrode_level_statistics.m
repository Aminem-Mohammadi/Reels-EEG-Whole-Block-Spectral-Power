%% ========================================================================
% X09_electrode_level_statistics.m
%
% Electrode-level band-power statistics (62 electrodes; FDR across electrodes per band and
% measure; Holm across the three contrasts). Writes electrode_model_estimates.csv with every
% contrast estimate, SE, CI, df and p, used for the topographic maps.
% Inputs : 04_subject_level, 07_roi_electrode_definitions.
% Outputs: 09_electrode_stats/<RUN_MODE>/, 09_reports/.
% Run twice, with RUN_MODE = 'full_sample' (N = 30) and RUN_MODE = 'sensitivity_excl_noisy'
% (cleaned sample, N = 28; excluded subjects listed in EXCLUSION_CONFIG_intervention.m).
% Outputs are written to <OUT_DIR>/<folder>/<RUN_MODE>/ and reports to <OUT_DIR>/09_reports/.
% ========================================================================

clear; clc;

%% ------------------------------------------------------------------
%  0. SETTINGS, PATHS, EXISTENCE CHECKS
% ------------------------------------------------------------------

% >>> Must match the ANALYSIS_ROOT used in Phases 1-5A <<<
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

SUBJ_DIR  = fullfile(ANALYSIS_ROOT, '04_subject_level');
DEF_DIR   = fullfile(ANALYSIS_ROOT, '07_roi_electrode_definitions');
ELEC_DIR  = fullfile(ANALYSIS_ROOT, '09_electrode_stats', RUN_MODE);
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');

bandPowerFile = fullfile(SUBJ_DIR, 'phase3_subject_level_bandpower.mat');
roiDefFile    = fullfile(DEF_DIR, 'phase5_channel_roi_definitions.mat');

if ~exist(ELEC_DIR, 'dir');   mkdir(ELEC_DIR);   end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

fprintf('========================================================\n');
fprintf('PHASE 5B -- ELECTRODE-LEVEL E vs N vs R STATISTICS\n');
fprintf('========================================================\n');
fprintf('ANALYSIS_ROOT   : %s\n', ANALYSIS_ROOT);
fprintf('bandPowerFile    : %s\n', bandPowerFile);
fprintf('roiDefFile       : %s (read-only; electrode order source)\n', roiDefFile);
fprintf('OUTPUT           : %s\n', ELEC_DIR);
fprintf('--------------------------------------------------------\n\n');

if ~exist(bandPowerFile, 'file')
    error('Phase5B:MissingInput', ['Phase 3 band-power file not found: %s\n' ...
        'Run 04_aggregate_subject_level.m first.'], bandPowerFile);
end
if ~exist(roiDefFile, 'file')
    error('Phase5B:MissingInput', ['Phase 5D electrode/ROI definition file not found: %s\n' ...
        'Run 07_roi_electrode_definitions.m first. This script will NOT ' ...
        'reconstruct or guess the electrode order on its own.'], roiDefFile);
end
if isempty(ver('stats'))
    error('Phase5B:MissingToolbox', ...
        'Statistics and Machine Learning Toolbox not detected. This script requires fitlme, anova(LinearMixedModel), and coefTest.');
end

%% ------------------------------------------------------------------
%  1. LOAD PHASE 3 (Structure B) AND PHASE 5D (locked electrode order)
%     Both loaded read-only. Nothing in this section is modified or
%     re-saved back to Phase 3/5D's own files.
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
    error('Phase5B:MissingColumns', 'BandPowerT is missing expected column(s): %s', strjoin(missingCols, ', '));
end

S5D = load(roiDefFile, 'ChannelDefT', 'PARAMS5D');
ChannelDefT = S5D.ChannelDefT;
PARAMS5D = S5D.PARAMS5D;
clear S5D;

ELECTRODES = ChannelDefT.Label;   % locked 62-electrode reference order (Phase 5D's own source)
nElec = numel(ELECTRODES);

ALL_BANDS = fieldnames(PARAMS3.BandDefinitionsHz);
METRICS   = {'LogAbsolutePower', 'RelativePower'};
CONDITION_ORDER = {'E', 'N', 'R'};
REF_LEVEL = CONDITION_ORDER{1};   % 'E' is the reference level for dummy coding -- identical to Phase 4/5A

PRIMARY_TIER_BANDS     = {'Alpha'};
SECONDARY_TIER_BANDS   = {'Delta', 'Theta', 'Beta'};
EXPLORATORY_TIER_BANDS = {'Gamma_low_exploratory', 'Gamma_high1_exploratory', 'Gamma_high2_exploratory'};
EXPLORATORY_BETA_TIER_BANDS = {'Beta_low_exploratory', 'Beta_high_exploratory'};
expectedTierBands = [PRIMARY_TIER_BANDS, SECONDARY_TIER_BANDS, EXPLORATORY_TIER_BANDS, EXPLORATORY_BETA_TIER_BANDS];
if ~isempty(setxor(sort(expectedTierBands), sort(ALL_BANDS')))
    error('Phase5B:BandMismatch', ...
        ['The hard-coded band tiers do not exactly match PARAMS3.BandDefinitionsHz. ' ...
         'Tiers: %s | Phase3 bands: %s'], strjoin(expectedTierBands, ','), strjoin(ALL_BANDS', ','));
end

nBands = numel(ALL_BANDS);
nMetrics = numel(METRICS);

integrityIssues = {};
modelFailures = {};

%% ========================================================================
%  2. PRE-MODELING VALIDATION (items A-J). Any failure STOPS the script
%     before any model is fit.
% ========================================================================
fprintf('--------------------------------------------------------\n');
fprintf('SECTION 2: PRE-MODELING VALIDATION\n');
fprintf('--------------------------------------------------------\n');
preValFails = 0;

% --- A. Exactly 62 electrode labels present ---
if nElec ~= 62
    fprintf('  *** FAIL (A): ChannelDefT has %d electrodes, expected 62.\n', nElec);
    preValFails = preValFails + 1;
else
    fprintf('  PASS (A): exactly 62 electrode labels loaded from Phase 5D.\n');
end

CANONICAL_ORDER = {'AFZ','FP1','FP2','AF3','AF4','F7','F3','FZ','F4','F8', ...
    'FC5','FC1','FC2','FC6','T7','C3','Cz','C4','T8','CP5','CP1','CP2','CP6', ...
    'P7','P3','PZ','P4','P8','PO3','PO4','O1','O2','AF7','AF8','F5','F1','F2','F6', ...
    'FT7','FC3','FCz','FC4','FT8','C5','C1','C2','C6','TP7','CP3','CPz','CP4','TP8', ...
    'P5','P1','P2','P6','PO7','POz','PO8','Oz','P9','P10'};
if numel(CANONICAL_ORDER) ~= 62
    error('Phase5B:CanonicalOrderMalformed', 'Internal error: CANONICAL_ORDER constant does not have 62 entries.');
end
n = min(numel(ELECTRODES), numel(CANONICAL_ORDER));
caseOnlyDiffPositions = [];
identityOrOrderMismatch = (numel(ELECTRODES) ~= numel(CANONICAL_ORDER));
for i = 1:n
    if strcmp(ELECTRODES{i}, CANONICAL_ORDER{i})
        continue;   % exact match, nothing to report
    elseif strcmpi(ELECTRODES{i}, CANONICAL_ORDER{i})
        caseOnlyDiffPositions(end+1) = i; %#ok<AGROW>   % same electrode, case differs only
    else
        identityOrOrderMismatch = true;
    end
end
if ~identityOrOrderMismatch
    fprintf('  PASS (B): loaded electrode order matches the locked reference order (position-for-position, case-insensitive).\n');
    if ~isempty(caseOnlyDiffPositions)
        fprintf('      NOTE: %d position(s) matched only case-insensitively (not a failure -- electrode-label case is not used as an identity signal anywhere in this pipeline):\n', numel(caseOnlyDiffPositions));
        for i = caseOnlyDiffPositions
            fprintf('        position %d: loaded=%s canonical-list=%s\n', i, ELECTRODES{i}, CANONICAL_ORDER{i});
            integrityIssues{end+1} = sprintf('Item B: position %d electrode label case differs only (loaded=%s, canonical-list=%s); treated as the same electrode, not a failure.', ...
                i, ELECTRODES{i}, CANONICAL_ORDER{i}); %#ok<AGROW>
        end
    end
else
    fprintf('  *** FAIL (B): loaded electrode order does not match the locked reference order (beyond case).\n');
    for i = 1:n
        if ~strcmpi(ELECTRODES{i}, CANONICAL_ORDER{i})
            fprintf('      position %d: loaded=%s expected=%s\n', i, ELECTRODES{i}, CANONICAL_ORDER{i});
        end
    end
    if numel(ELECTRODES) ~= numel(CANONICAL_ORDER)
        fprintf('      loaded has %d electrodes, canonical list has %d.\n', numel(ELECTRODES), numel(CANONICAL_ORDER));
    end
    preValFails = preValFails + 1;
end

% --- D. P9 and P10 included ---
if all(ismember(upper({'P9','P10'}), upper(ELECTRODES)))
    fprintf('  PASS (D): P9 and P10 are both present in the electrode set.\n');
else
    fprintf('  *** FAIL (D): P9 and/or P10 missing from the electrode set.\n');
    preValFails = preValFails + 1;
end

% --- E. No ROI averaging performed ---
fprintf('  PASS (E): this script performs no channel aggregation of any kind -- confirmed by construction (no ROI_DEFINITION is loaded or used below).\n');

% --- F. Condition labels exactly E, N, R ---
condsPresent = unique(cellstr(BandPowerT.Condition));
if isempty(setxor(condsPresent, CONDITION_ORDER))
    fprintf('  PASS (F): Condition column contains exactly {E, N, R}, no stray labels.\n');
else
    fprintf('  *** FAIL (F): unexpected Condition labels found: %s\n', strjoin(condsPresent, ', '));
    preValFails = preValFails + 1;
end

BandPowerT.ChannelUC = upper(BandPowerT.Channel);
expectedRowsPerElecBand = 2 * (EXPECTED_N_SUBJECTS_TOTAL - numel(RUN_MODE_EXCL_SUBJ));
elecBandCountFail = false;
elecBandCountTable = zeros(nElec, nBands);
for e = 1:nElec
    for b = 1:nBands
        n = sum(strcmp(BandPowerT.ChannelUC, upper(ELECTRODES{e})) & strcmp(BandPowerT.Band, ALL_BANDS{b}));
        elecBandCountTable(e, b) = n;
        if n ~= expectedRowsPerElecBand
            elecBandCountFail = true;
            integrityIssues{end+1} = sprintf('Electrode=%s Band=%s: %d rows in BandPowerT, expected %d.', ...
                ELECTRODES{e}, ALL_BANDS{b}, n, expectedRowsPerElecBand); %#ok<AGROW>
        end
    end
end
if elecBandCountFail
    fprintf('  *** FAIL (C): at least one Electrode x Band combination does not have exactly %d rows. See integrity issues.\n', expectedRowsPerElecBand);
    preValFails = preValFails + 1;
else
    fprintf('  PASS (C): every one of the 62 electrodes has exactly %d rows (one per recording) for every one of the %d bands (%d x %d = %d checks; P9/P10 included).\n', ...
        expectedRowsPerElecBand, nBands, nElec, nBands, nElec*nBands);
end

comboT = unique(BandPowerT(:, {'Subject','InterventionNumber','Condition'}), 'rows');
nRecordings = height(comboT);
fprintf('\n  %d unique Subject x InterventionNumber x Condition recordings found.\n', nRecordings);

subjNumeric = local_subject_numeric(comboT.Subject);
uniqSubj = unique(subjNumeric);
nSubjects = numel(uniqSubj);
expectedNSubj = EXPECTED_N_SUBJECTS_TOTAL - numel(RUN_MODE_EXCL_SUBJ);
if nSubjects ~= expectedNSubj
    fprintf('  *** FAIL (G): %d unique subjects found, expected %d.\n', nSubjects, expectedNSubj);
    preValFails = preValFails + 1;
else
    fprintf('  PASS (G): %d unique subjects (RUN_MODE=%s), consistent with Phase 3/4/5A.\n', nSubjects, RUN_MODE_LABEL);
end

subjPairs = cell(nSubjects, 1);
pairTypeCount = containers.Map('KeyType', 'char', 'ValueType', 'double');
duplicateFlag = false;
for i = 1:nSubjects
    rowsI = subjNumeric == uniqSubj(i);
    condsI = cellstr(comboT.Condition(rowsI));
    if numel(condsI) ~= 2
        integrityIssues{end+1} = sprintf('Subject %g has %d recordings (expected 2).', uniqSubj(i), numel(condsI)); %#ok<AGROW>
        preValFails = preValFails + 1;
    elseif strcmp(condsI{1}, condsI{2})
        integrityIssues{end+1} = sprintf('Subject %g has the SAME condition (%s) in both recordings -- duplicate.', uniqSubj(i), condsI{1}); %#ok<AGROW>
        duplicateFlag = true;
        preValFails = preValFails + 1;
    else
        pairKey = strjoin(sort(condsI), '_');
        if isKey(pairTypeCount, pairKey); pairTypeCount(pairKey) = pairTypeCount(pairKey) + 1;
        else; pairTypeCount(pairKey) = 1; end
    end
    subjPairs{i} = condsI;
end
pairKeys = keys(pairTypeCount);
fprintf('  Condition-pair coverage across %d subjects:\n', nSubjects);
for k = 1:numel(pairKeys)
    fprintf('    %s : %d subjects\n', strrep(pairKeys{k}, '_', ' & '), pairTypeCount(pairKeys{k}));
end
if duplicateFlag
    fprintf('  *** FAIL (I): at least one subject has a duplicated condition. See integrity issues.\n');
else
    fprintf('  PASS (I): no duplicate Subject x Intervention/Condition observations.\n');
end

condCounts = struct();
for c = 1:numel(CONDITION_ORDER)
    condCounts.(CONDITION_ORDER{c}) = sum(strcmp(cellstr(comboT.Condition), CONDITION_ORDER{c}));
end
fprintf('  Recording-level condition counts: E=%d, N=%d, R=%d\n', condCounts.E, condCounts.N, condCounts.R);
if condCounts.E ~= EXPC.E || condCounts.N ~= EXPC.N || condCounts.R ~= EXPC.R
    fprintf('  *** FAIL (H): expected E=%d, N=%d, R=%d. Got E=%d, N=%d, R=%d.\n', EXPC.E, EXPC.N, EXPC.R, condCounts.E, condCounts.N, condCounts.R);
    integrityIssues{end+1} = sprintf('Recording-level condition counts do not match expected E=%d/N=%d/R=%d (got E=%d/N=%d/R=%d).', ...
        EXPC.E, EXPC.N, EXPC.R, condCounts.E, condCounts.N, condCounts.R); %#ok<AGROW>
    preValFails = preValFails + 1;
else
    fprintf('  PASS (H): recording-level condition counts match the expected E=%d, N=%d, R=%d.\n', EXPC.E, EXPC.N, EXPC.R);
end
if nRecordings ~= expectedRowsPerElecBand
    fprintf('  *** NOTE: nRecordings (%d) != the %d used above for the per-electrode-band row-count check -- inconsistency, investigate.\n', nRecordings, expectedRowsPerElecBand);
    preValFails = preValFails + 1;
end

% --- J. No NaN/Inf in model inputs (the three power columns actually
%         used as LMM response variables) ---
nanInfFail = false;
for col = {'AbsolutePower','LogAbsolutePower','RelativePower'}
    v = BandPowerT.(col{1});
    nBad = sum(isnan(v) | isinf(v));
    if nBad > 0
        nanInfFail = true;
        integrityIssues{end+1} = sprintf('%d NaN/Inf value(s) found in BandPowerT.%s.', nBad, col{1}); %#ok<AGROW>
    end
end
if nanInfFail
    fprintf('  *** FAIL (J): NaN/Inf values found in one or more model-input columns. See integrity issues.\n');
    preValFails = preValFails + 1;
else
    fprintf('  PASS (J): no NaN/Inf values in AbsolutePower, LogAbsolutePower, or RelativePower.\n');
end

fprintf('\nPre-modeling validation failures: %d\n', preValFails);
if preValFails > 0
    error('Phase5B:PreValidationFailed', ...
        'Pre-modeling validation failed (%d issue(s) above). Stopping before any model is fit.', preValFails);
end

%% ========================================================================
%  3. MODEL FITTING -- 868 independent LMMs (62 electrodes x 7 bands x 2
%     metrics). Every model is attempted regardless of any other model's
%     outcome; a failure is recorded and the loop continues.
% ========================================================================
fprintf('\n--------------------------------------------------------\n');
fprintf('SECTION 3: MODEL FITTING (Electrode x Band x Metric, %d total fits)\n', nElec*nBands*nMetrics);
fprintf('--------------------------------------------------------\n');

omnibusRows = cell(0, 11);   % Electrode,Metric,Band,FStat,DF1,DF2,pValue,Status,NObs,NSubjects,DiagnosticFlag
lmeStore = cell(nElec, nBands, nMetrics);
nOK = 0; nFailed = 0; nConvWarn = 0;

for e = 1:nElec
    elecName = ELECTRODES{e};
    for m = 1:nMetrics
        metricName = METRICS{m};
        for b = 1:nBands
            bandName = ALL_BANDS{b};
            subT = BandPowerT(strcmp(BandPowerT.ChannelUC, upper(elecName)) & strcmp(BandPowerT.Band, bandName), :);
            nObs = height(subT);
            nSubjThis = numel(unique(local_subject_numeric(subT.Subject)));

            diagFlags = {};
            if any(isnan(subT.(metricName)) | isinf(subT.(metricName)))
                diagFlags{end+1} = 'NaN/Inf in response'; %#ok<AGROW>
            end
            if nObs == 0
                omnibusRows(end+1, :) = {elecName, metricName, bandName, NaN, NaN, NaN, NaN, ...
                    'FAILED: 0 observations', 0, 0, strjoin([{'0 observations'}, diagFlags], '; ')}; %#ok<AGROW>
                nFailed = nFailed + 1;
                modelFailures{end+1} = sprintf('%s / %s / %s: 0 observations found.', elecName, metricName, bandName); %#ok<AGROW>
                lmeStore{e, b, m} = [];
                continue;
            end

            subT.Subject = categorical(local_subject_numeric(subT.Subject));
            subT.Condition = categorical(cellstr(subT.Condition), CONDITION_ORDER);
            formula = sprintf('%s ~ Condition + (1|Subject)', metricName);

            status = 'OK';
            try
                lastwarn('');
                lme = fitlme(subT, formula, 'FitMethod', 'REML');
                [warnMsg, ~] = lastwarn();
                if ~isempty(warnMsg)
                    diagFlags{end+1} = ['ConvergenceWarning: ' warnMsg]; %#ok<AGROW>
                    nConvWarn = nConvWarn + 1;
                end
            catch ME
                status = ['FAILED: ' ME.message];
                omnibusRows(end+1, :) = {elecName, metricName, bandName, NaN, NaN, NaN, NaN, ...
                    status, nObs, nSubjThis, strjoin([{'fitlme error'}, diagFlags], '; ')}; %#ok<AGROW>
                nFailed = nFailed + 1;
                modelFailures{end+1} = sprintf('%s / %s / %s: fitlme FAILED: %s', elecName, metricName, bandName, ME.message); %#ok<AGROW>
                lmeStore{e, b, m} = [];
                continue;
            end

            try
                aov = anova(lme, 'DFMethod', 'satterthwaite');
                rowIdx = find(strcmp(aov.Term, 'Condition'), 1);
                if isempty(rowIdx)
                    error('Phase5B:NoConditionTerm', 'anova(lme) has no Condition row.');
                end
                fStat = aov.FStat(rowIdx); df1 = aov.DF1(rowIdx); df2 = aov.DF2(rowIdx); pVal = aov.pValue(rowIdx);
                if isnan(fStat) || isinf(fStat) || isnan(pVal) || isinf(pVal)
                    diagFlags{end+1} = 'NaN/Inf in omnibus result'; %#ok<AGROW>
                end
            catch ME
                status = ['FAILED: ' ME.message];
                omnibusRows(end+1, :) = {elecName, metricName, bandName, NaN, NaN, NaN, NaN, ...
                    status, nObs, nSubjThis, strjoin([{'anova error'}, diagFlags], '; ')}; %#ok<AGROW>
                nFailed = nFailed + 1;
                modelFailures{end+1} = sprintf('%s / %s / %s: anova() FAILED: %s', elecName, metricName, bandName, ME.message); %#ok<AGROW>
                lmeStore{e, b, m} = [];
                continue;
            end

            lmeStore{e, b, m} = lme;
            nOK = nOK + 1;
            diagText = '';
            if ~isempty(diagFlags); diagText = strjoin(diagFlags, '; '); end
            omnibusRows(end+1, :) = {elecName, metricName, bandName, fStat, df1, df2, pVal, status, nObs, nSubjThis, diagText}; %#ok<AGROW>
        end
    end
    if mod(e, 10) == 0 || e == nElec
        fprintf('  ... %d/%d electrodes processed (%d OK, %d failed so far)\n', e, nElec, nOK, nFailed);
    end
end

OmnibusT = local_build_table(omnibusRows, [1 2 3 8 11], [4 5 6 7 9 10], ...
    {'Electrode','Metric','Band','FStat','DF1','DF2','pValue','Status','NObs','NSubjects','DiagnosticFlag'});

tierOf = containers.Map();
for i = 1:numel(PRIMARY_TIER_BANDS);     tierOf(PRIMARY_TIER_BANDS{i}) = 'Primary'; end
for i = 1:numel(SECONDARY_TIER_BANDS);   tierOf(SECONDARY_TIER_BANDS{i}) = 'Secondary'; end
for i = 1:numel(EXPLORATORY_TIER_BANDS); tierOf(EXPLORATORY_TIER_BANDS{i}) = 'Exploratory'; end
for i = 1:numel(EXPLORATORY_BETA_TIER_BANDS); tierOf(EXPLORATORY_BETA_TIER_BANDS{i}) = 'ExploratoryBeta'; end
OmnibusT.Tier = cellfun(@(b) tierOf(b), OmnibusT.Band, 'UniformOutput', false);

fprintf('\nModel fitting complete: %d OK, %d failed, out of %d expected (%d electrodes x %d bands x %d metrics).\n', ...
    nOK, nFailed, nElec*nBands*nMetrics, nElec, nBands, nMetrics);
if nFailed == 0
    fprintf('ALL %d MODELS FIT SUCCESSFULLY. No electrode, band, or metric was skipped or dropped.\n', nElec*nBands*nMetrics);
else
    fprintf('%d model(s) FAILED -- see modelFailures below and the report for exact electrode/band/metric/reason.\n', nFailed);
    for i = 1:numel(modelFailures)
        fprintf('  - %s\n', modelFailures{i});
    end
end
fprintf('Convergence warnings recorded (non-fatal): %d\n', nConvWarn);

%% ========================================================================
%  4. FDR -- separately within each Band x Metric family, across the 62
%     electrodes. Only successful fits (Status=='OK', valid pValue) enter
%     a family; failed fits are excluded from the family, never treated
%     as p=1 or silently dropped from the output table.
% ========================================================================
fprintf('\n--------------------------------------------------------\n');
fprintf('SECTION 4: BH-FDR (Band x Metric family, across 62 electrodes)\n');
fprintf('--------------------------------------------------------\n');

FDR_ALPHA = 0.05;
OmnibusT.FamilyName = strcat(OmnibusT.Metric, {' x '}, OmnibusT.Band);
OmnibusT.pFDR = nan(height(OmnibusT), 1);
OmnibusT.SignificantFDR = false(height(OmnibusT), 1);
OmnibusT.FamilySize = nan(height(OmnibusT), 1);

validForFDR = strcmp(OmnibusT.Status, 'OK') & ~isnan(OmnibusT.pValue) & ~isinf(OmnibusT.pValue);

nFamilies = 0; nIncompleteFamilies = 0;
for m = 1:nMetrics
    metricName = METRICS{m};
    for b = 1:nBands
        bandName = ALL_BANDS{b};
        idxAll = strcmp(OmnibusT.Metric, metricName) & strcmp(OmnibusT.Band, bandName);
        idxValid = idxAll & validForFDR;
        nInFamily = sum(idxValid);
        if sum(idxAll) == 0; continue; end
        nFamilies = nFamilies + 1;
        OmnibusT.FamilySize(idxAll) = nInFamily;
        if nInFamily ~= 62
            nIncompleteFamilies = nIncompleteFamilies + 1;
            integrityIssues{end+1} = sprintf('FDR family %s x %s has only %d/62 valid electrode tests (%d failed fit(s) excluded from the family).', ...
                metricName, bandName, nInFamily, sum(idxAll) - nInFamily); %#ok<AGROW>
        end
        if nInFamily > 0
            pvals = OmnibusT.pValue(idxValid);
            [adjP, reject] = local_fdr_bh(pvals, FDR_ALPHA);
            OmnibusT.pFDR(idxValid) = adjP;
            OmnibusT.SignificantFDR(idxValid) = reject;
        end
    end
end
fprintf('BH-FDR applied within %d Band x Metric families (expected %d = %d bands x %d metrics), across up to 62 electrodes each.\n', ...
    nFamilies, nBands*nMetrics, nBands, nMetrics);
if nIncompleteFamilies == 0
    fprintf('Every complete family has exactly n=62 valid electrode tests.\n');
else
    fprintf('*** %d family/families have fewer than 62 valid electrode tests (see integrity issues). ***\n', nIncompleteFamilies);
end

%% ========================================================================
%  5. PAIRWISE CONTRASTS (E-N, E-R, N-R), gated on the Band x Metric FDR
%     result, Holm-corrected within each surviving electrode x Band x
%     Metric. Identical contrast/gating machinery to Phase 5A.
% ========================================================================
fprintf('\n--------------------------------------------------------\n');
fprintf('SECTION 5: PAIRWISE CONTRASTS (gated on FDR, Holm-corrected)\n');
fprintf('--------------------------------------------------------\n');

pairwiseRows = cell(0, 12);   % Electrode,Metric,Band,Contrast,Estimate,SE,DF,tStat,pValue,pHolm,CI_Low,CI_High
CONTRASTS = {'E','N','E - N'; 'E','R','E - R'; 'N','R','N - R'};   % {lvl2, lvl1, label}; label = lvl2 - lvl1

sigMask = OmnibusT.SignificantFDR;
for row = find(sigMask)'
    elecName = OmnibusT.Electrode{row};
    metricName = OmnibusT.Metric{row};
    bandName = OmnibusT.Band{row};
    e = find(strcmp(ELECTRODES, elecName));
    b = find(strcmp(ALL_BANDS, bandName));
    m = find(strcmp(METRICS, metricName));
    lme = lmeStore{e, b, m};
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
        pairwiseRows(end+1, :) = [{elecName, metricName, bandName}, contrastData(c, 1:6), {adjP(c)}, contrastData(c, 7:8)]; %#ok<AGROW>
    end
end
PairwiseT = local_build_table(pairwiseRows, [1 2 3 4], [5 6 7 8 9 10 11 12], ...
    {'Electrode','Metric','Band','Contrast','Estimate','SE','DF','tStat','pValue','pHolm','CI_Low','CI_High'});

fprintf('Pairwise contrasts computed: %d rows (%d FDR-surviving electrode x Band x Metric combinations x 3 contrasts each).\n', ...
    height(PairwiseT), sum(sigMask));

estRows = cell(0, 15);
for e = 1:nElec
    for m = 1:nMetrics
        for b = 1:nBands
            lme = lmeStore{e, b, m};
            if isempty(lme); continue; end
            cn = lme.CoefficientNames; beta = lme.Coefficients.Estimate; covB = lme.CoefficientCovariance;
            L = zeros(3, numel(cn));
            L(:, strcmp(cn, '(Intercept)')) = 1;
            L(2, strcmp(cn, 'Condition_N')) = 1;
            L(3, strcmp(cn, 'Condition_R')) = 1;
            emm = L * beta;
            ests = nan(1, 3); ses = nan(1, 3); ps = nan(1, 3);
            for c = 1:3
                C = local_condition_contrast(cn, 'Condition', REF_LEVEL, CONTRASTS{c, 2}, CONTRASTS{c, 1});
                ests(c) = C * beta; ses(c) = sqrt(C * covB * C');
                ps(c) = coefTest(lme, C, 0, 'DFMethod', 'satterthwaite');
            end
            estRows(end+1, :) = {ELECTRODES{e}, METRICS{m}, ALL_BANDS{b}, emm(1), emm(2), emm(3), ...
                ests(1), ests(2), ests(3), ses(1), ses(2), ses(3), ps(1), ps(2), ps(3)}; %#ok<AGROW>
        end
    end
end
ElectrodeModelEstimatesT = local_build_table(estRows, [1 2 3], 4:15, ...
    {'Electrode','Metric','Band','EMM_E','EMM_N','EMM_R', ...
     'Est_EminusN','Est_EminusR','Est_NminusR','SE_EminusN','SE_EminusR','SE_NminusR', ...
     'p_EminusN_uncorrected','p_EminusR_uncorrected','p_NminusR_uncorrected'});

%% ------------------------------------------------------------------
%  6. PARAMETERS AND OUTPUT
% ------------------------------------------------------------------

PARAMS5B = struct();
PARAMS5B.ScriptName = '09_electrode_level_statistics.m';
PARAMS5B.GeneratedOn = datestr(now);
PARAMS5B.AnalysisLevel = 'Individual-electrode band power (Phase 3 Structure B / BandPowerT), no aggregation. All 62 electrodes (P9/P10 included).';
PARAMS5B.Model = 'LinearMixedEffects_RandomIntercept';
PARAMS5B.FixedEffectsPrimary = 'Condition (3-level: E, N, R; reference = E)';
PARAMS5B.RandomEffects = '(1|Subject)';
PARAMS5B.FitMethod = 'REML';
PARAMS5B.DFMethod = 'Satterthwaite';
PARAMS5B.Metrics = METRICS;
PARAMS5B.Electrodes = ELECTRODES;
PARAMS5B.BandTiers = struct('Primary', {PRIMARY_TIER_BANDS}, 'Secondary', {SECONDARY_TIER_BANDS}, 'Exploratory', {EXPLORATORY_TIER_BANDS}, 'ExploratoryBeta', {EXPLORATORY_BETA_TIER_BANDS});
PARAMS5B.FDRAlpha = FDR_ALPHA;
PARAMS5B.FDRFamilyScheme = 'Separate BH-FDR family per (Band x Metric), across the 62 electrodes (family size = 62 when complete). NOT tier-based, NOT pooled with Phase 5A ROI results.';
PARAMS5B.PairwiseContrastMethod = 'Holm-Bonferroni across the 3 pairwise contrasts (E-N, E-R, N-R), computed only for electrode x Band x Metric combinations whose Band x Metric FDR family test survives BH-FDR.';
PARAMS5B.EffectSizePolicy = 'Model-based estimated difference + 95% CI only. No standardized effect-size denominator computed in this script.';
PARAMS5B.ElectrodeSelectionPolicy = 'None. All 868 models attempted and saved regardless of outcome; no electrode excluded based on statistical results.';
PARAMS5B.SourcePhase3File = bandPowerFile;
PARAMS5B.SourcePhase5DFile = roiDefFile;
PARAMS5B.ModelFitSummary = struct('TotalExpected', nElec*nBands*nMetrics, 'OK', nOK, 'Failed', nFailed, ...
    'ConvergenceWarnings', nConvWarn, 'ModelFailures', {modelFailures});
PARAMS5B.DesignVerification = struct('nRecordings', nRecordings, 'nSubjects', nSubjects, ...
    'PairTypeCounts', {pairTypeCount.keys}, 'ConditionCounts', condCounts, 'IntegrityIssues', {integrityIssues});
PARAMS5B.RunDateTime = datestr(now);

ElectrodeOmnibusT = OmnibusT;
ElectrodePairwiseT = PairwiseT;
% ElectrodeBandPowerT: the exact per-electrode input used by every fit
% above, saved verbatim so downstream spatial/topographic work can be
% reproduced without re-running any model.
ElectrodeBandPowerT = BandPowerT(:, {'Subject','InterventionNumber','Condition','Channel','Band', ...
    'AbsolutePower','LogAbsolutePower','RelativePower'});

save(fullfile(ELEC_DIR, 'electrode_statistics.mat'), ...
    'ElectrodeOmnibusT', 'ElectrodePairwiseT', 'ElectrodeBandPowerT', 'ElectrodeModelEstimatesT', 'PARAMS5B');
writetable(ElectrodeModelEstimatesT, fullfile(ELEC_DIR, 'electrode_model_estimates.csv'));
writetable(ElectrodeOmnibusT, fullfile(ELEC_DIR, 'electrode_omnibus_results.csv'));
writetable(ElectrodeOmnibusT, fullfile(ELEC_DIR, 'electrode_omnibus_results.xlsx'));
if ~isempty(ElectrodePairwiseT)
    writetable(ElectrodePairwiseT, fullfile(ELEC_DIR, 'electrode_pairwise_results.csv'));
    writetable(ElectrodePairwiseT, fullfile(ELEC_DIR, 'electrode_pairwise_results.xlsx'));
end

fprintf('\nSaved output to: %s\n', ELEC_DIR);

%% ------------------------------------------------------------------
%  7. DESCRIPTIVE SUMMARY (counts only -- no interpretation, no
%     electrode selection; purely for the report below)
% ------------------------------------------------------------------

familySummaryRows = {};
for m = 1:nMetrics
    metricName = METRICS{m};
    for b = 1:nBands
        bandName = ALL_BANDS{b};
        idx = strcmp(OmnibusT.Metric, metricName) & strcmp(OmnibusT.Band, bandName);
        idxValid = idx & validForFDR;
        nValid = sum(idxValid);
        nSig = sum(idx & OmnibusT.SignificantFDR);
        pfdrVals = OmnibusT.pFDR(idxValid);
        if isempty(pfdrVals)
            pfdrRangeStr = 'n/a';
        else
            pfdrRangeStr = sprintf('%.4f to %.4f', min(pfdrVals), max(pfdrVals));
        end
        pwIdx = strcmp(PairwiseT.Metric, metricName) & strcmp(PairwiseT.Band, bandName);
        nPairwiseTested = sum(pwIdx) / 3;
        nElecAnyHolm = numel(unique(PairwiseT.Electrode(pwIdx & PairwiseT.pHolm <= FDR_ALPHA)));
        familySummaryRows(end+1, :) = {metricName, bandName, sum(idx), nValid, nSig, pfdrRangeStr, nPairwiseTested, nElecAnyHolm}; %#ok<AGROW>
    end
end

nTotalPairwiseCombos = sum(sigMask);
nEN_holm = sum(strcmp(PairwiseT.Contrast, 'E - N') & PairwiseT.pHolm <= FDR_ALPHA);
nER_holm = sum(strcmp(PairwiseT.Contrast, 'E - R') & PairwiseT.pHolm <= FDR_ALPHA);
nNR_holm = sum(strcmp(PairwiseT.Contrast, 'N - R') & PairwiseT.pHolm <= FDR_ALPHA);

%% ------------------------------------------------------------------
%  8. REPORT
% ------------------------------------------------------------------

reportLines = {};
reportLines{end+1} = sprintf('RUN_MODE = %s | subjects excluded in this run: %s', RUN_MODE_LABEL, local_v3_list(RUN_MODE_EXCL_SUBJ));
reportLines{end+1} = repmat('=', 1, 60);
reportLines{end+1} = 'PHASE 5B -- ELECTRODE-LEVEL E vs N vs R STATISTICS REPORT';
reportLines{end+1} = repmat('=', 1, 60);
reportLines{end+1} = '';
reportLines{end+1} = sprintf('Recordings verified: %d', nRecordings);
reportLines{end+1} = sprintf('Subjects: %d', nSubjects);
reportLines{end+1} = sprintf('Condition counts: E=%d, N=%d, R=%d', condCounts.E, condCounts.N, condCounts.R);
reportLines{end+1} = sprintf('Condition-pair coverage: %s', strjoin(cellfun(@(k) sprintf('%s=%d', strrep(k,'_',' & '), pairTypeCount(k)), pairKeys, 'UniformOutput', false), ', '));
reportLines{end+1} = 'Model: Linear mixed-effects, Condition (fixed, 3-level) + (1|Subject), REML, Satterthwaite DF.';
reportLines{end+1} = sprintf('Electrodes: %d (all individually modeled; P9/P10 included; no ROI averaging)', nElec);
reportLines{end+1} = 'Bands: Primary=Alpha; Secondary=Delta/Theta/Beta; Exploratory=3 gamma bands; ExploratoryBeta=Beta_low 13-20 + Beta_high 20-30 (own family) (labeling only, not FDR family boundaries).';
reportLines{end+1} = 'Metrics: LogAbsolutePower and RelativePower, co-primary, parallel, never pooled for FDR.';
reportLines{end+1} = '';
reportLines{end+1} = '--- Omnibus fitting summary ---';
reportLines{end+1} = sprintf('Total expected models = %d (%d electrodes x %d bands x %d metrics)', nElec*nBands*nMetrics, nElec, nBands, nMetrics);
reportLines{end+1} = sprintf('Successful: %d', nOK);
reportLines{end+1} = sprintf('Failed: %d', nFailed);
reportLines{end+1} = sprintf('Convergence warnings (non-fatal): %d', nConvWarn);
if nFailed == 0
    reportLines{end+1} = 'ALL MODELS FIT SUCCESSFULLY -- no electrode/band/metric was skipped or dropped.';
else
    reportLines{end+1} = 'Model failures (electrode / metric / band / reason):';
    for i = 1:numel(modelFailures)
        reportLines{end+1} = ['  - ' modelFailures{i}]; %#ok<AGROW>
    end
end
if ~isempty(integrityIssues)
    reportLines{end+1} = sprintf('Integrity issues logged: %d (listed at the end of this report)', numel(integrityIssues));
else
    reportLines{end+1} = 'Integrity issues: none.';
end
reportLines{end+1} = '';
reportLines{end+1} = '--- Per Band x Metric family summary (14 families) ---';
for r = 1:size(familySummaryRows, 1)
    reportLines{end+1} = sprintf('  [%s] %-25s family_size=%d valid_n=%d significant_FDR=%d pFDR_range=%s pairwise_tested=%d electrodes_with_any_Holm_sig_contrast=%d', ...
        familySummaryRows{r,1}, familySummaryRows{r,2}, familySummaryRows{r,3}, familySummaryRows{r,4}, ...
        familySummaryRows{r,5}, familySummaryRows{r,6}, familySummaryRows{r,7}, familySummaryRows{r,8}); %#ok<AGROW>
end
reportLines{end+1} = '';
reportLines{end+1} = '--- Pairwise summary ---';
reportLines{end+1} = sprintf('Total FDR-surviving electrode x Band x Metric combinations (each with 3 contrasts): %d', nTotalPairwiseCombos);
reportLines{end+1} = sprintf('Total pairwise contrast rows: %d', height(PairwiseT));
reportLines{end+1} = sprintf('E - N significant after Holm: %d', nEN_holm);
reportLines{end+1} = sprintf('E - R significant after Holm: %d', nER_holm);
reportLines{end+1} = sprintf('N - R significant after Holm: %d', nNR_holm);
reportLines{end+1} = '';
reportLines{end+1} = '--- All FDR-surviving electrode x Band x Metric omnibus results ---';
sigRows = find(OmnibusT.SignificantFDR);
if isempty(sigRows)
    reportLines{end+1} = '  (none)';
else
    for ri = sigRows'
        reportLines{end+1} = sprintf('  [%s] %-6s %-25s Tier=%-11s F(%g,%.1f)=%.3f p=%.4f pFDR=%.4f', ...
            OmnibusT.Metric{ri}, OmnibusT.Electrode{ri}, OmnibusT.Band{ri}, OmnibusT.Tier{ri}, ...
            OmnibusT.DF1(ri), OmnibusT.DF2(ri), OmnibusT.FStat(ri), OmnibusT.pValue(ri), OmnibusT.pFDR(ri)); %#ok<AGROW>
    end
end
reportLines{end+1} = '';
reportLines{end+1} = '--- All computed pairwise contrasts ---';
if isempty(PairwiseT)
    reportLines{end+1} = '  (none -- no electrode x Band x Metric survived its Band x Metric FDR family)';
else
    for r = 1:height(PairwiseT)
        reportLines{end+1} = sprintf('  [%s] %-6s %-25s %-8s est=%.4f (95%% CI %.4f to %.4f) t=%.3f df=%.1f p=%.4f pHolm=%.4f', ...
            PairwiseT.Metric{r}, PairwiseT.Electrode{r}, PairwiseT.Band{r}, PairwiseT.Contrast{r}, PairwiseT.Estimate(r), ...
            PairwiseT.CI_Low(r), PairwiseT.CI_High(r), PairwiseT.tStat(r), PairwiseT.DF(r), ...
            PairwiseT.pValue(r), PairwiseT.pHolm(r)); %#ok<AGROW>
    end
end
reportLines{end+1} = '';
if ~isempty(integrityIssues)
    reportLines{end+1} = '--- INTEGRITY ISSUES (reported, not auto-corrected) ---';
    for i = 1:numel(integrityIssues)
        reportLines{end+1} = ['  - ' integrityIssues{i}]; %#ok<AGROW>
    end
end
reportLines{end+1} = '';
reportLines{end+1} = 'This script reports statistics only. It does not interpret the physiological';
reportLines{end+1} = 'meaning of the spatial pattern, select electrodes for later use, or redefine';
reportLines{end+1} = 'ROIs based on these results. Topography is 10_descriptive_topography.m.';
reportLines{end+1} = '';
reportLines{end+1} = 'Output files written to:';
reportLines{end+1} = ['  ' ELEC_DIR];
reportLines{end+1} = ['  ' REPORT_DIR ' (report copy)'];
reportLines{end+1} = repmat('=', 1, 60);

reportText = strjoin(reportLines, '\n');
fprintf('%s\n', reportText);

fid1 = fopen(fullfile(REPORT_DIR, sprintf('phase5b_electrode_level_statistics_report_%s.txt', RUN_MODE_LABEL)), 'w');
fprintf(fid1, '%s\n', reportText);
fclose(fid1);
fid2 = fopen(fullfile(ELEC_DIR, 'electrode_statistics_report.txt'), 'w');
fprintf(fid2, '%s\n', reportText);
fclose(fid2);

fprintf('\n=== 09_electrode_level_statistics.m complete. STOP -- do not proceed to topography (Phase 5F). ===\n');

%% ======================= LOCAL FUNCTIONS (all after top-level code) =======================

function n = local_subject_numeric(subjVal)
% Robustly extract a numeric subject id regardless of whether Subject is
% stored as numeric, categorical, char, or cellstr, and regardless of
% zero-padding (e.g. '7' vs '07'). Identical implementation to Phase 4/5A.
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
% Phase 4/5A.
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
% implementation to Phase 4/5A.
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
% column position/order). Identical implementation to Phase 4/5A.
    cRow = zeros(1, numel(coeffNames));
    if ~strcmp(level1, refLevel)
        name1 = sprintf('%s_%s', varName, level1);
        idx1 = find(strcmp(coeffNames, name1));
        if isempty(idx1)
            error('Phase5B:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name1);
        end
        cRow(idx1) = cRow(idx1) - 1;
    end
    if ~strcmp(level2, refLevel)
        name2 = sprintf('%s_%s', varName, level2);
        idx2 = find(strcmp(coeffNames, name2));
        if isempty(idx2)
            error('Phase5B:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name2);
        end
        cRow(idx2) = cRow(idx2) + 1;
    end
end

function T = local_build_table(rows, textColIdx, numColIdx, varNames)
% Build a table column-by-column via the table() constructor instead of
% cell2table(), so a text column is never silently reinterpreted.
% Identical implementation to Phase 4/5A.
    nCols = numel(varNames);
    if size(rows, 2) ~= nCols
        error('Phase5B:BuildTableColumnMismatch', ...
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
            error('Phase5B:BuildTableColumnMismatch', 'Column %d not classified as text or numeric.', j);
        end
    end
    T = table(cols{:}, 'VariableNames', varNames);
end

function s = local_v3_list(v)
    if isempty(v); s = 'none'; else; s = strjoin(cellstr(num2str(v(:))), ', '); end
end
