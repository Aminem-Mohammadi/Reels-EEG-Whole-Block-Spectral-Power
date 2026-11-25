%% ========================================================================
% X18_position_int1_int2_intervention.m
%
% Intervention-position model: Metric ~ Condition*Position + (1|Subject), with Position =
% first or second intervention (InterventionNumber 1 or 2).
% Inputs : 01_audit, 04_subject_level.
% Outputs: 18_position_int1_int2/<RUN_MODE>/, 09_reports/.
% Run twice, with RUN_MODE = 'full_sample' (N = 30) and RUN_MODE = 'sensitivity_excl_noisy'
% (cleaned sample, N = 28; excluded subjects listed in EXCLUSION_CONFIG_intervention.m).
% Outputs are written to <OUT_DIR>/<folder>/<RUN_MODE>/ and reports to <OUT_DIR>/09_reports/.
% ========================================================================

clear; clc;

%% ------------------------------------------------------------------
%  SECTION 0: PATHS AND EXISTENCE CHECKS (identical ANALYSIS_ROOT
%  derivation used throughout this entire project -- copied verbatim from
%  X01/X05, never re-derived differently)
% ------------------------------------------------------------------
config;
DATA_ROOT_FOR_SIBLING = DATA_DIR;
[DATA_ROOT_PARENT, ~] = fileparts(DATA_ROOT_FOR_SIBLING);
ANALYSIS_ROOT = OUT_DIR;

AUDIT_DIR  = fullfile(ANALYSIS_ROOT, '01_audit');
SUBJ_DIR   = fullfile(ANALYSIS_ROOT, '04_subject_level');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');

auditFile     = fullfile(AUDIT_DIR, 'audit_table_intervention.csv');
bandPowerFile = fullfile(SUBJ_DIR, 'phase3_subject_level_bandpower.mat');

if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

if ~exist(bandPowerFile, 'file')
    error('Phase18:MissingInput', ...
        'X04 output not found: %s\nRun X04_aggregate_subject_level.m first.', bandPowerFile);
end
if isempty(ver('stats'))
    error('Phase18:MissingToolbox', ...
        'Statistics and Machine Learning Toolbox not detected. This script requires fitlme, anova(LinearMixedModel), coefTest, and qqplot.');
end

addpath(fileparts(mfilename('fullpath')));
run(fullfile(fileparts(mfilename('fullpath')), 'EXCLUSION_CONFIG_intervention.m'));

RUN_MODE = 'full_sample';   % 'full_sample' (N = 30) or 'sensitivity_excl_noisy' (N = 28)
if ~ismember(RUN_MODE, RUN_MODES_intervention)
    error('Phase18:BadRunMode', 'RUN_MODE must be one of: %s', strjoin(RUN_MODES_intervention, ', '));
end
RUN_MODE_LABEL = RUN_MODES_DISPLAY_intervention.(RUN_MODE);
fprintf('RUN_MODE = %s\n', RUN_MODE_LABEL);
if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
    fprintf('Excluding subject(s) {%s} from every model in this run -- see EXCLUSION_CONFIG_intervention.m for rationale.\n', ...
        strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','));
end
% ---- END NEW ----

fprintf('========================================================\n');
fprintf('PHASE 18 (X18) -- POSITION (FIRST vs SECOND) ORDER/FATIGUE EFFECT\n');
fprintf('Exploratory. Position = Int1 (First) vs Int2 (Second).\n');
fprintf('========================================================\n\n');

%% ------------------------------------------------------------------
S4 = intervention_load(bandPowerFile, 'BandPowerT', 'GlobalBandPowerT', 'PARAMS3');
BandPowerT_ALL = S4.BandPowerT;
GlobalBandPowerT_ALL = S4.GlobalBandPowerT;
PARAMS3 = S4.PARAMS3;
clear S4;

requiredGlobalCols = {'Subject','TaskType','BlockNumber','Condition','Band', ...
    'LogAbsolutePower','RelativePower'};
missingG = requiredGlobalCols(~ismember(requiredGlobalCols, GlobalBandPowerT_ALL.Properties.VariableNames));
if ~isempty(missingG)
    error('Phase18:MissingColumns', 'GlobalBandPowerT is missing expected column(s): %s', strjoin(missingG, ', '));
end
requiredElecCols = {'Subject','TaskType','BlockNumber','Condition','Channel','Band', ...
    'LogAbsolutePower','RelativePower'};
missingE = requiredElecCols(~ismember(requiredElecCols, BandPowerT_ALL.Properties.VariableNames));
if ~isempty(missingE)
    error('Phase18:MissingColumns', 'BandPowerT is missing expected column(s): %s', strjoin(missingE, ', '));
end

ALL_BANDS = fieldnames(PARAMS3.BandDefinitionsHz);
METRICS   = {'LogAbsolutePower', 'RelativePower'};
CONDITION_ORDER = {'E', 'N', 'R'};
REF_LEVEL = CONDITION_ORDER{1};
POSITION_ORDER = {'First', 'Second'};   % reference = First (judgment call #1, see header)
POSITION_REF = POSITION_ORDER{1};
POSITION_SECOND = POSITION_ORDER{2};

PRIMARY_TIER_BANDS     = {'Alpha'};
SECONDARY_TIER_BANDS   = {'Delta', 'Theta', 'Beta'};
EXPLORATORY_TIER_BANDS = {'Gamma_low_exploratory', 'Gamma_high1_exploratory', 'Gamma_high2_exploratory'};
EXPLORATORY_BETA_TIER_BANDS = {'Beta_low_exploratory', 'Beta_high_exploratory'};
expectedTierBands = [PRIMARY_TIER_BANDS, SECONDARY_TIER_BANDS, EXPLORATORY_TIER_BANDS, EXPLORATORY_BETA_TIER_BANDS];
if ~isempty(setxor(sort(expectedTierBands), sort(ALL_BANDS')))
    error('Phase18:BandMismatch', ...
        ['The hard-coded band tiers do not exactly match PARAMS3.BandDefinitionsHz. ' ...
         'Tiers: %s | Phase3 bands: %s'], strjoin(expectedTierBands, ','), strjoin(ALL_BANDS', ','));
end
tierDefs = {'Primary', PRIMARY_TIER_BANDS; 'Secondary', SECONDARY_TIER_BANDS; 'Exploratory', EXPLORATORY_TIER_BANDS; 'ExploratoryBeta', EXPLORATORY_BETA_TIER_BANDS};
tierOf = containers.Map();
for i = 1:numel(PRIMARY_TIER_BANDS);     tierOf(PRIMARY_TIER_BANDS{i}) = 'Primary'; end     %#ok<NOCOL>
for i = 1:numel(SECONDARY_TIER_BANDS);   tierOf(SECONDARY_TIER_BANDS{i}) = 'Secondary'; end %#ok<NOCOL>
for i = 1:numel(EXPLORATORY_TIER_BANDS); tierOf(EXPLORATORY_TIER_BANDS{i}) = 'Exploratory'; end %#ok<NOCOL>
for i = 1:numel(EXPLORATORY_BETA_TIER_BANDS); tierOf(EXPLORATORY_BETA_TIER_BANDS{i}) = 'ExploratoryBeta'; end

TERMS = {'Condition', 'Position', 'Condition:Position'};
FDR_ALPHA = 0.05;
nBands = numel(ALL_BANDS);
nMetrics = numel(METRICS);

% Electrode reference list -- taken directly from X04's BandPowerT.
ELECTRODES = unique(BandPowerT_ALL.Channel);
nElectrodes = numel(ELECTRODES);

TASK_TYPES = {'intervention'};
TASK_DISPLAY = struct('intervention', 'Intervention');

%% ------------------------------------------------------------------

fprintf('========================================================\n');
fprintf('Single intervention analysis: Position = Int1 (First) vs Int2 (Second).\n');
fprintf('========================================================\n\n');

%% ========================================================================
for taskI = 1:numel(TASK_TYPES)
    TASK_KEY = TASK_TYPES{taskI};
    TASK_LABEL = TASK_DISPLAY.(TASK_KEY);

    STATS_DIR = fullfile(ANALYSIS_ROOT, '18_position_int1_int2', RUN_MODE);
    if ~exist(STATS_DIR, 'dir'); mkdir(STATS_DIR); end

    fprintf('\n########################################################\n');
    fprintf('### TASK = %s | RUN_MODE = %s\n', TASK_LABEL, RUN_MODE_LABEL);
    fprintf('########################################################\n\n');

    integrityIssues = {};

    %% ====================================================================
    thisTaskGlobal = GlobalBandPowerT_ALL(strcmp(GlobalBandPowerT_ALL.TaskType, TASK_KEY), :);
    nBaselineExcludedGlobal = sum(strcmp(thisTaskGlobal.Condition, 'Baseline'));
    GlobalBandPowerT = thisTaskGlobal(ismember(thisTaskGlobal.Condition, CONDITION_ORDER), :);

    thisTaskElec = BandPowerT_ALL(strcmp(BandPowerT_ALL.TaskType, TASK_KEY), :);
    nBaselineExcludedElec = sum(strcmp(thisTaskElec.Condition, 'Baseline'));
    BandPowerT = thisTaskElec(ismember(thisTaskElec.Condition, CONDITION_ORDER), :);

    fprintf('Task=%s: Global rows %d total, %d Baseline excluded, %d E/N/R retained.\n', ...
        TASK_LABEL, height(thisTaskGlobal), nBaselineExcludedGlobal, height(GlobalBandPowerT));
    fprintf('Task=%s: Electrode rows %d total, %d Baseline excluded, %d E/N/R retained.\n', ...
        TASK_LABEL, height(thisTaskElec), nBaselineExcludedElec, height(BandPowerT));

    if isempty(GlobalBandPowerT) || isempty(BandPowerT)
        warning('Task=%s: no E/N/R rows found after filtering -- skipping this Task entirely.', TASK_LABEL);
        continue;
    end

    %% ====================================================================
    recCombos = unique([GlobalBandPowerT(:, {'Subject','BlockNumber'}); BandPowerT(:, {'Subject','BlockNumber'})], 'rows');
    nRecCombos = height(recCombos);

    posLabelByCombo = cell(nRecCombos, 1);
    seqByCombo = cell(nRecCombos, 1);
    nAmbiguous = 0; nMissingAudit = 0; nUnparseable = 0;   % kept for the report; audit is not used

    for i = 1:nRecCombos
        subjStr = recCombos.Subject{i};
        blockNum = recCombos.BlockNumber(i);
        seqByCombo{i} = sprintf('Int%g', blockNum);
        if blockNum == 1
            posLabelByCombo{i} = 'First';
        elseif blockNum == 2
            posLabelByCombo{i} = 'Second';
        else
            posLabelByCombo{i} = '';
            nUnparseable = nUnparseable + 1;
            integrityIssues{end+1} = sprintf('Sub%s: InterventionNumber=%g is not 1 or 2 -- Position NOT guessed; recording EXCLUDED from modeling.', ...
                subjStr, blockNum); %#ok<AGROW>
        end
    end

    fprintf('Position derivation (Task=%s): %d recordings, %d missing-audit-match, %d ambiguous-match, %d unparseable SEQ -- all EXCLUDED, never guessed.\n', ...
        TASK_LABEL, nRecCombos, nMissingAudit, nAmbiguous, nUnparseable);

    % Attach Position to both tables (default '' -> excluded below)
    GlobalBandPowerT.Position = repmat({''}, height(GlobalBandPowerT), 1);
    BandPowerT.Position = repmat({''}, height(BandPowerT), 1);
    PositionByRecRows = cell(nRecCombos, 4);
    for i = 1:nRecCombos
        subjStr = recCombos.Subject{i};
        blockNum = recCombos.BlockNumber(i);
        posLabel = posLabelByCombo{i};
        gMask = strcmp(GlobalBandPowerT.Subject, subjStr) & GlobalBandPowerT.BlockNumber == blockNum;
        eMask = strcmp(BandPowerT.Subject, subjStr) & BandPowerT.BlockNumber == blockNum;
        if ~isempty(posLabel)
            GlobalBandPowerT.Position(gMask) = {posLabel};
            BandPowerT.Position(eMask) = {posLabel};
        end
        condsHere = unique(GlobalBandPowerT.Condition(gMask));
        condStr = strjoin(condsHere, ',');
        PositionByRecRows(i, :) = {subjStr, blockNum, condStr, iif_str(isempty(posLabel), 'EXCLUDED', posLabel)};
    end
    PositionByRecT = cell2table(PositionByRecRows, 'VariableNames', {'Subject','BlockNumber','Condition','Position'});
    writetable(PositionByRecT, fullfile(STATS_DIR, 'phase18_position_by_recording.csv'));

    nExcludedGlobal = sum(strcmp(GlobalBandPowerT.Position, ''));
    nExcludedElec   = sum(strcmp(BandPowerT.Position, ''));
    GlobalBandPowerT = GlobalBandPowerT(~strcmp(GlobalBandPowerT.Position, ''), :);
    BandPowerT       = BandPowerT(~strcmp(BandPowerT.Position, ''), :);
    fprintf('After excluding un-Position-labeled rows: Global %d rows retained (%d excluded); Electrode %d rows retained (%d excluded).\n\n', ...
        height(GlobalBandPowerT), nExcludedGlobal, height(BandPowerT), nExcludedElec);

    if isempty(GlobalBandPowerT) || isempty(BandPowerT)
        warning('Task=%s: no rows remain after Position derivation -- skipping this Task entirely.', TASK_LABEL);
        continue;
    end

    nExcludedRowsSensitivityGlobal = 0;
    nExcludedRowsSensitivityElec = 0;
    if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
        subjNumForExclGlobal = local_subject_numeric(GlobalBandPowerT.Subject);
        excludeMaskGlobal = ismember(subjNumForExclGlobal, EXCLUDE_NOISY_SUBJECTS);
        nExcludedRowsSensitivityGlobal = sum(excludeMaskGlobal);
        GlobalBandPowerT = GlobalBandPowerT(~excludeMaskGlobal, :);

        subjNumForExclElec = local_subject_numeric(BandPowerT.Subject);
        excludeMaskElec = ismember(subjNumForExclElec, EXCLUDE_NOISY_SUBJECTS);
        nExcludedRowsSensitivityElec = sum(excludeMaskElec);
        BandPowerT = BandPowerT(~excludeMaskElec, :);

        fprintf('Task=%s: RUN_MODE=%s -- excluded %d Global-level row(s) and %d Electrode-level row(s) belonging to subject(s) {%s}; %d Global-level rows remain.\n', ...
            TASK_LABEL, RUN_MODE_LABEL, nExcludedRowsSensitivityGlobal, nExcludedRowsSensitivityElec, ...
            strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','), height(GlobalBandPowerT));

        if isempty(GlobalBandPowerT) || isempty(BandPowerT)
            warning('Task=%s: no rows remain after the sensitivity exclusion -- skipping this Task entirely.', TASK_LABEL);
            continue;
        end
    end
    % ---- END NEW ----

    %% ====================================================================
    comboT = unique(GlobalBandPowerT(:, {'Subject','BlockNumber','Condition','Position'}), 'rows');
    nRecordings = height(comboT);
    fprintf('Section 3c: %d unique Subject x BlockNumber x Condition x Position combinations retained (Task=%s).\n', nRecordings, TASK_LABEL);

    posCounts = struct();
    for p = 1:numel(POSITION_ORDER)
        posCounts.(POSITION_ORDER{p}) = sum(strcmp(comboT.Position, POSITION_ORDER{p}));
    end
    fprintf('Position counts (Task=%s): First=%d, Second=%d.\n', TASK_LABEL, posCounts.First, posCounts.Second);

    condPosConfound = false;
    fprintf('Condition x Position cross-tabulation (Task=%s):\n', TASK_LABEL);
    for c = 1:numel(CONDITION_ORDER)
        condName = CONDITION_ORDER{c};
        nFirst = sum(strcmp(comboT.Condition, condName) & strcmp(comboT.Position, 'First'));
        nSecond = sum(strcmp(comboT.Condition, condName) & strcmp(comboT.Position, 'Second'));
        fprintf('  Condition %s: First=%d, Second=%d\n', condName, nFirst, nSecond);
        if nFirst == 0 || nSecond == 0
            condPosConfound = true;
            integrityIssues{end+1} = sprintf(['Task=%s: Condition %s has only one Position level represented ' ...
                '(First=%d, Second=%d) -- the Condition:Position interaction and Condition-marginal contrasts ' ...
                'for this level may be inestimable or unstable; not silently patched, flagged for review.'], ...
                TASK_LABEL, condName, nFirst, nSecond); %#ok<AGROW>
        end
    end

    subjNumeric = local_subject_numeric(comboT.Subject);
    nSubjects = numel(unique(subjNumeric));
    fprintf('Subjects retained (Task=%s): %d.\n', TASK_LABEL, nSubjects);

    if ~isempty(integrityIssues)
        fprintf(2, '\n*** DESIGN VERIFICATION FLAGGED %d ISSUE(S) SO FAR (Task=%s) -- see report. ***\n\n', numel(integrityIssues), TASK_LABEL);
    end

    %% ====================================================================
    %  SECTION 3d: GLOBAL-LEVEL OMNIBUS FITS (Band x Metric; each fit
    %  yields THREE terms: Condition, Position, Condition:Position)
    % ====================================================================
    omnibusGlobalRows = cell(0, 8);
    lmeMapGlobal = containers.Map();   % key: 'Band|Metric'

    for m = 1:nMetrics
        metricName = METRICS{m};
        for b = 1:nBands
            bandName = ALL_BANDS{b};
            subT = GlobalBandPowerT(strcmp(GlobalBandPowerT.Band, bandName), :);
            fr = local_fit_condition_position_omnibus(subT, metricName, CONDITION_ORDER, POSITION_ORDER);
            if ~fr.ok
                integrityIssues{end+1} = sprintf('Global: fitlme/anova FAILED for %s / %s (Task=%s): %s', ...
                    metricName, bandName, TASK_LABEL, fr.errMsg); %#ok<AGROW>
                continue;
            end
            lmeMapGlobal(sprintf('%s|%s', bandName, metricName)) = fr.lme;
            for t = 1:3
                omnibusGlobalRows(end+1, :) = {fr.rows(t).Term, metricName, bandName, tierOf(bandName), ...
                    fr.rows(t).FStat, fr.rows(t).DF1, fr.rows(t).DF2, fr.rows(t).pValue}; %#ok<AGROW>
            end
        end
    end
    OmnibusGlobalT = local_build_table(omnibusGlobalRows, [1 2 3 4], [5 6 7 8], ...
        {'Term','Metric','Band','Tier','FStat','DF1','DF2','pValue'});
    fprintf('Global-level omnibus fits complete (Task=%s): %d rows (%d Band x Metric fits x 3 terms, minus failures).\n', ...
        TASK_LABEL, height(OmnibusGlobalT), nBands*nMetrics);

    %% ====================================================================
    %  SECTION 3e: GLOBAL-LEVEL FDR -- THREE SEPARATE FAMILIES (one per
    %  Term), each using the SAME tiered (Tier x Metric) family structure
    %  as X05's own Global scheme.
    % ====================================================================
    OmnibusGlobalT.pFDR = nan(height(OmnibusGlobalT), 1);
    OmnibusGlobalT.SignificantFDR = false(height(OmnibusGlobalT), 1);
    OmnibusGlobalT.FamilySize = nan(height(OmnibusGlobalT), 1);

    for ti = 1:numel(TERMS)
        termName = TERMS{ti};
        for t = 1:size(tierDefs, 1)
            tierBands = tierDefs{t, 2};
            for m = 1:nMetrics
                metricName = METRICS{m};
                idx = strcmp(OmnibusGlobalT.Term, termName) & strcmp(OmnibusGlobalT.Metric, metricName) & ismember(OmnibusGlobalT.Band, tierBands);
                if ~any(idx); continue; end
                OmnibusGlobalT.FamilySize(idx) = sum(idx);
                pvals = OmnibusGlobalT.pValue(idx);
                if numel(pvals) == 1
                    OmnibusGlobalT.pFDR(idx) = pvals;
                    OmnibusGlobalT.SignificantFDR(idx) = pvals <= FDR_ALPHA;
                else
                    [adjP, reject] = local_fdr_bh(pvals, FDR_ALPHA);
                    OmnibusGlobalT.pFDR(idx) = adjP;
                    OmnibusGlobalT.SignificantFDR(idx) = reject;
                end
            end
        end
    end
    fprintf('Global-level FDR complete (Task=%s), 3 independent Term families: Condition sig=%d, Position sig=%d, Interaction sig=%d (out of %d each).\n', ...
        TASK_LABEL, sum(strcmp(OmnibusGlobalT.Term,'Condition') & OmnibusGlobalT.SignificantFDR), ...
        sum(strcmp(OmnibusGlobalT.Term,'Position') & OmnibusGlobalT.SignificantFDR), ...
        sum(strcmp(OmnibusGlobalT.Term,'Condition:Position') & OmnibusGlobalT.SignificantFDR), nBands*nMetrics);

    %% ====================================================================
    %  SECTION 3f: ELECTRODE-LEVEL OMNIBUS FITS (Channel x Band x Metric)
    % ====================================================================
    omnibusElecRows = cell(0, 8);
    lmeMapElec = containers.Map();   % key: 'Electrode|Band|Metric'
    nOKElec = 0; nFailedElec = 0;

    for m = 1:nMetrics
        metricName = METRICS{m};
        for b = 1:nBands
            bandName = ALL_BANDS{b};
            for e = 1:nElectrodes
                elecName = ELECTRODES{e};
                subT = BandPowerT(strcmp(BandPowerT.Channel, elecName) & strcmp(BandPowerT.Band, bandName), :);
                fr = local_fit_condition_position_omnibus(subT, metricName, CONDITION_ORDER, POSITION_ORDER);
                if ~fr.ok
                    integrityIssues{end+1} = sprintf('Electrode %s: fitlme/anova FAILED for %s / %s (Task=%s): %s', ...
                        elecName, metricName, bandName, TASK_LABEL, fr.errMsg); %#ok<AGROW>
                    nFailedElec = nFailedElec + 1;
                    continue;
                end
                nOKElec = nOKElec + 1;
                lmeMapElec(sprintf('%s|%s|%s', elecName, bandName, metricName)) = fr.lme;
                for t = 1:3
                    omnibusElecRows(end+1, :) = {fr.rows(t).Term, elecName, metricName, bandName, ...
                        fr.rows(t).FStat, fr.rows(t).DF1, fr.rows(t).DF2, fr.rows(t).pValue}; %#ok<AGROW>
                end
            end
        end
        fprintf('  ... metric %s done (Task=%s): %d electrode-fits OK so far, %d failed.\n', metricName, TASK_LABEL, nOKElec, nFailedElec);
    end
    OmnibusElectrodeT = local_build_table(omnibusElecRows, [1 2 3 4], [5 6 7 8], ...
        {'Term','Electrode','Metric','Band','FStat','DF1','DF2','pValue'});
    fprintf('Electrode-level omnibus fits complete (Task=%s): %d OK, %d failed, out of %d expected (%d electrodes x %d bands x %d metrics).\n', ...
        TASK_LABEL, nOKElec, nFailedElec, nElectrodes*nBands*nMetrics, nElectrodes, nBands, nMetrics);

    %% ====================================================================
    %  SECTION 3g: ELECTRODE-LEVEL FDR -- THREE SEPARATE FAMILIES (one per
    %  Term), each per (Band x Metric) across up to 62 electrodes, matching
    %  X09's own Electrode scheme.
    % ====================================================================
    OmnibusElectrodeT.pFDR = nan(height(OmnibusElectrodeT), 1);
    OmnibusElectrodeT.SignificantFDR = false(height(OmnibusElectrodeT), 1);
    OmnibusElectrodeT.FamilySize = nan(height(OmnibusElectrodeT), 1);

    for ti = 1:numel(TERMS)
        termName = TERMS{ti};
        for m = 1:nMetrics
            metricName = METRICS{m};
            for b = 1:nBands
                bandName = ALL_BANDS{b};
                idx = strcmp(OmnibusElectrodeT.Term, termName) & strcmp(OmnibusElectrodeT.Metric, metricName) & strcmp(OmnibusElectrodeT.Band, bandName);
                nInFamily = sum(idx);
                if nInFamily == 0; continue; end
                OmnibusElectrodeT.FamilySize(idx) = nInFamily;
                if nInFamily ~= nElectrodes
                    integrityIssues{end+1} = sprintf('Task=%s: Electrode FDR family %s / %s / %s has only %d/%d valid electrode tests.', ...
                        TASK_LABEL, termName, metricName, bandName, nInFamily, nElectrodes); %#ok<AGROW>
                end
                pvals = OmnibusElectrodeT.pValue(idx);
                if numel(pvals) == 1
                    OmnibusElectrodeT.pFDR(idx) = pvals;
                    OmnibusElectrodeT.SignificantFDR(idx) = pvals <= FDR_ALPHA;
                else
                    [adjP, reject] = local_fdr_bh(pvals, FDR_ALPHA);
                    OmnibusElectrodeT.pFDR(idx) = adjP;
                    OmnibusElectrodeT.SignificantFDR(idx) = reject;
                end
            end
        end
    end
    fprintf('Electrode-level FDR complete (Task=%s), 3 independent Term families: Condition sig=%d, Position sig=%d, Interaction sig=%d.\n', ...
        TASK_LABEL, sum(strcmp(OmnibusElectrodeT.Term,'Condition') & OmnibusElectrodeT.SignificantFDR), ...
        sum(strcmp(OmnibusElectrodeT.Term,'Position') & OmnibusElectrodeT.SignificantFDR), ...
        sum(strcmp(OmnibusElectrodeT.Term,'Condition:Position') & OmnibusElectrodeT.SignificantFDR));

    %% ====================================================================
    %  SECTION 3h: PAIRWISE / SIMPLE-EFFECTS CONTRASTS -- GLOBAL LEVEL
    %  (three independent gating types)
    % ====================================================================
    CONTRASTS = {'E','N','E - N'; 'E','R','E - R'; 'N','R','N - R'};   % {level2, level1, label}; estimate = level2 - level1

    % --- (i) Condition main effect survives -> 3 MARGINAL contrasts, Holm(3) ---
    pwCondGlobalRows = cell(0, 11);
    sigCond = strcmp(OmnibusGlobalT.Term, 'Condition') & OmnibusGlobalT.SignificantFDR;
    for r = find(sigCond)'
        metricName = OmnibusGlobalT.Metric{r}; bandName = OmnibusGlobalT.Band{r};
        key = sprintf('%s|%s', bandName, metricName);
        if ~isKey(lmeMapGlobal, key); continue; end
        lme = lmeMapGlobal(key);
        coeffNames = lme.CoefficientNames; beta = lme.Coefficients.Estimate; covB = lme.CoefficientCovariance;
        contrastP = nan(3, 1); contrastData = cell(3, 8);
        for c = 1:3
            lvl2 = CONTRASTS{c,1}; lvl1 = CONTRASTS{c,2}; label = CONTRASTS{c,3};
            C = local_condition_position_contrast(coeffNames, REF_LEVEL, lvl1, lvl2, POSITION_SECOND, 0.5, 0.5);
            [est, se, tstat, p, df2c, ciLow, ciHigh] = local_apply_contrast(lme, C, beta, covB);
            contrastP(c) = p;
            contrastData(c, :) = {label, est, se, df2c, tstat, p, ciLow, ciHigh};
        end
        [adjP, ~] = local_holm_bonferroni(contrastP, FDR_ALPHA);
        for c = 1:3
            pwCondGlobalRows(end+1, :) = [{metricName, bandName}, contrastData(c, 1:6), {adjP(c)}, contrastData(c, 7:8)]; %#ok<AGROW>
        end
    end
    PairwiseConditionMarginalGlobalT = local_build_table(pwCondGlobalRows, [1 2 3], [4 5 6 7 8 9 10 11], ...
        {'Metric','Band','Contrast','Estimate','SE','DF','tStat','pValue','pHolm','CI_Low','CI_High'});

    % --- (ii) Position main effect survives -> 1 marginal contrast, no Holm ---
    pwPosGlobalRows = cell(0, 11);
    sigPos = strcmp(OmnibusGlobalT.Term, 'Position') & OmnibusGlobalT.SignificantFDR;
    for r = find(sigPos)'
        metricName = OmnibusGlobalT.Metric{r}; bandName = OmnibusGlobalT.Band{r};
        key = sprintf('%s|%s', bandName, metricName);
        if ~isKey(lmeMapGlobal, key); continue; end
        lme = lmeMapGlobal(key);
        coeffNames = lme.CoefficientNames; beta = lme.Coefficients.Estimate; covB = lme.CoefficientCovariance;
        C = local_position_marginal_contrast(coeffNames, CONDITION_ORDER, POSITION_SECOND);
        [est, se, tstat, p, df2c, ciLow, ciHigh] = local_apply_contrast(lme, C, beta, covB);
        pwPosGlobalRows(end+1, :) = {metricName, bandName, 'Second - First', est, se, df2c, tstat, p, p, ciLow, ciHigh}; %#ok<AGROW>
    end
    PairwisePositionMarginalGlobalT = local_build_table(pwPosGlobalRows, [1 2 3], [4 5 6 7 8 9 10 11], ...
        {'Metric','Band','Contrast','Estimate','SE','DF','tStat','pValue','pHolm','CI_Low','CI_High'});

    % --- (iii) Interaction survives -> 6 simple-effects contrasts, Holm(6) ---
    pwIntGlobalRows = cell(0, 12);
    sigInt = strcmp(OmnibusGlobalT.Term, 'Condition:Position') & OmnibusGlobalT.SignificantFDR;
    for r = find(sigInt)'
        metricName = OmnibusGlobalT.Metric{r}; bandName = OmnibusGlobalT.Band{r};
        key = sprintf('%s|%s', bandName, metricName);
        if ~isKey(lmeMapGlobal, key); continue; end
        lme = lmeMapGlobal(key);
        coeffNames = lme.CoefficientNames; beta = lme.Coefficients.Estimate; covB = lme.CoefficientCovariance;
        contrastP = nan(6, 1); contrastData = cell(6, 9);
        rowIdx = 0;
        for p = 1:numel(POSITION_ORDER)
            posLvl = POSITION_ORDER{p};
            w = 0; if strcmp(posLvl, POSITION_SECOND); w = 1; end
            for c = 1:3
                lvl2 = CONTRASTS{c,1}; lvl1 = CONTRASTS{c,2}; label = CONTRASTS{c,3};
                rowIdx = rowIdx + 1;
                C = local_condition_position_contrast(coeffNames, REF_LEVEL, lvl1, lvl2, POSITION_SECOND, w, w);
                [est, se, tstat, p2, df2c, ciLow, ciHigh] = local_apply_contrast(lme, C, beta, covB);
                contrastP(rowIdx) = p2;
                contrastData(rowIdx, :) = {posLvl, label, est, se, df2c, tstat, p2, ciLow, ciHigh};
            end
        end
        [adjP, ~] = local_holm_bonferroni(contrastP, FDR_ALPHA);
        for rowIdx = 1:6
            pwIntGlobalRows(end+1, :) = [{metricName, bandName}, contrastData(rowIdx, 1:7), {adjP(rowIdx)}, contrastData(rowIdx, 8:9)]; %#ok<AGROW>
        end
    end
    PairwiseInteractionSimpleGlobalT = local_build_table(pwIntGlobalRows, [1 2 3 4], [5 6 7 8 9 10 11 12], ...
        {'Metric','Band','Position','Contrast','Estimate','SE','DF','tStat','pValue','pHolm','CI_Low','CI_High'});

    fprintf('Global pairwise/simple-effects (Task=%s): Condition-marginal=%d rows, Position-marginal=%d rows, Interaction-simple=%d rows.\n', ...
        TASK_LABEL, height(PairwiseConditionMarginalGlobalT), height(PairwisePositionMarginalGlobalT), height(PairwiseInteractionSimpleGlobalT));

    %% ====================================================================
    %  SECTION 3i: PAIRWISE / SIMPLE-EFFECTS CONTRASTS -- ELECTRODE LEVEL
    %  (identical logic to Section 3h, keyed by Electrode as well)
    % ====================================================================
    pwCondElecRows = cell(0, 12);
    sigCondE = strcmp(OmnibusElectrodeT.Term, 'Condition') & OmnibusElectrodeT.SignificantFDR;
    for r = find(sigCondE)'
        elecName = OmnibusElectrodeT.Electrode{r}; metricName = OmnibusElectrodeT.Metric{r}; bandName = OmnibusElectrodeT.Band{r};
        key = sprintf('%s|%s|%s', elecName, bandName, metricName);
        if ~isKey(lmeMapElec, key); continue; end
        lme = lmeMapElec(key);
        coeffNames = lme.CoefficientNames; beta = lme.Coefficients.Estimate; covB = lme.CoefficientCovariance;
        contrastP = nan(3, 1); contrastData = cell(3, 8);
        for c = 1:3
            lvl2 = CONTRASTS{c,1}; lvl1 = CONTRASTS{c,2}; label = CONTRASTS{c,3};
            C = local_condition_position_contrast(coeffNames, REF_LEVEL, lvl1, lvl2, POSITION_SECOND, 0.5, 0.5);
            [est, se, tstat, p, df2c, ciLow, ciHigh] = local_apply_contrast(lme, C, beta, covB);
            contrastP(c) = p;
            contrastData(c, :) = {label, est, se, df2c, tstat, p, ciLow, ciHigh};
        end
        [adjP, ~] = local_holm_bonferroni(contrastP, FDR_ALPHA);
        for c = 1:3
            pwCondElecRows(end+1, :) = [{elecName, metricName, bandName}, contrastData(c, 1:6), {adjP(c)}, contrastData(c, 7:8)]; %#ok<AGROW>
        end
    end
    PairwiseConditionMarginalElectrodeT = local_build_table(pwCondElecRows, [1 2 3 4], [5 6 7 8 9 10 11 12], ...
        {'Electrode','Metric','Band','Contrast','Estimate','SE','DF','tStat','pValue','pHolm','CI_Low','CI_High'});

    pwPosElecRows = cell(0, 12);
    sigPosE = strcmp(OmnibusElectrodeT.Term, 'Position') & OmnibusElectrodeT.SignificantFDR;
    for r = find(sigPosE)'
        elecName = OmnibusElectrodeT.Electrode{r}; metricName = OmnibusElectrodeT.Metric{r}; bandName = OmnibusElectrodeT.Band{r};
        key = sprintf('%s|%s|%s', elecName, bandName, metricName);
        if ~isKey(lmeMapElec, key); continue; end
        lme = lmeMapElec(key);
        coeffNames = lme.CoefficientNames; beta = lme.Coefficients.Estimate; covB = lme.CoefficientCovariance;
        C = local_position_marginal_contrast(coeffNames, CONDITION_ORDER, POSITION_SECOND);
        [est, se, tstat, p, df2c, ciLow, ciHigh] = local_apply_contrast(lme, C, beta, covB);
        pwPosElecRows(end+1, :) = {elecName, metricName, bandName, 'Second - First', est, se, df2c, tstat, p, p, ciLow, ciHigh}; %#ok<AGROW>
    end
    PairwisePositionMarginalElectrodeT = local_build_table(pwPosElecRows, [1 2 3 4], [5 6 7 8 9 10 11 12], ...
        {'Electrode','Metric','Band','Contrast','Estimate','SE','DF','tStat','pValue','pHolm','CI_Low','CI_High'});

    pwIntElecRows = cell(0, 13);
    sigIntE = strcmp(OmnibusElectrodeT.Term, 'Condition:Position') & OmnibusElectrodeT.SignificantFDR;
    for r = find(sigIntE)'
        elecName = OmnibusElectrodeT.Electrode{r}; metricName = OmnibusElectrodeT.Metric{r}; bandName = OmnibusElectrodeT.Band{r};
        key = sprintf('%s|%s|%s', elecName, bandName, metricName);
        if ~isKey(lmeMapElec, key); continue; end
        lme = lmeMapElec(key);
        coeffNames = lme.CoefficientNames; beta = lme.Coefficients.Estimate; covB = lme.CoefficientCovariance;
        contrastP = nan(6, 1); contrastData = cell(6, 9);
        rowIdx = 0;
        for p = 1:numel(POSITION_ORDER)
            posLvl = POSITION_ORDER{p};
            w = 0; if strcmp(posLvl, POSITION_SECOND); w = 1; end
            for c = 1:3
                lvl2 = CONTRASTS{c,1}; lvl1 = CONTRASTS{c,2}; label = CONTRASTS{c,3};
                rowIdx = rowIdx + 1;
                C = local_condition_position_contrast(coeffNames, REF_LEVEL, lvl1, lvl2, POSITION_SECOND, w, w);
                [est, se, tstat, p2, df2c, ciLow, ciHigh] = local_apply_contrast(lme, C, beta, covB);
                contrastP(rowIdx) = p2;
                contrastData(rowIdx, :) = {posLvl, label, est, se, df2c, tstat, p2, ciLow, ciHigh};
            end
        end
        [adjP, ~] = local_holm_bonferroni(contrastP, FDR_ALPHA);
        for rowIdx = 1:6
            pwIntElecRows(end+1, :) = [{elecName, metricName, bandName}, contrastData(rowIdx, 1:7), {adjP(rowIdx)}, contrastData(rowIdx, 8:9)]; %#ok<AGROW>
        end
    end
    PairwiseInteractionSimpleElectrodeT = local_build_table(pwIntElecRows, [1 2 3 4 5], [6 7 8 9 10 11 12 13], ...
        {'Electrode','Metric','Band','Position','Contrast','Estimate','SE','DF','tStat','pValue','pHolm','CI_Low','CI_High'});

    fprintf('Electrode pairwise/simple-effects (Task=%s): Condition-marginal=%d rows, Position-marginal=%d rows, Interaction-simple=%d rows.\n', ...
        TASK_LABEL, height(PairwiseConditionMarginalElectrodeT), height(PairwisePositionMarginalElectrodeT), height(PairwiseInteractionSimpleElectrodeT));

    %% ====================================================================
    %  SECTION 3j: RESIDUAL DIAGNOSTICS -- GLOBAL LEVEL ONLY, Primary +
    %  Secondary tier bands, both metrics
    % ====================================================================
    diagBands = [PRIMARY_TIER_BANDS, SECONDARY_TIER_BANDS];
    nDiagSaved = 0;
    for m = 1:nMetrics
        metricName = METRICS{m};
        for bi = 1:numel(diagBands)
            bandName = diagBands{bi};
            key = sprintf('%s|%s', bandName, metricName);
            if ~isKey(lmeMapGlobal, key); continue; end
            lme = lmeMapGlobal(key);
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
                    bandName, metricName, TASK_LABEL, ME.message); %#ok<AGROW>
            end
        end
    end
    fprintf('Saved %d residual diagnostic plots (Task=%s; Global level only, Primary+Secondary bands -- judgment call, see header).\n', nDiagSaved, TASK_LABEL);

    %% ====================================================================
    %  SECTION 3k: PARAMETERS
    % ====================================================================
    PARAMS18 = struct();
    PARAMS18.ScriptName = 'X18_position_int1_int2_intervention.m';
    PARAMS18.Task = TASK_LABEL;
    PARAMS18.GeneratedOn = datestr(now); %#ok<TNOW1,DATST>
    PARAMS18.Purpose = ['Exploratory analysis of a within-block ' ...
        'order/fatigue effect (Position: First vs Second), fit as Condition*Position with full interaction. ' ...
        'Distinct from DIAG01''s already-settled, purely arithmetic order-invariance finding -- see script header.'];
    PARAMS18.TaskSeparationRationale = 'Intervention recordings: single analysis; Position = InterventionNumber.';
    PARAMS18.ROIScopeDecision = 'ROI level intentionally omitted from this script (not required here); Global and Electrode levels only. See script header for how to extend to ROI level if needed.';
    PARAMS18.PositionDerivation = 'Intervention recordings: Position = InterventionNumber (1 -> First, 2 -> Second). Any other value excluded, never guessed.';
    PARAMS18.Model = 'LinearMixedEffects_RandomIntercept_FullInteraction';
    PARAMS18.FixedEffectsFormula = 'Metric ~ Condition*Position + (1|Subject)';
    PARAMS18.RandomEffects = '(1|Subject)';
    PARAMS18.FitMethod = 'REML';
    PARAMS18.DFMethod = 'Satterthwaite';
    PARAMS18.Metrics = METRICS;
    PARAMS18.ConditionOrder = CONDITION_ORDER;
    PARAMS18.PositionOrder = POSITION_ORDER;
    PARAMS18.BandTiers = struct('Primary', {PRIMARY_TIER_BANDS}, 'Secondary', {SECONDARY_TIER_BANDS}, 'Exploratory', {EXPLORATORY_TIER_BANDS}, 'ExploratoryBeta', {EXPLORATORY_BETA_TIER_BANDS});
    PARAMS18.FDRAlpha = FDR_ALPHA;
    PARAMS18.FDRFamilyScheme = ['THREE SEPARATE families per level (Global, Electrode), one per Term (Condition, ' ...
        'Position, Condition:Position), never pooled across Terms/levels/Tasks. Global: BH-FDR per (Tier x Metric), ' ...
        'across bands in tier (Primary=1, Secondary=3, Exploratory=5), matching X05''s tiered scheme. Electrode: ' ...
        'BH-FDR per (Band x Metric), across up to 62 electrodes, matching X09''s scheme.'];
    PARAMS18.PairwiseContrastMethod = ['Condition-main-effect survives -> 3 MARGINAL contrasts (E-N,E-R,N-R; equal-' ...
        'weighted over Position levels, including +/-0.5x the matching interaction coefficient), Holm(3). ' ...
        'Condition:Position interaction survives -> the same 3 contrasts computed SEPARATELY at each Position ' ...
        'level (6 simple-effects contrasts total), Holm-corrected together as one family of 6. Position-main-effect ' ...
        'survives -> 1 marginal First-vs-Second contrast (equal-weighted over Condition), no Holm correction (single comparison).'];
    PARAMS18.MarginalContrastCheck = ['Verified algebraically (see comment above local_condition_position_contrast): ' ...
        'under this model''s reference/dummy coding, the plain Condition_<level> coefficient equals the SIMPLE ' ...
        'Condition effect only at Position=First (the reference level); it does NOT equal the Position-marginal ' ...
        '(equal-weighted) Condition effect in general, which additionally requires +0.5x the relevant Condition: ' ...
        'Position interaction coefficient. The marginal contrast therefore explicitly includes that term rather ' ...
        'than assuming the reduction.'];
    PARAMS18.BaselineExclusionNote = sprintf('Global: %d Baseline row(s) excluded. Electrode: %d Baseline row(s) excluded. Baseline is never part of any inferential model in this pipeline.', nBaselineExcludedGlobal, nBaselineExcludedElec);
    PARAMS18.PositionExclusionCounts = struct('MissingAuditMatch', nMissingAudit, 'AmbiguousAuditMatch', nAmbiguous, 'UnparseableSEQ', nUnparseable, 'ExcludedGlobalRows', nExcludedGlobal, 'ExcludedElectrodeRows', nExcludedElec);
    PARAMS18.JudgmentCalls = { ...
        'Position reference level = First (arbitrary but explicit); Condition reference = E (locked, from X05).', ...
        'Marginal contrasts use EQUAL weighting across Position/Condition levels (standard EMM convention), not weighting by observed cell counts -- see header for rationale.', ...
        'Electrode reference list taken directly from X04''s BandPowerT.Channel (not from X07''s ROI/electrode-definitions file), consistent with the ROI-omission scope decision.', ...
        'Residual diagnostics scoped to Global level only, Primary+Secondary tier bands, both metrics (Electrode-level diagnostics not generated) -- mirrors X05/X11''s own scoping precedent.'};
    PARAMS18.SourcePhase0AuditFile = 'not used (Position from InterventionNumber)';
    PARAMS18.SourcePhase3BandPowerFile = bandPowerFile;
    PARAMS18.DesignVerification = struct('nRecordings', nRecordings, 'nSubjects', nSubjects, ...
        'PositionCounts', posCounts, 'CondPositionConfound', condPosConfound, 'IntegrityIssues', {integrityIssues});
    PARAMS18.RunDateTime = datestr(now); %#ok<TNOW1,DATST>
    PARAMS18.RunMode = RUN_MODE;
    if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
        PARAMS18.ExcludedSubjects = EXCLUDE_NOISY_SUBJECTS;
        PARAMS18.ExcludedSubjectsRationale = EXCLUSION_RATIONALE_intervention;
        PARAMS18.RowsExcludedForSensitivity = struct( ...
            'GlobalLevel', nExcludedRowsSensitivityGlobal, ...
            'ElectrodeLevel', nExcludedRowsSensitivityElec);
    else
        PARAMS18.ExcludedSubjects = [];
    end

    save(fullfile(STATS_DIR, 'phase18_parameters.mat'), 'PARAMS18');

    %% ====================================================================
    %  SECTION 3l: SAVE OUTPUTS
    % ====================================================================
    save(fullfile(STATS_DIR, 'phase18_omnibus_global.mat'), 'OmnibusGlobalT');
    save(fullfile(STATS_DIR, 'phase18_omnibus_electrode.mat'), 'OmnibusElectrodeT');
    save(fullfile(STATS_DIR, 'phase18_pairwise_condition_marginal_global.mat'), 'PairwiseConditionMarginalGlobalT');
    save(fullfile(STATS_DIR, 'phase18_pairwise_condition_marginal_electrode.mat'), 'PairwiseConditionMarginalElectrodeT');
    save(fullfile(STATS_DIR, 'phase18_pairwise_position_marginal_global.mat'), 'PairwisePositionMarginalGlobalT');
    save(fullfile(STATS_DIR, 'phase18_pairwise_position_marginal_electrode.mat'), 'PairwisePositionMarginalElectrodeT');
    save(fullfile(STATS_DIR, 'phase18_pairwise_interaction_simple_global.mat'), 'PairwiseInteractionSimpleGlobalT');
    save(fullfile(STATS_DIR, 'phase18_pairwise_interaction_simple_electrode.mat'), 'PairwiseInteractionSimpleElectrodeT');

    writetable(OmnibusGlobalT, fullfile(STATS_DIR, 'phase18_omnibus_global.csv'));
    writetable(OmnibusElectrodeT, fullfile(STATS_DIR, 'phase18_omnibus_electrode.csv'));
    if ~isempty(PairwiseConditionMarginalGlobalT); writetable(PairwiseConditionMarginalGlobalT, fullfile(STATS_DIR, 'phase18_pairwise_condition_marginal_global.csv')); end
    if ~isempty(PairwiseConditionMarginalElectrodeT); writetable(PairwiseConditionMarginalElectrodeT, fullfile(STATS_DIR, 'phase18_pairwise_condition_marginal_electrode.csv')); end
    if ~isempty(PairwisePositionMarginalGlobalT); writetable(PairwisePositionMarginalGlobalT, fullfile(STATS_DIR, 'phase18_pairwise_position_marginal_global.csv')); end
    if ~isempty(PairwisePositionMarginalElectrodeT); writetable(PairwisePositionMarginalElectrodeT, fullfile(STATS_DIR, 'phase18_pairwise_position_marginal_electrode.csv')); end
    if ~isempty(PairwiseInteractionSimpleGlobalT); writetable(PairwiseInteractionSimpleGlobalT, fullfile(STATS_DIR, 'phase18_pairwise_interaction_simple_global.csv')); end
    if ~isempty(PairwiseInteractionSimpleElectrodeT); writetable(PairwiseInteractionSimpleElectrodeT, fullfile(STATS_DIR, 'phase18_pairwise_interaction_simple_electrode.csv')); end

    fprintf('\nSaved output to: %s\n', STATS_DIR);

    %% ====================================================================
    %  SECTION 3m: REPORT
    % ====================================================================
    reportLines = {};
    reportLines{end+1} = repmat('=', 1, 60);
    reportLines{end+1} = sprintf('PHASE 18 (X18) -- POSITION (FIRST vs SECOND) EFFECT REPORT (Task=%s)', TASK_LABEL);
    reportLines{end+1} = repmat('=', 1, 60);
    reportLines{end+1} = '';
    reportLines{end+1} = 'Exploratory. Position = intervention slot (Int1 = First, Int2 = Second). Tests whether the';
    reportLines{end+1} = 'Condition effect depends on whether the condition was the 1st or 2nd intervention (order/carry-over).';
    reportLines{end+1} = '';
    reportLines{end+1} = sprintf('Baseline excluded: %d Global rows, %d Electrode rows.', nBaselineExcludedGlobal, nBaselineExcludedElec);
    reportLines{end+1} = sprintf('Position derivation (from InterventionNumber): %d recordings checked; excluded (InterventionNumber not 1/2) = %d.', ...
        nRecCombos, nUnparseable);
    reportLines{end+1} = sprintf('Retained after Position derivation: %d Global rows, %d Electrode rows.', height(GlobalBandPowerT), height(BandPowerT));
    reportLines{end+1} = sprintf('Recordings: %d | Subjects: %d | Position counts: First=%d, Second=%d.', nRecordings, nSubjects, posCounts.First, posCounts.Second);
    reportLines{end+1} = sprintf('Condition x Position confound detected: %d', condPosConfound);
    reportLines{end+1} = '';
    reportLines{end+1} = 'Model: Metric ~ Condition*Position + (1|Subject), REML, Satterthwaite DF (full interaction).';
    reportLines{end+1} = 'Metrics: LogAbsolutePower and RelativePower, co-primary, never pooled for FDR.';
    reportLines{end+1} = 'Levels analyzed: Global (whole-scalp) and Electrode (up to 62 channels). ROI level intentionally omitted (see header).';
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Global-level omnibus results, by Term (3 independent FDR families) ---';
    for ti = 1:numel(TERMS)
        termName = TERMS{ti};
        idxT = strcmp(OmnibusGlobalT.Term, termName);
        reportLines{end+1} = sprintf('  Term=%s: %d significant / %d tested', termName, sum(idxT & OmnibusGlobalT.SignificantFDR), sum(idxT)); %#ok<AGROW>
        subRows = find(idxT);
        for r = subRows'
            reportLines{end+1} = sprintf('    [%s] %-25s Tier=%-11s F(%g,%.1f)=%.3f p=%.4f pFDR=%.4f family_n=%d %s', ...
                OmnibusGlobalT.Metric{r}, OmnibusGlobalT.Band{r}, OmnibusGlobalT.Tier{r}, OmnibusGlobalT.DF1(r), OmnibusGlobalT.DF2(r), ...
                OmnibusGlobalT.FStat(r), OmnibusGlobalT.pValue(r), OmnibusGlobalT.pFDR(r), OmnibusGlobalT.FamilySize(r), ...
                local_yn(OmnibusGlobalT.SignificantFDR(r))); %#ok<AGROW>
        end
    end
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Electrode-level omnibus results, by Term (3 independent FDR families; summary only, see .csv for all electrodes) ---';
    for ti = 1:numel(TERMS)
        termName = TERMS{ti};
        idxT = strcmp(OmnibusElectrodeT.Term, termName);
        reportLines{end+1} = sprintf('  Term=%s: %d significant / %d tested', termName, sum(idxT & OmnibusElectrodeT.SignificantFDR), sum(idxT)); %#ok<AGROW>
    end
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Pairwise / simple-effects contrasts (Global level) ---';
    reportLines{end+1} = sprintf('  Condition-marginal (from Condition FDR survival, Holm-3): %d rows', height(PairwiseConditionMarginalGlobalT));
    for r = 1:height(PairwiseConditionMarginalGlobalT)
        reportLines{end+1} = sprintf('    [%s] %-25s %-8s est=%.4f (95%% CI %.4f to %.4f) t=%.3f df=%.1f p=%.4f pHolm=%.4f', ...
            PairwiseConditionMarginalGlobalT.Metric{r}, PairwiseConditionMarginalGlobalT.Band{r}, PairwiseConditionMarginalGlobalT.Contrast{r}, ...
            PairwiseConditionMarginalGlobalT.Estimate(r), PairwiseConditionMarginalGlobalT.CI_Low(r), PairwiseConditionMarginalGlobalT.CI_High(r), ...
            PairwiseConditionMarginalGlobalT.tStat(r), PairwiseConditionMarginalGlobalT.DF(r), PairwiseConditionMarginalGlobalT.pValue(r), PairwiseConditionMarginalGlobalT.pHolm(r)); %#ok<AGROW>
    end
    reportLines{end+1} = sprintf('  Position-marginal (from Position FDR survival, no correction, single comparison): %d rows', height(PairwisePositionMarginalGlobalT));
    for r = 1:height(PairwisePositionMarginalGlobalT)
        reportLines{end+1} = sprintf('    [%s] %-25s %-16s est=%.4f (95%% CI %.4f to %.4f) t=%.3f df=%.1f p=%.4f', ...
            PairwisePositionMarginalGlobalT.Metric{r}, PairwisePositionMarginalGlobalT.Band{r}, PairwisePositionMarginalGlobalT.Contrast{r}, ...
            PairwisePositionMarginalGlobalT.Estimate(r), PairwisePositionMarginalGlobalT.CI_Low(r), PairwisePositionMarginalGlobalT.CI_High(r), ...
            PairwisePositionMarginalGlobalT.tStat(r), PairwisePositionMarginalGlobalT.DF(r), PairwisePositionMarginalGlobalT.pValue(r)); %#ok<AGROW>
    end
    reportLines{end+1} = sprintf('  Interaction-simple-effects (from Interaction FDR survival, Holm-6 per unit): %d rows', height(PairwiseInteractionSimpleGlobalT));
    for r = 1:height(PairwiseInteractionSimpleGlobalT)
        reportLines{end+1} = sprintf('    [%s] %-25s Pos=%-6s %-8s est=%.4f (95%% CI %.4f to %.4f) t=%.3f df=%.1f p=%.4f pHolm=%.4f', ...
            PairwiseInteractionSimpleGlobalT.Metric{r}, PairwiseInteractionSimpleGlobalT.Band{r}, PairwiseInteractionSimpleGlobalT.Position{r}, PairwiseInteractionSimpleGlobalT.Contrast{r}, ...
            PairwiseInteractionSimpleGlobalT.Estimate(r), PairwiseInteractionSimpleGlobalT.CI_Low(r), PairwiseInteractionSimpleGlobalT.CI_High(r), ...
            PairwiseInteractionSimpleGlobalT.tStat(r), PairwiseInteractionSimpleGlobalT.DF(r), PairwiseInteractionSimpleGlobalT.pValue(r), PairwiseInteractionSimpleGlobalT.pHolm(r)); %#ok<AGROW>
    end
    reportLines{end+1} = '';
    reportLines{end+1} = sprintf('--- Pairwise / simple-effects contrasts (Electrode level; counts only, see .csv for full detail) ---');
    reportLines{end+1} = sprintf('  Condition-marginal rows: %d | Position-marginal rows: %d | Interaction-simple rows: %d', ...
        height(PairwiseConditionMarginalElectrodeT), height(PairwisePositionMarginalElectrodeT), height(PairwiseInteractionSimpleElectrodeT));
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Judgment calls (disclosed) ---';
    for i = 1:numel(PARAMS18.JudgmentCalls)
        reportLines{end+1} = ['  - ' PARAMS18.JudgmentCalls{i}]; %#ok<AGROW>
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
    fid = fopen(fullfile(REPORT_DIR, sprintf('X18_position_report_%s_%s.txt', TASK_LABEL, RUN_MODE_LABEL)), 'w');
    fprintf(fid, '%s\n', reportText);
    fclose(fid);

end % TASK_TYPES loop

fprintf('\n=== X18 complete for both Tasks. Purely additive exploratory script -- X01-X17 outputs are untouched. ===\n');

%% ======================= LOCAL FUNCTIONS (all after top-level code) =======================
%  local_subject_numeric, local_fdr_bh, local_holm_bonferroni,
%  local_condition_contrast, local_build_table are copied VERBATIM from
%  X05_primary_statistics.m (self-contained-file convention,
%  same as every other script in this pipeline). local_robust_cellstr,
%  local_position_from_seq, local_find_interaction_name,
%  local_condition_position_contrast, local_position_marginal_contrast,
%  local_fit_condition_position_omnibus, local_apply_contrast, local_yn,
%  iif_str are NEW helpers written for this script's Condition*Position
%  interaction model, which X05 (Condition-only) has no equivalent for.

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
            error('Phase18:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name1);
        end
        cRow(idx1) = cRow(idx1) - 1;
    end
    if ~strcmp(level2, refLevel)
        name2 = sprintf('%s_%s', varName, level2);
        idx2 = find(strcmp(coeffNames, name2));
        if isempty(idx2)
            error('Phase18:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name2);
        end
        cRow(idx2) = cRow(idx2) + 1;
    end
end

function T = local_build_table(rows, textColIdx, numColIdx, varNames)
    nCols = numel(varNames);
    if size(rows, 2) ~= nCols
        error('Phase18:BuildTableColumnMismatch', ...
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
            error('Phase18:BuildTableColumnMismatch', 'Column %d not classified as text or numeric.', j);
        end
    end
    T = table(cols{:}, 'VariableNames', varNames);
end

function s = local_yn(tf)
    if tf
        s = 'SIGNIFICANT';
    else
        s = 'n.s.';
    end
end

function s = iif_str(cond, a, b)
    if cond
        s = a;
    else
        s = b;
    end
end

function c = local_robust_cellstr(col)
% Normalizes a readtable column of unknown inferred type (cell/string/
% char/logical/numeric) to cellstr of char, so downstream strcmp/strcmpi
% comparisons are reliable regardless of how readtable happened to infer
% that column's type from the .csv round-trip (defensive; the pipeline's
% own writetable/normalize_text_columns convention writes text as
% cellstr, but a logical column such as ParseOK may come back as
% logical, numeric 0/1, or text depending on the CSV serialization).
    if iscell(col)
        c = col;
    elseif isstring(col)
        c = cellstr(col);
    elseif ischar(col)
        c = cellstr(col);
    elseif islogical(col)
        c = cellstr(string(col));
    elseif isnumeric(col)
        c = cellstr(num2str(col(:)));
    else
        c = cellstr(col);
    end
end

function posLabel = local_position_from_seq(taskKey, seqU)
    posLabel = '';
    switch taskKey
        case 'nback'
            if strcmp(seqU, 'NS')
                posLabel = 'First';
            elseif strcmp(seqU, 'SN')
                posLabel = 'Second';
            end
        case 'sart'
            if strcmp(seqU, 'SN')
                posLabel = 'First';
            elseif strcmp(seqU, 'NS')
                posLabel = 'Second';
            end
    end
end

function nm = local_find_interaction_name(coeffNames, condLevel, posLevel)
    cand1 = sprintf('Condition_%s:Position_%s', condLevel, posLevel);
    cand2 = sprintf('Position_%s:Condition_%s', posLevel, condLevel);
    if any(strcmp(coeffNames, cand1))
        nm = cand1;
    elseif any(strcmp(coeffNames, cand2))
        nm = cand2;
    else
        error('Phase18:InteractionNameNotFound', ...
            'Neither "%s" nor "%s" found in CoefficientNames.', cand1, cand2);
    end
end

function cRow = local_condition_position_contrast(coeffNames, refLevel, level1, level2, posSecondLabel, w1, w2)
% Builds the contrast vector for (Condition=level2 - Condition=level1),
% at a Position "mixture" set by w1/w2 -- the weight applied to each
% level's OWN Condition:Position interaction coefficient (interaction
% coefficients only exist for non-reference Condition levels; the
% reference level E contributes exactly 0 by construction under
% reference/dummy coding).
%
% ALGEBRAIC DERIVATION (documented here, not re-derived elsewhere): for
% model Y = b0 + bLevel*I(Condition=level) + bPos*I(Position=Second) +
% bLevel_Pos*I(Condition=level)*I(Position=Second) + ..., the simple
% Condition effect of level vs. refLevel AT Position=p (p=0 for First,
% p=1 for Second) is (bLevel2-bLevel1) + p*(bLevel2_Pos - bLevel1_Pos).
% At p=0 this is exactly (bLevel2-bLevel1), i.e. local_condition_contrast
% with no interaction term -- CONFIRMING the "reduces to the plain
% Condition-coefficient contrast" case holds ONLY at Position's
% reference level. The equal-weighted MARGINAL effect (average of p=0
% and p=1) is (bLevel2-bLevel1) + 0.5*(bLevel2_Pos - bLevel1_Pos), which
% does NOT drop the interaction term -- CONFIRMING it does NOT reduce to
% the plain coefficient in general (only would if the interaction
% coefficients were exactly zero). This is exactly what w1=w2=0 vs.
% w1=w2=0.5 implement below.
    cRow = local_condition_contrast(coeffNames, 'Condition', refLevel, level1, level2);
    if ~strcmp(level1, refLevel) && w1 ~= 0
        name1 = local_find_interaction_name(coeffNames, level1, posSecondLabel);
        idx1 = find(strcmp(coeffNames, name1));
        cRow(idx1) = cRow(idx1) - w1;
    end
    if ~strcmp(level2, refLevel) && w2 ~= 0
        name2 = local_find_interaction_name(coeffNames, level2, posSecondLabel);
        idx2 = find(strcmp(coeffNames, name2));
        cRow(idx2) = cRow(idx2) + w2;
    end
end

function cRow = local_position_marginal_contrast(coeffNames, conditionOrder, posSecondLabel)
% Builds the contrast vector for the Position MARGINAL effect
% (Second - First), equal-weighted over ALL Condition levels (including
% the reference level E, whose interaction contribution is exactly 0
% since Condition=E has no Condition:Position interaction coefficient
% under reference/dummy coding).
    cRow = zeros(1, numel(coeffNames));
    idxPos = find(strcmp(coeffNames, sprintf('Position_%s', posSecondLabel)));
    if isempty(idxPos)
        error('Phase18:ContrastNameNotFound', 'Coefficient name "Position_%s" not found in model.', posSecondLabel);
    end
    cRow(idxPos) = cRow(idxPos) + 1;
    nCond = numel(conditionOrder);
    refLevel = conditionOrder{1};
    for i = 1:nCond
        lvl = conditionOrder{i};
        if strcmp(lvl, refLevel); continue; end
        nm = local_find_interaction_name(coeffNames, lvl, posSecondLabel);
        idx = find(strcmp(coeffNames, nm));
        cRow(idx) = cRow(idx) + 1/nCond;
    end
end

function fr = local_fit_condition_position_omnibus(subT, metricName, conditionOrder, positionOrder)
% Fits Metric ~ Condition*Position + (1|Subject), REML, and extracts the
% THREE omnibus Satterthwaite terms (Condition, Position, and the
% Condition:Position interaction). Returns a struct rather than erroring
% out, so the calling loops can log a failure (via fr.errMsg) and
% continue, exactly like X05/X09/X11's own try/catch pattern around
% fitlme.
    fr = struct('ok', false, 'lme', [], 'rows', struct('Term', {}, 'FStat', {}, 'DF1', {}, 'DF2', {}, 'pValue', {}), ...
        'warnMsg', '', 'errMsg', '');
    if isempty(subT) || ~ismember(metricName, subT.Properties.VariableNames)
        fr.errMsg = 'Empty subT or missing metric column.';
        return;
    end
    subT.Subject = categorical(local_subject_numeric(subT.Subject));
    subT.Condition = categorical(cellstr(subT.Condition), conditionOrder);
    subT.Position = categorical(cellstr(subT.Position), positionOrder);
    formula = sprintf('%s ~ Condition*Position + (1|Subject)', metricName);
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
        fr.errMsg = ['anova() failed: ' ME.message];
        return;
    end
    idxCond = find(strcmp(aov.Term, 'Condition'), 1);
    idxPos  = find(strcmp(aov.Term, 'Position'), 1);
    idxInt  = find(strcmp(aov.Term, 'Condition:Position'), 1);
    if isempty(idxInt)
        idxInt = find(strcmp(aov.Term, 'Position:Condition'), 1);
    end
    if isempty(idxCond) || isempty(idxPos) || isempty(idxInt)
        fr.errMsg = sprintf('anova(lme) missing one or more expected terms (Condition/Position/interaction). Terms found: %s', strjoin(aov.Term, ', '));
        return;
    end
    fr.ok = true;
    fr.lme = lme;
    fr.rows = struct('Term', {'Condition', 'Position', 'Condition:Position'}, ...
        'FStat', {aov.FStat(idxCond), aov.FStat(idxPos), aov.FStat(idxInt)}, ...
        'DF1', {aov.DF1(idxCond), aov.DF1(idxPos), aov.DF1(idxInt)}, ...
        'DF2', {aov.DF2(idxCond), aov.DF2(idxPos), aov.DF2(idxInt)}, ...
        'pValue', {aov.pValue(idxCond), aov.pValue(idxPos), aov.pValue(idxInt)});
end

function [est, se, tstat, p, df2c, ciLow, ciHigh] = local_apply_contrast(lme, C, beta, covB)
% Applies a 1 x nCoeff contrast row vector C to a fitted lme, returning
% estimate/SE/t/p/df/CI -- identical arithmetic to X05/X09/X11's own
% inline pairwise-contrast block, factored out since this script builds
% many more distinct contrast vectors (marginal / simple-effects /
% position) than those scripts' single Condition-only contrast.
    est = C * beta;
    se = sqrt(C * covB * C');
    [p, ~, ~, df2c] = coefTest(lme, C, 0, 'DFMethod', 'satterthwaite');
    tstat = est / se;
    tcrit = tinv(0.975, df2c);
    ciLow = est - tcrit * se;
    ciHigh = est + tcrit * se;
end
