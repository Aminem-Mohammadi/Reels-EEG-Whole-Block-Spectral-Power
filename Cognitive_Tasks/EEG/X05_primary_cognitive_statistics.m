%% ========================================================================
% X05_primary_cognitive_statistics.m
%
% Purpose: primary E vs N vs R statistical analysis of global (whole-
% scalp) band power. Runs as three separate, parallel analyses, one per
% Task (Full / N-back / SART); Task is never modeled as a factor and
% Task x Condition is never fit -- each Task's data is filtered out and
% analyzed with its own independent copy of the same model. Condition is
% filtered to {E, N, R} before modeling in every Task (Baseline rows are
% excluded; the excluded count is logged per Task, never silently
% dropped). Design verification (full-factorial/incomplete-design check)
% re-runs independently inside each Task's own E/N/R-filtered data.
%
% Method: LMM (Condition fixed 3-level + (1|Subject), REML, Satterthwaite
% DF) fit on 14 outcome units (Primary = Alpha; Secondary = Delta, Theta,
% Beta; Exploratory = 3 gamma bands) x 2 co-primary metrics
% (LogAbsolutePower, RelativePower, never pooled for FDR). BH-FDR within
% each Tier x Metric family (the Primary tier has family size 1 and is
% reported as-is, never treated as pre-corrected). Pairwise contrasts
% (E-N, E-R, N-R) are Holm-corrected for units that survive FDR. Four
% pre-specified sensitivity analyses: (1) all-bands-pooled FDR; (2) a
% single-subject exclusion check, not run here -- no equivalent
% subject-specific exclusion has been identified for this dataset, and
% this slot is logged as skipped in sensitivityLog rather than silently
% dropped or reused from a different dataset; (3) BlockNumber (session
% order, 2 or 3, since Block 1 is Baseline and is excluded here) added as
% an extra fixed-effect covariate, with a check for whether Condition and
% order are fully confounded; (4) a paired-subset consistency check;
% plus residual diagnostics.
%
% Analysis level: global scalp band power (X04's GlobalBandPowerT) -- one
% row per Subject x Task x BlockNumber x Condition x Band. Channel-level
% and ROI/topographic analysis are separate scripts (X07 onward).
%
% Requirements: base MATLAB + Statistics and Machine Learning Toolbox
% (fitlme, anova on LinearMixedModel, coefTest, ttest, signrank, qqplot).
%
% Inputs (read-only, not modified):
%   <ANALYSIS_ROOT>/04_subject_level/phase3_subject_level_bandpower.mat
%     (GlobalBandPowerT, PARAMS3)
%
% Outputs (one complete set per Task):
%   <ANALYSIS_ROOT>/<Full|Nback|SART>/05_primary_stats/
%       phase4_omnibus_results.mat / .xlsx
%       phase4_pairwise_contrasts.mat / .xlsx
%       phase4_sensitivity_results.mat / .xlsx
%       phase4_paired_subset_checks.xlsx
%       phase4_parameters.mat
%       diagnostic_residuals_<band>_<metric>.png
%   <ANALYSIS_ROOT>/09_reports/
%       phase4_primary_statistics_report_<Full|Nback|SART>.txt
%
% RUN_MODE: full_sample / sensitivity_excl_noisy (see
% EXCLUSION_CONFIG_cognitive.m); run this script once per mode.
% ========================================================================

clear; clc;

%% Section 0: Paths and existence checks
config;
DATA_ROOT_FOR_SIBLING = DATA_DIR;
[DATA_ROOT_PARENT, ~] = fileparts(DATA_ROOT_FOR_SIBLING);
ANALYSIS_ROOT = OUT_DIR;

SUBJ_DIR   = fullfile(ANALYSIS_ROOT, '04_subject_level');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');

if ~exist(SUBJ_DIR, 'dir')
    error('Phase3:MissingInput', 'X04 output directory not found: %s', SUBJ_DIR);
end
bandPowerFile = fullfile(SUBJ_DIR, 'phase3_subject_level_bandpower.mat');
if ~exist(bandPowerFile, 'file')
    error('Phase3:MissingInput', 'Required X04 file not found: %s', bandPowerFile);
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
    error('Phase4:BadRunMode', 'RUN_MODE must be one of: %s', strjoin(RUN_MODES_cognitive, ', '));
end
RUN_MODE_LABEL = RUN_MODES_DISPLAY_cognitive.(RUN_MODE);
fprintf('RUN_MODE = %s\n', RUN_MODE_LABEL);
if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
    fprintf('Excluding subject(s) {%s} from every model in this run -- see EXCLUSION_CONFIG_cognitive.m for rationale.\n', ...
        strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','));
end
% ---- END NEW ----

if isempty(ver('stats'))
    error('Phase4:MissingToolbox', ...
        'Statistics and Machine Learning Toolbox not detected. This script requires fitlme, anova(LinearMixedModel), coefTest, ttest, signrank, and qqplot.');
end

%% Section 1: Load X04 output (shared across all three Task analyses)
S = load(bandPowerFile, 'GlobalBandPowerT', 'PARAMS3');
GlobalBandPowerT_ALL = S.GlobalBandPowerT;
PARAMS3 = S.PARAMS3;
clear S;

requiredCols = {'Subject','TaskType','BlockNumber','Condition','Band', ...
                'AbsolutePower','LogAbsolutePower','RelativePower'};
missingCols = requiredCols(~ismember(requiredCols, GlobalBandPowerT_ALL.Properties.VariableNames));
if ~isempty(missingCols)
    error('Phase4:MissingColumns', 'GlobalBandPowerT is missing expected column(s): %s', strjoin(missingCols, ', '));
end

ALL_BANDS = fieldnames(PARAMS3.BandDefinitionsHz);
METRICS   = {'LogAbsolutePower', 'RelativePower'};
CONDITION_ORDER = {'E', 'N', 'R'};   % Baseline excluded from this analysis (see header item 2)
REF_LEVEL = CONDITION_ORDER{1};

PRIMARY_TIER_BANDS     = {'Alpha'};
SECONDARY_TIER_BANDS   = {'Delta', 'Theta', 'Beta'};
EXPLORATORY_TIER_BANDS = {'Beta_low_exploratory', 'Beta_high_exploratory', ...
    'Gamma_low_exploratory', 'Gamma_high1_exploratory', 'Gamma_high2_exploratory'};
% Beta_low/Beta_high added (a sub-division of the unchanged locked
% Beta[13,30] band); each gets its own separate BH-FDR
% family exactly like the pre-existing Gamma exploratory bands.

expectedTierBands = [PRIMARY_TIER_BANDS, SECONDARY_TIER_BANDS, EXPLORATORY_TIER_BANDS];
if ~isempty(setxor(sort(expectedTierBands), sort(ALL_BANDS')))
    error('Phase4:BandMismatch', ...
        ['The hard-coded band tiers do not exactly match PARAMS3.BandDefinitionsHz. ' ...
         'Tiers: %s | Phase3 bands: %s'], strjoin(expectedTierBands, ','), strjoin(ALL_BANDS', ','));
end

TASK_TYPES = {'full', 'nback', 'sart'};
TASK_DISPLAY = struct('full', 'Full', 'nback', 'Nback', 'sart', 'SART');

fprintf('========================================================\n');
fprintf('PHASE 4 -- PRIMARY E vs N vs R STATISTICS, COGNITIVE-TASK DATA\n');
fprintf('Running as THREE INDEPENDENT, PARALLEL analyses (Full / Nback / SART).\n');
fprintf('========================================================\n\n');

%% ========================================================================
%  TASK LOOP -- everything below is repeated independently per Task
% ========================================================================
for taskI = 1:numel(TASK_TYPES)
    TASK_KEY = TASK_TYPES{taskI};
    TASK_LABEL = TASK_DISPLAY.(TASK_KEY);

    STATS_DIR = fullfile(ANALYSIS_ROOT, TASK_LABEL, '05_primary_stats', RUN_MODE);
    if ~exist(STATS_DIR, 'dir'); mkdir(STATS_DIR); end

    fprintf('\n########################################################\n');
    fprintf('### TASK = %s | RUN_MODE = %s\n', TASK_LABEL, RUN_MODE_LABEL);
    fprintf('########################################################\n\n');

    %% Section 1b: filter to this Task, then to Condition in {E,N,R}
    thisTaskT = GlobalBandPowerT_ALL(strcmp(GlobalBandPowerT_ALL.TaskType, TASK_KEY), :);
    nBaselineExcluded = sum(strcmp(thisTaskT.Condition, 'Baseline'));
    GlobalBandPowerT = thisTaskT(ismember(thisTaskT.Condition, CONDITION_ORDER), :);
    fprintf('Task=%s: %d total rows loaded, %d Baseline rows excluded, %d E/N/R rows retained for modeling.\n', ...
        TASK_LABEL, height(thisTaskT), nBaselineExcluded, height(GlobalBandPowerT));

    % ---- NEW: apply the shared noisy-subject exclusion list (see
    % EXCLUSION_CONFIG_cognitive.m), only when RUN_MODE=sensitivity_excl_noisy.
    % Never applied to the full_sample run. ----
    nExcludedRowsSensitivity = 0;
    if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
        subjNumForExcl = local_subject_numeric(GlobalBandPowerT.Subject);
        excludeMaskSensitivity = ismember(subjNumForExcl, EXCLUDE_NOISY_SUBJECTS);
        nExcludedRowsSensitivity = sum(excludeMaskSensitivity);
        GlobalBandPowerT = GlobalBandPowerT(~excludeMaskSensitivity, :);
        fprintf('Task=%s: RUN_MODE=%s -- excluded %d row(s) belonging to subject(s) {%s}; %d rows remain for modeling.\n', ...
            TASK_LABEL, RUN_MODE_LABEL, nExcludedRowsSensitivity, ...
            strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','), height(GlobalBandPowerT));
    end

    if isempty(GlobalBandPowerT)
        warning('Task=%s: no E/N/R rows found after filtering -- skipping this Task entirely.', TASK_LABEL);
        continue;
    end

    integrityIssues = {};
    sensitivityLog  = {};

    %% Section 2: Data-driven design verification (never assumed)
    comboT = unique(GlobalBandPowerT(:, {'Subject','BlockNumber','Condition'}), 'rows');
    nRecordings = height(comboT);
    fprintf('Section 2: %d unique Subject x BlockNumber x Condition combinations found.\n', nRecordings);

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
    if numel(pairKeys) > 3
        integrityIssues{end+1} = sprintf('Unexpected pair-type key encountered during design verification for Task=%s.', TASK_LABEL);
    end
    allSubjectsComplete3 = false;
    condsPerSubjCounts = cellfun(@numel, subjPairs);
    if all(condsPerSubjCounts >= 3)
        allSubjectsComplete3 = true;
        integrityIssues{end+1} = sprintf(['Task=%s: every subject appears to have 3+ recordings/conditions -- the incomplete-design ' ...
            'justification for using an LMM instead of RM-ANOVA no longer applies. STOP and re-examine.'], TASK_LABEL);
    end

    % Condition x BlockNumber confound check (population level) -- the
    % structural analog of the intervention pipeline's Condition x
    % InterventionNumber check (see header item 4).
    condList = cellstr(comboT.Condition);
    blockList = comboT.BlockNumber;
    confoundFlag = false;
    for c = 1:numel(CONDITION_ORDER)
        thisCond = CONDITION_ORDER{c};
        blocksForCond = blockList(strcmp(condList, thisCond));
        uBlocks = unique(blocksForCond);
        blockCounts = arrayfun(@(u) sum(blocksForCond == u), uBlocks);
        fprintf('  Condition %s: BlockNumber counts = %s\n', thisCond, mat2str(blockCounts'));
        if numel(uBlocks) < 2
            confoundFlag = true;
            integrityIssues{end+1} = sprintf(['Task=%s: Condition %s occurs at only one BlockNumber value -- ' ...
                'Condition and block order are FULLY confounded for this condition. The order-covariate ' ...
                'sensitivity analysis (Section 9) cannot separate them and will be skipped/flagged.'], TASK_LABEL, thisCond); %#ok<AGROW>
        end
    end

    if ~isempty(integrityIssues)
        fprintf(2, '\n*** DESIGN VERIFICATION FLAGGED %d ISSUE(S) (Task=%s) -- see report. ***\n\n', numel(integrityIssues), TASK_LABEL);
        for i = 1:numel(integrityIssues)
            fprintf(2, '  - %s\n', integrityIssues{i});
        end
    end

    %% Section 3: PARAMS4 (all agreed methodological decisions, explicit)
    PARAMS4 = struct();
    PARAMS4.Task = TASK_LABEL;
    PARAMS4.AnalysisLevel = sprintf('GlobalScalpBandPower (X04 GlobalBandPowerT), Task=%s only, Condition in {E,N,R} only (Baseline excluded)', TASK_LABEL);
    PARAMS4.Model = 'LinearMixedEffects_RandomIntercept';
    PARAMS4.FixedEffectsPrimary = 'Condition (3-level: E, N, R; reference = E)';
    PARAMS4.RandomEffects = '(1|Subject)';
    PARAMS4.FitMethod = 'REML';
    PARAMS4.DFMethod = 'Satterthwaite';
    PARAMS4.Metrics = METRICS;
    PARAMS4.BandTiers = struct('Primary', {PRIMARY_TIER_BANDS}, 'Secondary', {SECONDARY_TIER_BANDS}, ...
                                'Exploratory', {EXPLORATORY_TIER_BANDS});
    PARAMS4.TierHierarchyPreSpecified = true;
    PARAMS4.FDRAlpha = 0.05;
    PARAMS4.FDRFamilyScheme = 'Separate BH-FDR family per (Tier x Metric); metrics NEVER pooled; Primary tier (n=1) has no correction applied.';
    PARAMS4.PairwiseContrastMethod = 'Holm-Bonferroni across the 3 pairwise contrasts (N-E, R-E, R-N), computed only for bands whose tiered-family omnibus test survives FDR.';
    PARAMS4.BaselineExclusionNote = sprintf('%d Baseline row(s) excluded before any modeling; Baseline is audited/aggregated upstream but never enters this or any other inferential model in this pipeline.', nBaselineExcluded);
    PARAMS4.Sensitivity2Note = 'Sensitivity 2 (exclude a specific known-outlier subject) is NOT run for this dataset -- see script header item 5.';
    PARAMS4.SensitivityAnalyses = {'AllBandsPooledFDR_perMetric', 'BlockNumberCovariateAdjusted', 'PairedSubsetConsistencyCheck'};
    PARAMS4.DesignVerification = struct('nRecordings', nRecordings, 'nSubjects', nSubjects, ...
        'PairTypeCounts', {pairTypeCount.keys}, 'ConfoundFlag', confoundFlag, 'IntegrityIssues', {integrityIssues});
    PARAMS4.SourcePhase3File = bandPowerFile;
    PARAMS4.RunDateTime = datestr(now); %#ok<TNOW1,DATST>
    PARAMS4.RunMode = RUN_MODE;
    if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
        PARAMS4.ExcludedSubjects = EXCLUDE_NOISY_SUBJECTS;
        PARAMS4.ExcludedSubjectsRationale = EXCLUSION_RATIONALE_cognitive;
        PARAMS4.RowsExcludedForSensitivity = nExcludedRowsSensitivity;
    else
        PARAMS4.ExcludedSubjects = [];
    end

    save(fullfile(STATS_DIR, 'phase4_parameters.mat'), 'PARAMS4');

    if allSubjectsComplete3
        error('Phase4:DesignAssumptionViolated', ...
            'Task=%s: design verification found a complete 3-condition-per-subject structure. Stopping before fitting any model -- see integrityIssues.', TASK_LABEL);
    end

    %% Section 4: Primary tiered analysis -- 14 core LMM fits (7 bands x 2 metrics)
    nBands = numel(ALL_BANDS);
    nMetrics = numel(METRICS);
    omnibusRows = cell(0, 6);
    lmeStore = cell(nBands, nMetrics);
    convergenceWarnings = {};

    for m = 1:nMetrics
        metricName = METRICS{m};
        for b = 1:nBands
            bandName = ALL_BANDS{b};
            subT = GlobalBandPowerT(strcmp(GlobalBandPowerT.Band, bandName), :);
            subT.Subject = categorical(local_subject_numeric(subT.Subject));
            subT.Condition = categorical(cellstr(subT.Condition), CONDITION_ORDER);
            formula = sprintf('%s ~ Condition + (1|Subject)', metricName);
            try
                lastwarn('');
                lme = fitlme(subT, formula, 'FitMethod', 'REML');
                [warnMsg, ~] = lastwarn();
                if ~isempty(warnMsg)
                    convergenceWarnings{end+1} = sprintf('%s / %s: %s', metricName, bandName, warnMsg); %#ok<AGROW>
                end
            catch ME
                integrityIssues{end+1} = sprintf('fitlme FAILED for %s / %s (Task=%s): %s', metricName, bandName, TASK_LABEL, ME.message); %#ok<AGROW>
                lmeStore{b, m} = [];
                continue;
            end
            lmeStore{b, m} = lme;
            aov = anova(lme, 'DFMethod', 'satterthwaite');
            rowIdx = find(strcmp(aov.Term, 'Condition'), 1);
            if isempty(rowIdx)
                integrityIssues{end+1} = sprintf('anova(lme) has no Condition row for %s / %s (Task=%s).', metricName, bandName, TASK_LABEL); %#ok<AGROW>
                continue;
            end
            omnibusRows(end+1, :) = {metricName, bandName, aov.FStat(rowIdx), aov.DF1(rowIdx), aov.DF2(rowIdx), aov.pValue(rowIdx)}; %#ok<AGROW>
        end
    end

    OmnibusT = local_build_table(omnibusRows, [1 2], [3 4 5 6], {'Metric','Band','FStat','DF1','DF2','pValue'});

    %% Section 5: FDR within pre-specified tier x metric families
    OmnibusT.Tier = repmat({'?'}, height(OmnibusT), 1);
    OmnibusT.pFDR = nan(height(OmnibusT), 1);
    OmnibusT.SignificantFDR = false(height(OmnibusT), 1);
    OmnibusT.FamilySize = nan(height(OmnibusT), 1);

    tierDefs = {'Primary', PRIMARY_TIER_BANDS; 'Secondary', SECONDARY_TIER_BANDS; 'Exploratory', EXPLORATORY_TIER_BANDS};
    for t = 1:size(tierDefs, 1)
        tierName = tierDefs{t, 1};
        tierBands = tierDefs{t, 2};
        for m = 1:nMetrics
            metricName = METRICS{m};
            idx = strcmp(OmnibusT.Metric, metricName) & ismember(OmnibusT.Band, tierBands);
            OmnibusT.Tier(idx) = {tierName};
            OmnibusT.FamilySize(idx) = sum(idx);
            pvals = OmnibusT.pValue(idx);
            if numel(pvals) == 1
                OmnibusT.pFDR(idx) = pvals;
                OmnibusT.SignificantFDR(idx) = pvals <= PARAMS4.FDRAlpha;
            else
                [adjP, reject] = local_fdr_bh(pvals, PARAMS4.FDRAlpha);
                OmnibusT.pFDR(idx) = adjP;
                OmnibusT.SignificantFDR(idx) = reject;
            end
        end
    end

    fprintf('\n--- Primary tiered omnibus results (Task=%s) ---\n', TASK_LABEL);
    disp(OmnibusT);

    %% Section 6: Pairwise contrasts for bands significant under the primary tiered scheme
    pairwiseRows = cell(0, 11);
    CONTRASTS = {'N','E','N - E'; 'R','E','R - E'; 'R','N','R - N'};

    sigMask = OmnibusT.SignificantFDR;
    for r = find(sigMask)'
        metricName = OmnibusT.Metric{r};
        bandName = OmnibusT.Band{r};
        b = find(strcmp(ALL_BANDS, bandName));
        m = find(strcmp(METRICS, metricName));
        lme = lmeStore{b, m};
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
            [p, ~, df1c, df2c] = coefTest(lme, C, 0, 'DFMethod', 'satterthwaite');
            tstat = est / se;
            tcrit = tinv(0.975, df2c);
            ciLow = est - tcrit * se;
            ciHigh = est + tcrit * se;
            contrastP(c) = p;
            contrastData(c, :) = {label, est, se, df2c, tstat, p, ciLow, ciHigh};
        end
        [adjP, ~] = local_holm_bonferroni(contrastP, PARAMS4.FDRAlpha);
        for c = 1:3
            pairwiseRows(end+1, :) = [{metricName, bandName}, contrastData(c, 1:6), {adjP(c)}, contrastData(c, 7:8)]; %#ok<AGROW>
        end
    end
    PairwiseT = local_build_table(pairwiseRows, [1 2 3], [4 5 6 7 8 9 10 11], ...
        {'Metric','Band','Contrast','Estimate','SE','DF','tStat','pValue','pHolm','CI_Low','CI_High'});

    fprintf('\n--- Pairwise contrasts (Task=%s; bands significant under primary tiered scheme) ---\n', TASK_LABEL);
    disp(PairwiseT);

    %% Section 7: Sensitivity 1 -- all-7-bands pooled FDR family per metric
    SensAllBandsT = OmnibusT(:, {'Metric','Band','pValue'});
    SensAllBandsT.pFDR_pooled = nan(height(SensAllBandsT), 1);
    SensAllBandsT.SignificantFDR_pooled = false(height(SensAllBandsT), 1);
    for m = 1:nMetrics
        metricName = METRICS{m};
        idx = strcmp(SensAllBandsT.Metric, metricName);
        [adjP, reject] = local_fdr_bh(SensAllBandsT.pValue(idx), PARAMS4.FDRAlpha);
        SensAllBandsT.pFDR_pooled(idx) = adjP;
        SensAllBandsT.SignificantFDR_pooled(idx) = reject;
    end
    SensAllBandsT.SignificantFDR_tiered = OmnibusT.SignificantFDR;
    SensAllBandsT.ConclusionMatches = SensAllBandsT.SignificantFDR_pooled == SensAllBandsT.SignificantFDR_tiered;
    nFlipped = sum(~SensAllBandsT.ConclusionMatches);
    sensitivityLog{end+1} = sprintf('Sensitivity 1 (all-7-bands pooled FDR): %d/%d band-metric conclusions changed vs. the primary tiered scheme.', nFlipped, height(SensAllBandsT));

    %% Section 8: Sensitivity 2 -- NOT RUN for this dataset (see header item 5)
    sensitivityLog{end+1} = 'Sensitivity 2 (exclude a specific known-outlier subject, e.g. "sub07" in the intervention pipeline): SKIPPED for the cognitive-task dataset -- sub07 was specific to the intervention dataset''s own prior diagnostics and has no established meaning here. See script header item 5.';
    SensNoSub7T = table();

    %% Section 9: Sensitivity 3 -- order (BlockNumber) covariate (was
    %  InterventionNumber in the intervention pipeline; see header item 4)
    if confoundFlag
        sensitivityLog{end+1} = 'Sensitivity 3 (BlockNumber covariate): SKIPPED -- Condition and BlockNumber are fully confounded for at least one condition (see Section 2).';
        SensOrderT = table();
    else
        sensOrderRows = cell(0, 6);
        for m = 1:nMetrics
            metricName = METRICS{m};
            for b = 1:nBands
                bandName = ALL_BANDS{b};
                subT = GlobalBandPowerT(strcmp(GlobalBandPowerT.Band, bandName), :);
                subT.Subject = categorical(local_subject_numeric(subT.Subject));
                subT.Condition = categorical(cellstr(subT.Condition), CONDITION_ORDER);
                subT.BlockNumber = categorical(subT.BlockNumber);
                formula = sprintf('%s ~ Condition + BlockNumber + (1|Subject)', metricName);
                try
                    lme3 = fitlme(subT, formula, 'FitMethod', 'REML');
                    aov3 = anova(lme3, 'DFMethod', 'satterthwaite');
                    rowIdx = find(strcmp(aov3.Term, 'Condition'), 1);
                    sensOrderRows(end+1, :) = {metricName, bandName, aov3.FStat(rowIdx), aov3.DF1(rowIdx), aov3.DF2(rowIdx), aov3.pValue(rowIdx)}; %#ok<AGROW>
                catch ME
                    sensitivityLog{end+1} = sprintf('Sensitivity 3 fit failed for %s / %s (Task=%s): %s', metricName, bandName, TASK_LABEL, ME.message); %#ok<AGROW>
                end
            end
        end
        SensOrderT = local_build_table(sensOrderRows, [1 2], [3 4 5 6], {'Metric','Band','FStat','DF1','DF2','pValue'});
        SensOrderT.SignificantFDR = false(height(SensOrderT), 1);
        for t = 1:size(tierDefs, 1)
            tierBands = tierDefs{t, 2};
            for m = 1:nMetrics
                metricName = METRICS{m};
                idx = strcmp(SensOrderT.Metric, metricName) & ismember(SensOrderT.Band, tierBands);
                pvals = SensOrderT.pValue(idx);
                if numel(pvals) == 1
                    SensOrderT.SignificantFDR(idx) = pvals <= PARAMS4.FDRAlpha;
                else
                    [~, reject] = local_fdr_bh(pvals, PARAMS4.FDRAlpha);
                    SensOrderT.SignificantFDR(idx) = reject;
                end
            end
        end
        matchMask = strcat(SensOrderT.Metric, '|', SensOrderT.Band);
        primMask  = strcat(OmnibusT.Metric, '|', OmnibusT.Band);
        [~, locB] = ismember(matchMask, primMask);
        SensOrderT.SignificantFDR_primary = OmnibusT.SignificantFDR(locB);
        SensOrderT.ConclusionMatches = SensOrderT.SignificantFDR == SensOrderT.SignificantFDR_primary;
        nFlipped3 = sum(~SensOrderT.ConclusionMatches);
        sensitivityLog{end+1} = sprintf('Sensitivity 3 (BlockNumber covariate added): %d/%d band-metric conclusions changed vs. the primary (order-free) analysis.', nFlipped3, height(SensOrderT));
    end

    %% Section 10: Sensitivity 4 -- paired-subset consistency check
    pairedCheckRows = cell(0, 9);
    for r = 1:height(PairwiseT)
        metricName = PairwiseT.Metric{r};
        bandName = PairwiseT.Band{r};
        contrastLabel = PairwiseT.Contrast{r};
        parts = strsplit(contrastLabel, ' - ');
        lvl2 = parts{1}; lvl1 = parts{2};

        subjInPair = uniqSubj(cellfun(@(x) all(ismember(sort({lvl1,lvl2}), sort(x))), subjPairs));
        bandT = GlobalBandPowerT(strcmp(GlobalBandPowerT.Band, bandName), :);
        bandT.SubjNum = local_subject_numeric(bandT.Subject);

        x1 = nan(numel(subjInPair), 1);
        x2 = nan(numel(subjInPair), 1);
        for i = 1:numel(subjInPair)
            v1 = bandT.(metricName)(bandT.SubjNum == subjInPair(i) & strcmp(cellstr(bandT.Condition), lvl1));
            v2 = bandT.(metricName)(bandT.SubjNum == subjInPair(i) & strcmp(cellstr(bandT.Condition), lvl2));
            if numel(v1) == 1 && numel(v2) == 1
                x1(i) = v1; x2(i) = v2;
            end
        end
        valid = ~isnan(x1) & ~isnan(x2);
        x1 = x1(valid); x2 = x2(valid);
        if numel(x1) < 3
            pairedCheckRows(end+1, :) = {metricName, bandName, contrastLabel, numel(x1), PairwiseT.Estimate(r), NaN, NaN, NaN, 'Too few paired subjects'}; %#ok<AGROW>
            continue;
        end
        diffVec = x2 - x1;
        [~, pT] = ttest(diffVec);
        pW = signrank(diffVec);
        pairedCheckRows(end+1, :) = {metricName, bandName, contrastLabel, numel(x1), PairwiseT.Estimate(r), mean(diffVec), pT, pW, ''}; %#ok<AGROW>
    end
    PairedCheckT = local_build_table(pairedCheckRows, [1 2 3 9], [4 5 6 7 8], ...
        {'Metric','Band','Contrast','N_PairedSubjects','LMM_Estimate','Paired_MeanDiff','PairedT_p','Wilcoxon_p','Note'});

    %% Section 11: Residual diagnostics (normality/homoscedasticity of the LMM)
    diagBands = unique([PRIMARY_TIER_BANDS, SECONDARY_TIER_BANDS]);
    for m = 1:nMetrics
        metricName = METRICS{m};
        for bi = 1:numel(diagBands)
            bandName = diagBands{bi};
            b = find(strcmp(ALL_BANDS, bandName));
            lme = lmeStore{b, m};
            if isempty(lme); continue; end
            try
                res = residuals(lme, 'ResidualType', 'Pearson');
                fig = figure('Visible', 'off', 'Position', [100 100 900 400]);
                subplot(1,2,1);
                histogram(res, 15);
                title(sprintf('Residuals: %s (%s), Task=%s', bandName, metricName, TASK_LABEL), 'Interpreter', 'none');
                xlabel('Pearson residual'); ylabel('Count');
                subplot(1,2,2);
                qqplot(res);
                title('Normal Q-Q');
                outPng = fullfile(STATS_DIR, sprintf('diagnostic_residuals_%s_%s.png', bandName, metricName));
                saveas(fig, outPng);
                close(fig);
            catch ME
                sensitivityLog{end+1} = sprintf('Residual diagnostic failed for %s / %s (Task=%s): %s', metricName, bandName, TASK_LABEL, ME.message); %#ok<AGROW>
            end
        end
    end

    %% Section 12: Save outputs
    save(fullfile(STATS_DIR, 'phase4_omnibus_results.mat'), 'OmnibusT', 'SensAllBandsT');
    save(fullfile(STATS_DIR, 'phase4_pairwise_contrasts.mat'), 'PairwiseT');
    save(fullfile(STATS_DIR, 'phase4_sensitivity_results.mat'), 'SensAllBandsT', 'SensNoSub7T', 'SensOrderT', 'sensitivityLog');
    writetable(OmnibusT, fullfile(STATS_DIR, 'phase4_omnibus_results.xlsx'));
    if ~isempty(PairwiseT)
        writetable(PairwiseT, fullfile(STATS_DIR, 'phase4_pairwise_contrasts.xlsx'));
    end
    writetable(SensAllBandsT, fullfile(STATS_DIR, 'phase4_sensitivity_allbands.xlsx'));
    if ~isempty(SensOrderT);  writetable(SensOrderT,  fullfile(STATS_DIR, 'phase4_sensitivity_order.xlsx'));  end
    writetable(PairedCheckT, fullfile(STATS_DIR, 'phase4_paired_subset_checks.xlsx'));

    %% Section 13: Report
    reportLines = {};
    reportLines{end+1} = repmat('=', 1, 60);
    reportLines{end+1} = sprintf('PHASE 4 -- PRIMARY E vs N vs R STATISTICS REPORT (Task=%s)', TASK_LABEL);
    reportLines{end+1} = repmat('=', 1, 60);
    reportLines{end+1} = '';
    reportLines{end+1} = sprintf('Baseline rows excluded before modeling: %d', nBaselineExcluded);
    reportLines{end+1} = sprintf('Recordings verified: %d | Subjects: %d', nRecordings, nSubjects);
    reportLines{end+1} = sprintf('Condition-pair coverage: %s', strjoin(cellfun(@(k) sprintf('%s=%d', strrep(k,'_',' & '), pairTypeCount(k)), pairKeys, 'UniformOutput', false), ', '));
    reportLines{end+1} = sprintf('Condition x BlockNumber confound detected: %d', confoundFlag);
    reportLines{end+1} = '';
    reportLines{end+1} = 'Model: Linear mixed-effects, Condition (fixed, 3-level) + (1|Subject), REML, Satterthwaite DF.';
    reportLines{end+1} = sprintf('Analysis level: Global scalp band power (X04 Structure C), Task=%s only.', TASK_LABEL);
    reportLines{end+1} = 'Tiers: Primary=Alpha; Secondary=Delta/Theta/Beta; Exploratory=3 gamma bands (pre-specified, not data-driven).';
    reportLines{end+1} = 'Metrics: LogAbsolutePower and RelativePower, co-primary, parallel, never pooled for FDR.';
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Primary tiered omnibus results ---';
    for r = 1:height(OmnibusT)
        reportLines{end+1} = sprintf('  [%s] %-25s Tier=%-11s F(%g,%.1f)=%.3f p=%.4f pFDR=%.4f family_n=%d %s', ...
            OmnibusT.Metric{r}, OmnibusT.Band{r}, OmnibusT.Tier{r}, OmnibusT.DF1(r), OmnibusT.DF2(r), ...
            OmnibusT.FStat(r), OmnibusT.pValue(r), OmnibusT.pFDR(r), OmnibusT.FamilySize(r), ...
            local_yn(OmnibusT.SignificantFDR(r))); %#ok<AGROW>
    end
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Pairwise contrasts (bands significant under primary tiered scheme; Holm-corrected within band) ---';
    if isempty(PairwiseT)
        reportLines{end+1} = '  (none -- no band/metric survived its tier''s FDR correction)';
    else
        for r = 1:height(PairwiseT)
            reportLines{end+1} = sprintf('  [%s] %-25s %-8s est=%.4f (95%% CI %.4f to %.4f) t=%.3f df=%.1f p=%.4f pHolm=%.4f', ...
                PairwiseT.Metric{r}, PairwiseT.Band{r}, PairwiseT.Contrast{r}, PairwiseT.Estimate(r), ...
                PairwiseT.CI_Low(r), PairwiseT.CI_High(r), PairwiseT.tStat(r), PairwiseT.DF(r), ...
                PairwiseT.pValue(r), PairwiseT.pHolm(r)); %#ok<AGROW>
        end
    end
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Sensitivity analyses ---';
    for i = 1:numel(sensitivityLog)
        reportLines{end+1} = ['  - ' sensitivityLog{i}]; %#ok<AGROW>
    end
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Paired-subset consistency checks (descriptive; not independently multiplicity-corrected) ---';
    for r = 1:height(PairedCheckT)
        reportLines{end+1} = sprintf('  [%s] %-25s %-8s N=%d LMM_est=%.4f paired_meandiff=%s pairedT_p=%s Wilcoxon_p=%s %s', ...
            PairedCheckT.Metric{r}, PairedCheckT.Band{r}, PairedCheckT.Contrast{r}, PairedCheckT.N_PairedSubjects(r), ...
            PairedCheckT.LMM_Estimate(r), local_fmt(PairedCheckT.Paired_MeanDiff(r)), local_fmt(PairedCheckT.PairedT_p(r)), ...
            local_fmt(PairedCheckT.Wilcoxon_p(r)), PairedCheckT.Note{r}); %#ok<AGROW>
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
    fid = fopen(fullfile(REPORT_DIR, sprintf('phase4_primary_statistics_report_%s_%s.txt', TASK_LABEL, RUN_MODE_LABEL)), 'w');
    fprintf(fid, '%s\n', reportText);
    fclose(fid);

end % TASK_TYPES loop

fprintf('\n=== X05 complete for all three Tasks. Do not proceed to X06 until every Task''s output has been reviewed. ===\n');

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
            error('Phase4:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name1);
        end
        cRow(idx1) = cRow(idx1) - 1;
    end
    if ~strcmp(level2, refLevel)
        name2 = sprintf('%s_%s', varName, level2);
        idx2 = find(strcmp(coeffNames, name2));
        if isempty(idx2)
            error('Phase4:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name2);
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
        error('Phase4:BuildTableColumnMismatch', ...
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
            error('Phase4:BuildTableColumnMismatch', 'Column %d not classified as text or numeric.', j);
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
