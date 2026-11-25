%% ========================================================================
% X05_primary_statistics.m
%
% Whole-scalp (global) band-power statistics.
% Linear mixed model Metric ~ Condition + (1|Subject), REML, Satterthwaite df, for
% ln(absolute) and relative power. Omnibus tests are BH-FDR corrected within families
% (alpha alone; delta/theta/beta; beta subdivisions; gamma bands); pairwise contrasts
% E - N, E - R, N - R are Holm-corrected across the three contrasts. Also writes all
% contrasts without gating, an order-covariate model, and within-person paired checks for
% participants who received both conditions (paired t and Wilcoxon signed-rank tests).
% Inputs : 04_subject_level.
% Outputs: 05_primary_stats/<RUN_MODE>/ (phase4_*.xlsx), 09_reports/.
% Run twice, with RUN_MODE = 'full_sample' (N = 30) and RUN_MODE = 'sensitivity_excl_noisy'
% (cleaned sample, N = 28; excluded subjects listed in EXCLUSION_CONFIG_intervention.m).
% Outputs are written to <OUT_DIR>/<folder>/<RUN_MODE>/ and reports to <OUT_DIR>/09_reports/.
% ========================================================================

clear; clc;

%% Section 0: Paths and existence checks
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
STATS_DIR  = fullfile(ANALYSIS_ROOT, '05_primary_stats', RUN_MODE);
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');

if ~exist(SUBJ_DIR, 'dir')
    error('Phase3:MissingInput', 'Phase 3 output directory not found: %s', SUBJ_DIR);
end
bandPowerFile = fullfile(SUBJ_DIR, 'phase3_subject_level_bandpower.mat');
if ~exist(bandPowerFile, 'file')
    error('Phase3:MissingInput', 'Required Phase 3 file not found: %s', bandPowerFile);
end
if ~exist(STATS_DIR, 'dir');  mkdir(STATS_DIR);  end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

if isempty(ver('stats'))
    error('Phase4:MissingToolbox', ...
        'Statistics and Machine Learning Toolbox not detected. This script requires fitlme, anova(LinearMixedModel), coefTest, ttest, signrank, and qqplot.');
end

%% Section 1: Load Phase 3 output
S = load(bandPowerFile, 'GlobalBandPowerT', 'PARAMS3');
GlobalBandPowerT = S.GlobalBandPowerT;
[GlobalBandPowerT, RUN_MODE_N_EXCL_ROWS, RUN_MODE_EXCL_SUBJ] = apply_run_mode_filter(GlobalBandPowerT, RUN_MODE, EXCLUDE_NOISY_SUBJECTS);
if RUN_MODE_N_EXCL_ROWS > 0; fprintf('RUN_MODE=%s: removed %d row(s) of subject(s) {%s} from GlobalBandPowerT.\n', RUN_MODE_LABEL, RUN_MODE_N_EXCL_ROWS, num2str(RUN_MODE_EXCL_SUBJ)); end

PARAMS3 = S.PARAMS3;
clear S;

requiredCols = {'Subject','InterventionNumber','Condition','Band', ...
                'AbsolutePower','LogAbsolutePower','RelativePower'};
missingCols = requiredCols(~ismember(requiredCols, GlobalBandPowerT.Properties.VariableNames));
if ~isempty(missingCols)
    error('Phase4:MissingColumns', 'GlobalBandPowerT is missing expected column(s): %s', strjoin(missingCols, ', '));
end

ALL_BANDS = fieldnames(PARAMS3.BandDefinitionsHz);
METRICS   = {'LogAbsolutePower', 'RelativePower'};
CONDITION_ORDER = {'E', 'N', 'R'};
REF_LEVEL = CONDITION_ORDER{1}; % 'E' is the reference level for dummy coding

PRIMARY_TIER_BANDS     = {'Alpha'};
SECONDARY_TIER_BANDS   = {'Delta', 'Theta', 'Beta'};
EXPLORATORY_TIER_BANDS = {'Gamma_low_exploratory', 'Gamma_high1_exploratory', 'Gamma_high2_exploratory'};
EXPLORATORY_BETA_TIER_BANDS = {'Beta_low_exploratory', 'Beta_high_exploratory'};

expectedTierBands = [PRIMARY_TIER_BANDS, SECONDARY_TIER_BANDS, EXPLORATORY_TIER_BANDS, EXPLORATORY_BETA_TIER_BANDS];
if ~isempty(setxor(sort(expectedTierBands), sort(ALL_BANDS')))
    error('Phase4:BandMismatch', ...
        ['The hard-coded band tiers do not exactly match PARAMS3.BandDefinitionsHz. ' ...
         'Tiers: %s | Phase3 bands: %s'], strjoin(expectedTierBands, ','), strjoin(ALL_BANDS', ','));
end

integrityIssues = {};
sensitivityLog  = {};

%% Section 2: Data-driven design check (never assumed)
comboT = unique(GlobalBandPowerT(:, {'Subject','InterventionNumber','Condition'}), 'rows');
nRecordings = height(comboT);
fprintf('Section 2: %d unique Subject x InterventionNumber x Condition combinations found.\n', nRecordings);

subjNumeric = local_subject_numeric(comboT.Subject);
uniqSubj = unique(subjNumeric);
nSubjects = numel(uniqSubj);

subjPairs = cell(nSubjects, 1);
pairTypeCount = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:nSubjects
    rowsI = subjNumeric == uniqSubj(i);
    condsI = cellstr(comboT.Condition(rowsI));
    if numel(condsI) ~= 2
        integrityIssues{end+1} = sprintf('Subject %g has %d recordings (expected 2).', uniqSubj(i), numel(condsI)); %#ok<AGROW>
    elseif strcmp(condsI{1}, condsI{2})
        integrityIssues{end+1} = sprintf('Subject %g has the SAME condition (%s) in both recordings -- not a valid pair.', uniqSubj(i), condsI{1}); %#ok<AGROW>
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

fprintf('Condition-pair coverage across %d subjects:\n', nSubjects);
pairKeys = keys(pairTypeCount);
for k = 1:numel(pairKeys)
    fprintf('  %s : %d subjects\n', strrep(pairKeys{k}, '_', ' & '), pairTypeCount(pairKeys{k}));
end
if numel(pairKeys) > 3
    integrityIssues{end+1} = 'Unexpected pair-type key encountered during design verification.';
end
allSubjectsComplete3 = false;
condsPerSubjCounts = cellfun(@numel, subjPairs);
if all(condsPerSubjCounts >= 3)
    allSubjectsComplete3 = true;
    integrityIssues{end+1} = ['Every subject appears to have 3+ recordings/conditions -- the incomplete-design ' ...
        'justification for using an LMM instead of RM-ANOVA no longer applies to this dataset. STOP and re-examine.'];
end

% Condition x InterventionNumber confound check (population level)
condList = cellstr(comboT.Condition);
intList  = comboT.InterventionNumber;
confoundFlag = false;
for c = 1:numel(CONDITION_ORDER)
    thisCond = CONDITION_ORDER{c};
    intsForCond = intList(strcmp(condList, thisCond));
    uInts = unique(intsForCond);
    intCounts = arrayfun(@(u) sum(intsForCond == u), uInts);
    fprintf('  Condition %s: InterventionNumber counts = %s\n', thisCond, mat2str(intCounts'));
    if numel(uInts) < 2
        confoundFlag = true;
        integrityIssues{end+1} = sprintf(['Condition %s occurs at only one InterventionNumber value -- ' ...
            'Condition and session order are FULLY confounded for this condition. The order-covariate ' ...
            'sensitivity analysis (Section 9) cannot separate them and will be skipped/flagged.'], thisCond); %#ok<AGROW>
    end
end

if ~isempty(integrityIssues)
    fprintf(2, '\n*** DESIGN VERIFICATION FLAGGED %d ISSUE(S) -- see report. ***\n\n', numel(integrityIssues));
    for i = 1:numel(integrityIssues)
        fprintf(2, '  - %s\n', integrityIssues{i});
    end
end

%% Section 3: PARAMS4 (methodological settings)
PARAMS4 = struct();
PARAMS4.AnalysisLevel = 'GlobalScalpBandPower (Phase 3 Structure C, GlobalBandPowerT)';
PARAMS4.Model = 'LinearMixedEffects_RandomIntercept';
PARAMS4.FixedEffectsPrimary = 'Condition (3-level: E, N, R; reference = E)';
PARAMS4.RandomEffects = '(1|Subject)';
PARAMS4.FitMethod = 'REML';
PARAMS4.DFMethod = 'Satterthwaite';
PARAMS4.Metrics = METRICS; % co-primary, parallel, never pooled for FDR
PARAMS4.BandTiers = struct('Primary', {PRIMARY_TIER_BANDS}, 'Secondary', {SECONDARY_TIER_BANDS}, ...
                            'Exploratory', {EXPLORATORY_TIER_BANDS}, 'ExploratoryBeta', {EXPLORATORY_BETA_TIER_BANDS});
PARAMS4.TierHierarchyPreSpecified = true;
PARAMS4.FDRAlpha = 0.05;
PARAMS4.FDRFamilyScheme = 'Separate BH-FDR family per (Tier x Metric); metrics NEVER pooled; Primary tier (n=1) has no correction applied.';
PARAMS4.PairwiseContrastMethod = 'Holm-Bonferroni across the 3 pairwise contrasts (E-N, E-R, N-R), computed only for bands whose tiered-family omnibus test survives FDR.';
PARAMS4.Sub07Status = 'Included in the primary analysis per explicit instruction. Excluded only in a dedicated sensitivity re-fit (Section 8).';
PARAMS4.SensitivityAnalyses = {'AllBandsPooledFDR_perMetric', 'ExcludeSub07', 'OrderCovariateAdjusted', 'PairedSubsetConsistencyCheck'};
PARAMS4.DesignVerification = struct('nRecordings', nRecordings, 'nSubjects', nSubjects, ...
    'PairTypeCounts', {pairTypeCount.keys}, 'ConfoundFlag', confoundFlag, 'IntegrityIssues', {integrityIssues});
PARAMS4.SourcePhase3File = bandPowerFile;
PARAMS4.RunDateTime = datestr(now); %#ok<TNOW1,DATST>

save(fullfile(STATS_DIR, 'phase4_parameters.mat'), 'PARAMS4');

if allSubjectsComplete3
    error('Phase4:DesignAssumptionViolated', ...
        'Design verification found a complete 3-condition-per-subject structure. Stopping before fitting any model -- see integrityIssues.');
end

%% Section 4: Primary tiered analysis -- 14 core LMM fits (7 bands x 2 metrics)
nBands = numel(ALL_BANDS);
nMetrics = numel(METRICS);
omnibusRows = cell(0, 6);   % Metric, Band, F, DF1, DF2, pValue
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
            integrityIssues{end+1} = sprintf('fitlme FAILED for %s / %s: %s', metricName, bandName, ME.message); %#ok<AGROW>
            lmeStore{b, m} = [];
            continue;
        end
        lmeStore{b, m} = lme;
        aov = anova(lme, 'DFMethod', 'satterthwaite');
        rowIdx = find(strcmp(aov.Term, 'Condition'), 1);
        if isempty(rowIdx)
            integrityIssues{end+1} = sprintf('anova(lme) has no Condition row for %s / %s.', metricName, bandName); %#ok<AGROW>
            continue;
        end
        omnibusRows(end+1, :) = {metricName, bandName, aov.FStat(rowIdx), aov.DF1(rowIdx), aov.DF2(rowIdx), aov.pValue(rowIdx)}; %#ok<AGROW>
    end
end

OmnibusT = local_build_table(omnibusRows, [1 2], [3 4 5 6], {'Metric','Band','FStat','DF1','DF2','pValue'});

%% Section 5: FDR within pre-specified tier x metric families (primary analysis)
OmnibusT.Tier = repmat({'?'}, height(OmnibusT), 1);
OmnibusT.pFDR = nan(height(OmnibusT), 1);
OmnibusT.SignificantFDR = false(height(OmnibusT), 1);
OmnibusT.FamilySize = nan(height(OmnibusT), 1);

tierDefs = {'Primary', PRIMARY_TIER_BANDS; 'Secondary', SECONDARY_TIER_BANDS; 'Exploratory', EXPLORATORY_TIER_BANDS; 'ExploratoryBeta', EXPLORATORY_BETA_TIER_BANDS};
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

fprintf('\n--- Primary tiered omnibus results ---\n');
disp(OmnibusT);

%% Section 6: Pairwise contrasts for bands significant under the primary tiered scheme
pairwiseRows = cell(0, 11); % Metric,Band,Contrast,Estimate,SE,DF,t,p,pHolm,CI_low,CI_high
CONTRASTS = {'E','N','E - N'; 'E','R','E - R'; 'N','R','N - R'};   % {level2, level1, label}; estimate = level2 - level1

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

fprintf('\n--- Pairwise contrasts (bands significant under primary tiered scheme) ---\n');
disp(PairwiseT);

allRows = cell(0, 11);
for b = 1:nBands
    for m = 1:nMetrics
        lme = lmeStore{b, m};
        if isempty(lme); continue; end
        cp = nan(3, 1); cdat = cell(3, 8);
        for c = 1:3
            C = local_condition_contrast(lme.CoefficientNames, 'Condition', REF_LEVEL, CONTRASTS{c, 2}, CONTRASTS{c, 1});
            est = C * lme.Coefficients.Estimate; se = sqrt(C * lme.CoefficientCovariance * C');
            [p, ~, ~, df2c] = coefTest(lme, C, 0, 'DFMethod', 'satterthwaite');
            tc = tinv(0.975, df2c);
            cp(c) = p; cdat(c, :) = {CONTRASTS{c, 3}, est, se, df2c, est / se, p, est - tc*se, est + tc*se};
        end
        adj = local_holm_bonferroni(cp, PARAMS4.FDRAlpha);
        for c = 1:3
            allRows(end+1, :) = [{METRICS{m}, ALL_BANDS{b}}, cdat(c, 1:6), {adj(c)}, cdat(c, 7:8)]; %#ok<AGROW>
        end
    end
end
AllContrastsT = local_build_table(allRows, [1 2 3], [4 5 6 7 8 9 10 11], ...
    {'Metric','Band','Contrast','Estimate','SE','DF','tStat','pValue','pHolm_withinModel','CI_Low','CI_High'});

%% Section 7: Sensitivity 1 -- all-bands pooled FDR family per metric
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
sensitivityLog{end+1} = sprintf('Sensitivity 1 (all-bands pooled FDR): %d/%d band-metric conclusions changed vs. the primary tiered scheme.', nFlipped, height(SensAllBandsT));

%% Section 8: Sensitivity 2 -- exclude sub07
sub07Numeric = local_subject_numeric(GlobalBandPowerT.Subject);
noSub7T = GlobalBandPowerT(sub07Numeric ~= 7, :);
nExcludedRows = height(GlobalBandPowerT) - height(noSub7T);
if nExcludedRows == 0
    sensitivityLog{end+1} = 'Sensitivity 2 (exclude sub07): WARNING -- filter matched 0 rows. Subject numbering convention may differ from assumed; sensitivity analysis NOT performed.';
    SensNoSub7T = table();
else
    sensNoSub7Rows = cell(0, 6);
    for m = 1:nMetrics
        metricName = METRICS{m};
        for b = 1:nBands
            bandName = ALL_BANDS{b};
            subT = noSub7T(strcmp(noSub7T.Band, bandName), :);
            subT.Subject = categorical(local_subject_numeric(subT.Subject));
            subT.Condition = categorical(cellstr(subT.Condition), CONDITION_ORDER);
            formula = sprintf('%s ~ Condition + (1|Subject)', metricName);
            try
                lme2 = fitlme(subT, formula, 'FitMethod', 'REML');
                aov2 = anova(lme2, 'DFMethod', 'satterthwaite');
                rowIdx = find(strcmp(aov2.Term, 'Condition'), 1);
                sensNoSub7Rows(end+1, :) = {metricName, bandName, aov2.FStat(rowIdx), aov2.DF1(rowIdx), aov2.DF2(rowIdx), aov2.pValue(rowIdx)}; %#ok<AGROW>
            catch ME
                sensitivityLog{end+1} = sprintf('Sensitivity 2 fit failed for %s / %s: %s', metricName, bandName, ME.message); %#ok<AGROW>
            end
        end
    end
    SensNoSub7T = local_build_table(sensNoSub7Rows, [1 2], [3 4 5 6], {'Metric','Band','FStat','DF1','DF2','pValue'});
    SensNoSub7T.SignificantFDR = false(height(SensNoSub7T), 1);
    for t = 1:size(tierDefs, 1)
        tierBands = tierDefs{t, 2};
        for m = 1:nMetrics
            metricName = METRICS{m};
            idx = strcmp(SensNoSub7T.Metric, metricName) & ismember(SensNoSub7T.Band, tierBands);
            pvals = SensNoSub7T.pValue(idx);
            if numel(pvals) == 1
                SensNoSub7T.SignificantFDR(idx) = pvals <= PARAMS4.FDRAlpha;
            else
                [~, reject] = local_fdr_bh(pvals, PARAMS4.FDRAlpha);
                SensNoSub7T.SignificantFDR(idx) = reject;
            end
        end
    end
    matchMask = strcat(SensNoSub7T.Metric, '|', SensNoSub7T.Band);
    primMask  = strcat(OmnibusT.Metric, '|', OmnibusT.Band);
    [~, locB] = ismember(matchMask, primMask);
    if any(locB == 0); error('Phase4:SensitivityRowMismatch', 'Sensitivity 2 produced a band/metric row with no primary counterpart.'); end
    SensNoSub7T.SignificantFDR_primary = OmnibusT.SignificantFDR(locB);
    SensNoSub7T.ConclusionMatches = SensNoSub7T.SignificantFDR == SensNoSub7T.SignificantFDR_primary;
    nFlipped2 = sum(~SensNoSub7T.ConclusionMatches);
    sensitivityLog{end+1} = sprintf('Sensitivity 2 (exclude sub07, N=%d rows removed): %d/%d band-metric conclusions changed vs. the primary (sub07-included) analysis.', nExcludedRows, nFlipped2, height(SensNoSub7T));
end

%% Section 9: Sensitivity 3 -- order (InterventionNumber) covariate
if confoundFlag
    sensitivityLog{end+1} = 'Sensitivity 3 (order covariate): SKIPPED -- Condition and InterventionNumber are fully confounded for at least one condition (see Section 2).';
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
            subT.InterventionNumber = categorical(subT.InterventionNumber);
            formula = sprintf('%s ~ Condition + InterventionNumber + (1|Subject)', metricName);
            try
                lme3 = fitlme(subT, formula, 'FitMethod', 'REML');
                aov3 = anova(lme3, 'DFMethod', 'satterthwaite');
                rowIdx = find(strcmp(aov3.Term, 'Condition'), 1);
                sensOrderRows(end+1, :) = {metricName, bandName, aov3.FStat(rowIdx), aov3.DF1(rowIdx), aov3.DF2(rowIdx), aov3.pValue(rowIdx)}; %#ok<AGROW>
            catch ME
                sensitivityLog{end+1} = sprintf('Sensitivity 3 fit failed for %s / %s: %s', metricName, bandName, ME.message); %#ok<AGROW>
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
    if any(locB == 0); error('Phase4:SensitivityRowMismatch', 'Sensitivity 3 produced a band/metric row with no primary counterpart.'); end
    SensOrderT.SignificantFDR_primary = OmnibusT.SignificantFDR(locB);
    SensOrderT.ConclusionMatches = SensOrderT.SignificantFDR == SensOrderT.SignificantFDR_primary;
    nFlipped3 = sum(~SensOrderT.ConclusionMatches);
    sensitivityLog{end+1} = sprintf('Sensitivity 3 (order covariate added): %d/%d band-metric conclusions changed vs. the primary (order-free) analysis.', nFlipped3, height(SensOrderT));
end

%% Section 10: Sensitivity 4 -- paired-subset consistency check
pairedCheckRows = cell(0, 9); % Metric,Band,Contrast,N,LMM_Estimate,Paired_MeanDiff,PairedT_p,Wilcoxon_p
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
            title(sprintf('Residuals: %s (%s)', bandName, metricName));
            xlabel('Pearson residual'); ylabel('Count');
            subplot(1,2,2);
            qqplot(res);
            title('Normal Q-Q');
            outPng = fullfile(STATS_DIR, sprintf('diagnostic_residuals_%s_%s.png', bandName, metricName));
            saveas(fig, outPng);
            close(fig);
        catch ME
            sensitivityLog{end+1} = sprintf('Residual diagnostic failed for %s / %s: %s', metricName, bandName, ME.message); %#ok<AGROW>
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
if ~isempty(SensNoSub7T); writetable(SensNoSub7T, fullfile(STATS_DIR, 'phase4_sensitivity_nosub07.xlsx')); end
if ~isempty(SensOrderT);  writetable(SensOrderT,  fullfile(STATS_DIR, 'phase4_sensitivity_order.xlsx'));  end
writetable(PairedCheckT, fullfile(STATS_DIR, 'phase4_paired_subset_checks.xlsx'));
writetable(AllContrastsT, fullfile(STATS_DIR, 'phase4_all_contrasts_ungated.xlsx'));
save(fullfile(STATS_DIR, 'phase4_all_contrasts_ungated.mat'), 'AllContrastsT');

%% Section 13: Report
reportLines = {};
reportLines{end+1} = sprintf('RUN_MODE = %s | subjects excluded in this run: %s', RUN_MODE_LABEL, local_v3_list(RUN_MODE_EXCL_SUBJ));
reportLines{end+1} = repmat('=', 1, 60);
reportLines{end+1} = 'PHASE 4 -- PRIMARY E vs N vs R STATISTICS REPORT';
reportLines{end+1} = repmat('=', 1, 60);
reportLines{end+1} = '';
reportLines{end+1} = sprintf('Recordings verified: %d | Subjects: %d', nRecordings, nSubjects);
reportLines{end+1} = sprintf('Condition-pair coverage: %s', strjoin(cellfun(@(k) sprintf('%s=%d', strrep(k,'_',' & '), pairTypeCount(k)), pairKeys, 'UniformOutput', false), ', '));
reportLines{end+1} = sprintf('Condition x Order confound detected: %d', confoundFlag);
reportLines{end+1} = '';
reportLines{end+1} = 'Model: Linear mixed-effects, Condition (fixed, 3-level) + (1|Subject), REML, Satterthwaite DF.';
reportLines{end+1} = 'Analysis level: Global scalp band power (Phase 3 Structure C).';
reportLines{end+1} = 'Tiers: Primary=Alpha; Secondary=Delta/Theta/Beta; Exploratory=3 gamma bands; ExploratoryBeta=Beta_low 13-20 + Beta_high 20-30 (own family) (pre-specified, not data-driven).';
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
reportLines{end+1} = '--- [v2] ALL contrasts, every band/metric (UNGATED, descriptive only; sign = first minus second) ---';
for r = 1:height(AllContrastsT)
    reportLines{end+1} = sprintf('  [%s] %-25s %-6s est=%.4f (95%% CI %.4f to %.4f) p=%.4f pHolm(within model)=%.4f', ...
        AllContrastsT.Metric{r}, AllContrastsT.Band{r}, AllContrastsT.Contrast{r}, AllContrastsT.Estimate(r), ...
        AllContrastsT.CI_Low(r), AllContrastsT.CI_High(r), AllContrastsT.pValue(r), AllContrastsT.pHolm_withinModel(r)); %#ok<AGROW>
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
fid = fopen(fullfile(REPORT_DIR, sprintf('phase4_primary_statistics_report_%s.txt', RUN_MODE_LABEL)), 'w');
fprintf(fid, '%s\n', reportText);
fclose(fid);

%% ======================= LOCAL FUNCTIONS (all after top-level code) =======================

function n = local_subject_numeric(subjVal)
% Robustly extract a numeric subject id regardless of whether Subject is
% stored as numeric, categorical, char, or cellstr, and regardless of
% zero-padding (e.g. '7' vs '07').
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
% Benjamini-Hochberg step-up FDR procedure.
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
% Holm step-down procedure (controls family-wise error rate).
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
% column position/order).
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
            cols{j} = rows(:, j); % keep as Nx1 cell array of char -- never collapsed
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

function s = local_v3_list(v)
    if isempty(v); s = 'none'; else; s = strjoin(cellstr(num2str(v(:))), ', '); end
end
