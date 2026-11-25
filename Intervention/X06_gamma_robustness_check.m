%% ========================================================================
% X06_gamma_robustness_check.m
%
% Robustness check for the gamma-band (30-100 Hz) results: refits the X05 models after
% excluding the recordings flagged for a gamma band by the X04 extreme-value diagnostic
% (04_subject_level/diagnostic5_extreme_value_flags.csv) and compares the conclusions
% with X05.
% Inputs : 04_subject_level, 05_primary_stats.
% Outputs: 05_primary_stats/<RUN_MODE>/gamma_robustness_*.xlsx, 09_reports/.
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
GAMMA_DIR  = fullfile(STATS_DIR, 'gamma_robustness_check');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');

bandPowerFile   = fullfile(SUBJ_DIR, 'phase3_subject_level_bandpower.mat');
flagsFile       = fullfile(SUBJ_DIR, 'diagnostic5_extreme_value_flags.csv');
omnibusFile     = fullfile(STATS_DIR, 'phase4_omnibus_results.mat');
pairwiseFile    = fullfile(STATS_DIR, 'phase4_pairwise_contrasts.mat');

for f = {bandPowerFile, flagsFile, omnibusFile, pairwiseFile}
    if ~exist(f{1}, 'file')
        error('Phase4Gamma:MissingInput', 'Required input file not found: %s', f{1});
    end
end
if ~exist(GAMMA_DIR, 'dir');  mkdir(GAMMA_DIR);  end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end
if isempty(ver('stats'))
    error('Phase4Gamma:MissingToolbox', 'Statistics and Machine Learning Toolbox not detected.');
end

%% Section 1: Load Phase 3 band power + Phase 4 original results (read-only)
S = load(bandPowerFile, 'GlobalBandPowerT', 'PARAMS3');
GlobalBandPowerT = S.GlobalBandPowerT;
[GlobalBandPowerT, RUN_MODE_N_EXCL_ROWS, RUN_MODE_EXCL_SUBJ] = apply_run_mode_filter(GlobalBandPowerT, RUN_MODE, EXCLUDE_NOISY_SUBJECTS);
if RUN_MODE_N_EXCL_ROWS > 0; fprintf('RUN_MODE=%s: removed %d row(s) of subject(s) {%s} from GlobalBandPowerT.\n', RUN_MODE_LABEL, RUN_MODE_N_EXCL_ROWS, num2str(RUN_MODE_EXCL_SUBJ)); end

PARAMS3 = S.PARAMS3;
clear S;

S2 = load(omnibusFile, 'OmnibusT');
OriginalOmnibusT = S2.OmnibusT;
clear S2;

S3 = load(pairwiseFile, 'PairwiseT');
OriginalPairwiseT = S3.PairwiseT;
clear S3;

GAMMA_BANDS = {'Gamma_low_exploratory', 'Gamma_high1_exploratory', 'Gamma_high2_exploratory'};
allDefinedBands = fieldnames(PARAMS3.BandDefinitionsHz);
missingGammaBands = GAMMA_BANDS(~ismember(GAMMA_BANDS, allDefinedBands));
if ~isempty(missingGammaBands)
    error('Phase4Gamma:BandMismatch', 'Expected gamma band(s) not found in PARAMS3.BandDefinitionsHz: %s', strjoin(missingGammaBands, ', '));
end
METRICS = {'LogAbsolutePower', 'RelativePower'};
CONDITION_ORDER = {'E', 'N', 'R'};
REF_LEVEL = CONDITION_ORDER{1};
FDR_ALPHA = 0.05;

%% Section 2: Determine the exact gamma-flagged recordings from Phase 3's own diagnostic
flagsT = readtable(flagsFile);
gammaFlagRows = flagsT(ismember(flagsT.Band, GAMMA_BANDS), :);
if isempty(gammaFlagRows)
    msg = sprintf(['No gamma-band recordings were flagged in %s under the v2 rule ' ...
        '(log scale, within condition). There is nothing to exclude, so the gamma ' ...
        'robustness check is not needed; the Phase 4 gamma results stand as reported.'], flagsFile);
    fprintf('%s\n', msg);
    fid = fopen(fullfile(REPORT_DIR, sprintf('phase4_gamma_robustness_report_%s.txt', RUN_MODE_LABEL)), 'w');
    fprintf(fid, 'PHASE 4 SENSITIVITY ADDENDUM -- GAMMA ROBUSTNESS CHECK\n\n%s\n', msg);
    fclose(fid);
    return;
end
excludeSubj = local_subject_numeric(gammaFlagRows.Subject);
excludeInt  = gammaFlagRows.InterventionNumber;
excludeKeyTable = unique(table(excludeSubj, excludeInt, 'VariableNames', {'Subject', 'InterventionNumber'}), 'rows');

excludeCond = cell(height(excludeKeyTable), 1);
gSubjNum = local_subject_numeric(GlobalBandPowerT.Subject);
for k = 1:height(excludeKeyTable)
    hit = find(gSubjNum == excludeKeyTable.Subject(k) & GlobalBandPowerT.InterventionNumber == excludeKeyTable.InterventionNumber(k), 1);
    if isempty(hit); excludeCond{k} = '?'; else; excludeCond{k} = char(GlobalBandPowerT.Condition{hit}); end
end
excludeKeyTable.Condition = excludeCond;
exclCondCounts = struct('E', sum(strcmp(excludeCond, 'E')), 'N', sum(strcmp(excludeCond, 'N')), 'R', sum(strcmp(excludeCond, 'R')));
exclDominated = height(excludeKeyTable) >= 3 && ...
    max([exclCondCounts.E, exclCondCounts.N, exclCondCounts.R]) > (2/3) * height(excludeKeyTable);

fprintf('Gamma-related extreme-value flags found in Phase 3 diagnostic (%d rows, %d unique recordings):\n', height(gammaFlagRows), height(excludeKeyTable));
for i = 1:height(gammaFlagRows)
    fprintf('  Subject %g, Int%g, %s -- Band=%s, RobustZ=%.2f\n', local_subject_numeric(gammaFlagRows.Subject(i)), ...
        gammaFlagRows.InterventionNumber(i), char(gammaFlagRows.Condition(i)), char(gammaFlagRows.Band(i)), gammaFlagRows.RobustZ(i));
end
fprintf('\nUnique recordings excluded from THIS ROBUSTNESS CHECK ONLY (union across all 3 gamma bands, applied identically to every gamma band/metric fit below):\n');
for i = 1:height(excludeKeyTable)
    fprintf('  Subject %g, Int%g\n', excludeKeyTable.Subject(i), excludeKeyTable.InterventionNumber(i));
end

%% Section 3: Fit the 6 gamma-band x metric models on the exclusion-adjusted data
nBands = numel(GAMMA_BANDS);
nMetrics = numel(METRICS);
excludedOmnibusRows = cell(0, 6);
excludedLmeStore = cell(nBands, nMetrics);
excludedN = nan(nBands, nMetrics);
condCountsExcluded = cell(nBands, nMetrics);
integrityIssues = {};

for m = 1:nMetrics
    metricName = METRICS{m};
    for b = 1:nBands
        bandName = GAMMA_BANDS{b};
        subT = GlobalBandPowerT(strcmp(GlobalBandPowerT.Band, bandName), :);
        subjNum = local_subject_numeric(subT.Subject);

        excludeMask = false(height(subT), 1);
        for k = 1:height(excludeKeyTable)
            excludeMask = excludeMask | (subjNum == excludeKeyTable.Subject(k) & subT.InterventionNumber == excludeKeyTable.InterventionNumber(k));
        end
        subT_excl = subT(~excludeMask, :);
        excludedN(b, m) = sum(excludeMask);

        subT_excl.Subject = categorical(local_subject_numeric(subT_excl.Subject));
        subT_excl.Condition = categorical(cellstr(subT_excl.Condition), CONDITION_ORDER);
        condCountsExcluded{b, m} = countcats(subT_excl.Condition);

        formula = sprintf('%s ~ Condition + (1|Subject)', metricName);
        try
            lme = fitlme(subT_excl, formula, 'FitMethod', 'REML');
            aov = anova(lme, 'DFMethod', 'satterthwaite');
            rowIdx = find(strcmp(aov.Term, 'Condition'), 1);
            excludedOmnibusRows(end+1, :) = {metricName, bandName, aov.FStat(rowIdx), aov.DF1(rowIdx), aov.DF2(rowIdx), aov.pValue(rowIdx)}; %#ok<AGROW>
            excludedLmeStore{b, m} = lme;
        catch ME
            integrityIssues{end+1} = sprintf('fitlme FAILED (exclusion-adjusted) for %s / %s: %s', metricName, bandName, ME.message); %#ok<AGROW>
        end
    end
end

ExcludedOmnibusT = local_build_table(excludedOmnibusRows, [1 2], [3 4 5 6], {'Metric','Band','FStat','DF1','DF2','pValue'});

%% Section 4: FDR within the SAME per-metric gamma family (3 bands), matching Phase 4's exploratory tier
ExcludedOmnibusT.pFDR = nan(height(ExcludedOmnibusT), 1);
ExcludedOmnibusT.SignificantFDR = false(height(ExcludedOmnibusT), 1);
for m = 1:nMetrics
    metricName = METRICS{m};
    idx = strcmp(ExcludedOmnibusT.Metric, metricName);
    [adjP, reject] = local_fdr_bh(ExcludedOmnibusT.pValue(idx), FDR_ALPHA);
    ExcludedOmnibusT.pFDR(idx) = adjP;
    ExcludedOmnibusT.SignificantFDR(idx) = reject;
end

%% Section 5: Build the comparison table (Original vs. exclusion-adjusted)
compRows = cell(0, 12);
for r = 1:height(ExcludedOmnibusT)
    metricName = ExcludedOmnibusT.Metric{r};
    bandName = ExcludedOmnibusT.Band{r};
    origIdx = strcmp(OriginalOmnibusT.Metric, metricName) & strcmp(OriginalOmnibusT.Band, bandName);
    if sum(origIdx) ~= 1
        integrityIssues{end+1} = sprintf('Could not find a unique original Phase 4 omnibus row for %s / %s.', metricName, bandName); %#ok<AGROW>
        continue;
    end
    origF = OriginalOmnibusT.FStat(origIdx); origDF2 = OriginalOmnibusT.DF2(origIdx);
    origP = OriginalOmnibusT.pValue(origIdx); origPFDR = OriginalOmnibusT.pFDR(origIdx);
    origSig = OriginalOmnibusT.SignificantFDR(origIdx);

    exclF = ExcludedOmnibusT.FStat(r); exclDF2 = ExcludedOmnibusT.DF2(r);
    exclP = ExcludedOmnibusT.pValue(r); exclPFDR = ExcludedOmnibusT.pFDR(r);
    exclSig = ExcludedOmnibusT.SignificantFDR(r);

    compRows(end+1, :) = {metricName, bandName, origF, origDF2, origP, origPFDR, origSig, exclF, exclDF2, exclP, exclPFDR, exclSig}; %#ok<AGROW>
end
ComparisonT = local_build_table(compRows, [1 2], [3 4 5 6 7 8 9 10 11 12], ...
    {'Metric','Band','Orig_F','Orig_DF2','Orig_p','Orig_pFDR','Orig_Significant', ...
     'Excl_F','Excl_DF2','Excl_p','Excl_pFDR','Excl_Significant'});
ComparisonT.ConclusionChanged = ComparisonT.Orig_Significant ~= ComparisonT.Excl_Significant;
ComparisonT.N_Excluded = nan(height(ComparisonT), 1);
for r = 1:height(ComparisonT)
    b = find(strcmp(GAMMA_BANDS, ComparisonT.Band{r}));
    m = find(strcmp(METRICS, ComparisonT.Metric{r}));
    ComparisonT.N_Excluded(r) = excludedN(b, m);
end
% Interpretation column, computed cleanly against the FINAL table rows
% (never against the pre-skip loop index above, to avoid any risk of
% misalignment if a band/metric combination were ever missing).
ComparisonT.Interpretation = cell(height(ComparisonT), 1);
for r = 1:height(ComparisonT)
    if ComparisonT.ConclusionChanged(r)
        ComparisonT.Interpretation{r} = 'CONCLUSION CHANGED after excluding gamma-flagged recordings -- treat with caution, see jackknife.';
    elseif ComparisonT.Excl_Significant(r)
        ComparisonT.Interpretation{r} = sprintf('Remains significant after excluding %d flagged recording(s); conclusion unchanged.', ComparisonT.N_Excluded(r));
    else
        ComparisonT.Interpretation{r} = sprintf('Remains non-significant after excluding %d flagged recording(s); conclusion unchanged.', ComparisonT.N_Excluded(r));
    end
end

fprintf('\n--- Comparison: original vs. gamma-flagged-excluded omnibus results ---\n');
disp(ComparisonT);

%% Section 6: Pairwise contrasts for exclusion-adjusted models that remain significant
CONTRASTS = {'E','N','E - N'; 'E','R','E - R'; 'N','R','N - R'};   % {level2, level1, label}; estimate = level2 - level1
exclPairwiseRows = cell(0, 11);
sigMask = ExcludedOmnibusT.SignificantFDR;
for r = find(sigMask)'
    metricName = ExcludedOmnibusT.Metric{r};
    bandName = ExcludedOmnibusT.Band{r};
    b = find(strcmp(GAMMA_BANDS, bandName));
    m = find(strcmp(METRICS, metricName));
    lme = excludedLmeStore{b, m};
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
        contrastP(c) = p;
        contrastData(c, :) = {label, est, se, df2c, tstat, p, est - tcrit*se, est + tcrit*se};
    end
    [adjP, ~] = local_holm_bonferroni(contrastP, FDR_ALPHA);
    for c = 1:3
        exclPairwiseRows(end+1, :) = [{metricName, bandName}, contrastData(c, 1:6), {adjP(c)}, contrastData(c, 7:8)]; %#ok<AGROW>
    end
end
ExcludedPairwiseT = local_build_table(exclPairwiseRows, [1 2 3], [4 5 6 7 8 9 10 11], ...
    {'Metric','Band','Contrast','Estimate','SE','DF','tStat','pValue','pHolm','CI_Low','CI_High'});

fprintf('\n--- Pairwise contrasts, exclusion-adjusted models (only where still significant) ---\n');
disp(ExcludedPairwiseT);

OriginalGammaPairwiseT = OriginalPairwiseT(ismember(OriginalPairwiseT.Band, GAMMA_BANDS), :);
fprintf('\n--- Original (Phase 4, no exclusion) pairwise contrasts for the gamma bands, for reference ---\n');
disp(OriginalGammaPairwiseT);

%% Section 7: Diagnostics -- residual histogram + QQ + residual-vs-fitted, per gamma band x metric
for m = 1:nMetrics
    metricName = METRICS{m};
    for b = 1:nBands
        bandName = GAMMA_BANDS{b};
        lme = excludedLmeStore{b, m};
        if isempty(lme); continue; end
        try
            res = residuals(lme, 'ResidualType', 'Pearson');
            fittedVals = fitted(lme);
            fig = figure('Visible', 'off', 'Position', [100 100 1300 400]);
            subplot(1,3,1);
            histogram(res, 12);
            title(sprintf('Residuals: %s (%s)', bandName, metricName), 'Interpreter', 'none');
            xlabel('Pearson residual'); ylabel('Count');
            subplot(1,3,2);
            qqplot(res);
            title('Normal Q-Q (n small -- do not over-interpret)');
            subplot(1,3,3);
            scatter(fittedVals, res, 25, 'filled');
            hold on; yline(0, 'k--'); hold off;
            xlabel('Fitted value'); ylabel('Pearson residual');
            title('Residual vs. Fitted');
            outPng = fullfile(GAMMA_DIR, sprintf('diagnostic_gamma_%s_%s.png', bandName, metricName));
            saveas(fig, outPng);
            close(fig);
        catch ME
            integrityIssues{end+1} = sprintf('Diagnostic plot failed for %s / %s: %s', metricName, bandName, ME.message); %#ok<AGROW>
        end
    end
end

%% Section 8: Leave-one-subject-out jackknife on the exclusion-adjusted dataset
jackRows = cell(0, 5); % Metric, Band, SubjectDropped, pValue_withoutSubject, Delta_from_ExclModel
mostInfluential = struct();
for m = 1:nMetrics
    metricName = METRICS{m};
    for b = 1:nBands
        bandName = GAMMA_BANDS{b};
        subT = GlobalBandPowerT(strcmp(GlobalBandPowerT.Band, bandName), :);
        subjNum = local_subject_numeric(subT.Subject);
        excludeMask = false(height(subT), 1);
        for k = 1:height(excludeKeyTable)
            excludeMask = excludeMask | (subjNum == excludeKeyTable.Subject(k) & subT.InterventionNumber == excludeKeyTable.InterventionNumber(k));
        end
        subT_excl = subT(~excludeMask, :);
        subjNum_excl = local_subject_numeric(subT_excl.Subject);
        remainingSubj = unique(subjNum_excl);

        pFullExcl = ExcludedOmnibusT.pValue(strcmp(ExcludedOmnibusT.Metric, metricName) & strcmp(ExcludedOmnibusT.Band, bandName));
        if numel(pFullExcl) ~= 1
            integrityIssues{end+1} = sprintf('Skipping jackknife for %s / %s: no unique exclusion-adjusted omnibus p-value available (model likely failed to fit in Section 3).', metricName, bandName); %#ok<AGROW>
            key = sprintf('%s_%s', metricName, bandName);
            mostInfluential.(key) = struct('Subject', NaN, 'pWithoutSubject', NaN, 'pFullExcl', NaN, 'Delta', NaN);
            continue;
        end
        maxIncrease = -Inf; maxIncreaseSubj = NaN; maxIncreaseP = NaN;

        for si = 1:numel(remainingSubj)
            thisSubj = remainingSubj(si);
            jkT = subT_excl(subjNum_excl ~= thisSubj, :);
            jkT.Subject = categorical(local_subject_numeric(jkT.Subject));
            jkT.Condition = categorical(cellstr(jkT.Condition), CONDITION_ORDER);
            formula = sprintf('%s ~ Condition + (1|Subject)', metricName);
            try
                lmeJK = fitlme(jkT, formula, 'FitMethod', 'REML');
                aovJK = anova(lmeJK, 'DFMethod', 'satterthwaite');
                rowIdx = find(strcmp(aovJK.Term, 'Condition'), 1);
                pJK = aovJK.pValue(rowIdx);
            catch
                pJK = NaN;
            end
            jackRows(end+1, :) = {metricName, bandName, thisSubj, pJK, pJK - pFullExcl}; %#ok<AGROW>
            if ~isnan(pJK) && (pJK - pFullExcl) > maxIncrease
                maxIncrease = pJK - pFullExcl;
                maxIncreaseSubj = thisSubj;
                maxIncreaseP = pJK;
            end
        end
        key = sprintf('%s_%s', metricName, bandName);
        mostInfluential.(key) = struct('Subject', maxIncreaseSubj, 'pWithoutSubject', maxIncreaseP, ...
            'pFullExcl', pFullExcl, 'Delta', maxIncrease);
    end
end
JackknifeT = local_build_table(jackRows, [1 2], [3 4 5], {'Metric','Band','SubjectDropped','pValue_withoutSubject','Delta_from_ExclModel'});

%% Section 9: Save outputs
writetable(ComparisonT, fullfile(GAMMA_DIR, 'gamma_robustness_comparison.xlsx'));
save(fullfile(GAMMA_DIR, 'gamma_robustness_comparison.mat'), 'ComparisonT', 'ExcludedOmnibusT', 'OriginalOmnibusT');
if ~isempty(ExcludedPairwiseT)
    writetable(ExcludedPairwiseT, fullfile(GAMMA_DIR, 'gamma_robustness_pairwise.xlsx'));
end
writetable(JackknifeT, fullfile(GAMMA_DIR, 'gamma_robustness_jackknife.xlsx'));

PARAMS_GAMMA = struct();
PARAMS_GAMMA.Purpose = 'Phase 4 sensitivity addendum: robustness of exploratory gamma-band findings to Phase-3-flagged recordings.';
PARAMS_GAMMA.ExclusionSet = excludeKeyTable;
PARAMS_GAMMA.ExclusionRationale = ['Union of all recordings flagged for ANY of the 3 gamma bands in ' ...
    'diagnostic5_extreme_value_flags.csv, applied identically to all gamma band/metric fits. ' ...
    'Not applied to Alpha, Delta, Theta, or Beta -- only used for this gamma-focused check.'];
PARAMS_GAMMA.Model = '(unchanged from Phase 4) Condition (fixed, 3-level, ref=E) + (1|Subject), REML, Satterthwaite DF';
PARAMS_GAMMA.FDRFamilyScheme = 'Same as Phase 4 exploratory tier: BH-FDR across the 3 gamma bands, separately per metric.';
PARAMS_GAMMA.FDRAlpha = FDR_ALPHA;
PARAMS_GAMMA.JackknifeNote = 'Leave-one-subject-out on the exclusion-adjusted data; descriptive proxy for influence, not a formal test.';
PARAMS_GAMMA.RunDateTime = datestr(now); %#ok<TNOW1,DATST>
save(fullfile(GAMMA_DIR, 'gamma_robustness_parameters.mat'), 'PARAMS_GAMMA');

%% Section 10: Report
reportLines = {};
reportLines{end+1} = sprintf('RUN_MODE = %s | subjects excluded in this run: %s', RUN_MODE_LABEL, local_v3_list(RUN_MODE_EXCL_SUBJ));
reportLines{end+1} = repmat('=', 1, 60);
reportLines{end+1} = 'PHASE 4 SENSITIVITY ADDENDUM -- GAMMA ROBUSTNESS CHECK';
reportLines{end+1} = '(This is NOT a new primary analysis and NOT Phase 5.)';
reportLines{end+1} = repmat('=', 1, 60);
reportLines{end+1} = '';
reportLines{end+1} = sprintf('Gamma-flagged recordings excluded from THIS CHECK ONLY (%d unique recordings):', height(excludeKeyTable));
for i = 1:height(excludeKeyTable)
    reportLines{end+1} = sprintf('  - Subject %g, InterventionNumber %g, Condition %s', excludeKeyTable.Subject(i), ...
        excludeKeyTable.InterventionNumber(i), excludeKeyTable.Condition{i}); %#ok<AGROW>
end
reportLines{end+1} = sprintf('Exclusion set by condition: E=%d, N=%d, R=%d', exclCondCounts.E, exclCondCounts.N, exclCondCounts.R);
if exclDominated
    reportLines{end+1} = ['  *** WARNING: one condition makes up more than 2/3 of the exclusion set. Excluding these ' ...
        'recordings changes that condition''s mean more than the others, so a change in a contrast after ' ...
        'exclusion is partly produced by the exclusion itself -- interpret with caution.'];
end
reportLines{end+1} = '';
reportLines{end+1} = 'Comparison table: Band x Metric | Original | Exclusion-adjusted | Changed?';
for r = 1:height(ComparisonT)
    reportLines{end+1} = sprintf(['  [%s] %-25s N_excluded=%d | orig F(2,%.1f)=%.3f p=%.4f pFDR=%.4f sig=%d | ' ...
        'excl F(2,%.1f)=%.3f p=%.4f pFDR=%.4f sig=%d | CHANGED=%d'], ...
        ComparisonT.Metric{r}, ComparisonT.Band{r}, ComparisonT.N_Excluded(r), ...
        ComparisonT.Orig_DF2(r), ComparisonT.Orig_F(r), ComparisonT.Orig_p(r), ComparisonT.Orig_pFDR(r), ComparisonT.Orig_Significant(r), ...
        ComparisonT.Excl_DF2(r), ComparisonT.Excl_F(r), ComparisonT.Excl_p(r), ComparisonT.Excl_pFDR(r), ComparisonT.Excl_Significant(r), ...
        ComparisonT.ConclusionChanged(r)); %#ok<AGROW>
end
reportLines{end+1} = '';
reportLines{end+1} = 'Pairwise contrasts, exclusion-adjusted models (computed only where still FDR-significant):';
if isempty(ExcludedPairwiseT)
    reportLines{end+1} = '  (none survived FDR correction post-exclusion)';
else
    for r = 1:height(ExcludedPairwiseT)
        reportLines{end+1} = sprintf('  [%s] %-25s %-8s est=%.4f t=%.3f df=%.1f p=%.4f pHolm=%.4f', ...
            ExcludedPairwiseT.Metric{r}, ExcludedPairwiseT.Band{r}, ExcludedPairwiseT.Contrast{r}, ...
            ExcludedPairwiseT.Estimate(r), ExcludedPairwiseT.tStat(r), ExcludedPairwiseT.DF(r), ...
            ExcludedPairwiseT.pValue(r), ExcludedPairwiseT.pHolm(r)); %#ok<AGROW>
    end
end
reportLines{end+1} = '';
reportLines{end+1} = 'Original (Phase 4, no exclusion) pairwise contrasts for the gamma bands, for direct comparison:';
for r = 1:height(OriginalGammaPairwiseT)
    reportLines{end+1} = sprintf('  [%s] %-25s %-8s est=%.4f p=%.4f pHolm=%.4f', ...
        OriginalGammaPairwiseT.Metric{r}, OriginalGammaPairwiseT.Band{r}, OriginalGammaPairwiseT.Contrast{r}, ...
        OriginalGammaPairwiseT.Estimate(r), OriginalGammaPairwiseT.pValue(r), OriginalGammaPairwiseT.pHolm(r)); %#ok<AGROW>
end
reportLines{end+1} = '';
reportLines{end+1} = 'Leave-one-subject-out jackknife summary (descriptive; NOT a formal influence test):';
for m = 1:nMetrics
    for b = 1:nBands
        key = sprintf('%s_%s', METRICS{m}, GAMMA_BANDS{b});
        mi = mostInfluential.(key);
        reportLines{end+1} = sprintf(['  [%s] %-25s exclusion-adjusted p=%.4f | largest p-increase on dropping Subject %g: ' ...
            'p=%.4f (Delta=+%.4f)'], METRICS{m}, GAMMA_BANDS{b}, mi.pFullExcl, mi.Subject, mi.pWithoutSubject, mi.Delta); %#ok<AGROW>
    end
end
reportLines{end+1} = '';
if ~isempty(integrityIssues)
    reportLines{end+1} = '--- Integrity issues ---';
    for i = 1:numel(integrityIssues)
        reportLines{end+1} = ['  - ' integrityIssues{i}]; %#ok<AGROW>
    end
else
    reportLines{end+1} = 'Integrity issues: none.';
end
reportLines{end+1} = '';
reportLines{end+1} = '--- Answers to the four specific questions ---';
nChanged = sum(ComparisonT.ConclusionChanged);
reportLines{end+1} = sprintf('1. Do the gamma findings remain significant after exclusion? %d/%d band-metric combinations kept their original conclusion; %d changed. See table above for which.', ...
    height(ComparisonT) - nChanged, height(ComparisonT), nChanged);
reportLines{end+1} = '2. Does the E-vs-N (E - N) conclusion change? Per band/metric, original vs. exclusion-adjusted E - N:';
for r = 1:height(OriginalGammaPairwiseT)
    if ~strcmp(OriginalGammaPairwiseT.Contrast{r}, 'E - N')
        continue;
    end
    metricName = OriginalGammaPairwiseT.Metric{r};
    bandName = OriginalGammaPairwiseT.Band{r};
    origP_EN = OriginalGammaPairwiseT.pValue(r);
    if OriginalGammaPairwiseT.pHolm(r) <= FDR_ALPHA; origSigStr = 'significant'; else; origSigStr = 'n.s.'; end
    exclIdx = strcmp(ExcludedPairwiseT.Metric, metricName) & strcmp(ExcludedPairwiseT.Band, bandName) & strcmp(ExcludedPairwiseT.Contrast, 'E - N');
    if any(exclIdx)
        if ExcludedPairwiseT.pHolm(exclIdx) <= FDR_ALPHA; exclSigStr = 'significant'; else; exclSigStr = 'n.s.'; end
        reportLines{end+1} = sprintf('   [%s] %-25s original E-N est=%.4f p=%.4f (%s) -> exclusion-adjusted E-N est=%.4f p=%.4f (%s)', ...
            metricName, bandName, OriginalGammaPairwiseT.Estimate(r), origP_EN, origSigStr, ...
            ExcludedPairwiseT.Estimate(exclIdx), ExcludedPairwiseT.pValue(exclIdx), exclSigStr); %#ok<AGROW>
    else
        reportLines{end+1} = sprintf(['   [%s] %-25s original E-N p=%.4f (%s) -> exclusion-adjusted: omnibus no longer FDR-significant, ' ...
            'so E-N was not (re)tested (same gatekeeping rule as the primary analysis).'], ...
            metricName, bandName, origP_EN, origSigStr); %#ok<AGROW>
    end
end
reportLines{end+1} = '3. Evidence of a small number of recordings driving the effect: see the jackknife summary above -- the single largest p-value increase from dropping one additional subject is reported per band/metric; a large jump there (especially one that would flip significance on its own) would indicate fragility.';
reportLines{end+1} = '4. Should the gamma findings remain exploratory? This script does not decide that -- it only supplies the robustness evidence above; the exploratory-tier status itself is a pre-specified, primary-analysis decision this script does not touch.';
reportLines{end+1} = '';
reportLines{end+1} = 'Output files written to:';
reportLines{end+1} = ['  ' GAMMA_DIR];
reportLines{end+1} = repmat('=', 1, 60);

reportText = strjoin(reportLines, '\n');
fprintf('%s\n', reportText);
fid = fopen(fullfile(REPORT_DIR, sprintf('phase4_gamma_robustness_report_%s.txt', RUN_MODE_LABEL)), 'w');
fprintf(fid, '%s\n', reportText);
fclose(fid);

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
            error('Phase4Gamma:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name1);
        end
        cRow(idx1) = cRow(idx1) - 1;
    end
    if ~strcmp(level2, refLevel)
        name2 = sprintf('%s_%s', varName, level2);
        idx2 = find(strcmp(coeffNames, name2));
        if isempty(idx2)
            error('Phase4Gamma:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name2);
        end
        cRow(idx2) = cRow(idx2) + 1;
    end
end

function T = local_build_table(rows, textColIdx, numColIdx, varNames)
% See 05_primary_statistics.m for the full rationale: table() is used
% directly (never cell2table) so text columns are never silently
% reinterpreted.
    nCols = numel(varNames);
    if size(rows, 2) ~= nCols
        error('Phase4Gamma:BuildTableColumnMismatch', ...
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
            error('Phase4Gamma:BuildTableColumnMismatch', 'Column %d not classified as text or numeric.', j);
        end
    end
    T = table(cols{:}, 'VariableNames', varNames);
end

function s = local_v3_list(v)
    if isempty(v); s = 'none'; else; s = strjoin(cellstr(num2str(v(:))), ', '); end
end
