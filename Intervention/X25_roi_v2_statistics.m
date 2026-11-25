%% ========================================================================
% X25_roi_v2_statistics.m
%
% Planned contrasts and regional statistics for the X24 regions.
% Planned contrasts (two-sided, Holm across the four tests = 2 hypotheses x absolute and
% relative power): frontal-midline theta E - N and occipital alpha E - N. Each planned
% contrast is also reported with an order-adjusted model, a within-person paired check and
% leave-one-out. Regional omnibus tests are FDR corrected; contrasts Holm corrected.
% Inputs : 02_metadata, 04_subject_level, 24_roi_v2_definitions.
% Outputs: 25_roi_v2_stats/<RUN_MODE>/ (roi_v2_confirmatory.csv, roi_v2_all_contrasts.csv,
%          roi_v2_emm.csv), 09_reports/.
% Run twice, with RUN_MODE = 'full_sample' (N = 30) and RUN_MODE = 'sensitivity_excl_noisy'
% (cleaned sample, N = 28; excluded subjects listed in EXCLUSION_CONFIG_intervention.m).
% Outputs are written to <OUT_DIR>/<folder>/<RUN_MODE>/ and reports to <OUT_DIR>/09_reports/.
% ========================================================================

clear; clc;

%% 0. SETTINGS
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
ROI2_DIR   = fullfile(ANALYSIS_ROOT, '24_roi_v2_definitions');
OUT_DIR    = fullfile(ANALYSIS_ROOT, '25_roi_v2_stats', RUN_MODE);
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');
bandPowerFile = fullfile(SUBJ_DIR, 'phase3_subject_level_bandpower.mat');
roi2File      = fullfile(ROI2_DIR, 'roi_v2_definitions.mat');
metaFile      = fullfile(ANALYSIS_ROOT, '02_metadata', 'intervention_metadata.mat');

ALPHA = 0.05;
CONDITION_ORDER = {'E', 'N', 'R'};
REF_LEVEL = 'E';
CONTRASTS = {'E','N','E - N'; 'E','R','E - R'; 'N','R','N - R'};   % {lvl2, lvl1, label}; label = lvl2 - lvl1
METRICS = {'LogAbsolutePower', 'RelativePower'};

for f = {bandPowerFile, roi2File, metaFile}
    if ~exist(f{1}, 'file')
        error('Phase6B:MissingInput', 'Required input not found: %s', f{1});
    end
end
if isempty(ver('stats'))
    error('Phase6B:MissingToolbox', 'Statistics and Machine Learning Toolbox is required.');
end
if ~exist(OUT_DIR, 'dir');    mkdir(OUT_DIR);    end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

fprintf('========================================================\n');
fprintf('PHASE 6B -- ROI v2 STATISTICS\n');
fprintf('========================================================\n');

%% 1. LOAD
S3 = load(bandPowerFile, 'BandPowerT', 'PARAMS3');
BandPowerT = S3.BandPowerT; PARAMS3 = S3.PARAMS3; clear S3;
[BandPowerT, RUN_MODE_N_EXCL_ROWS, RUN_MODE_EXCL_SUBJ] = apply_run_mode_filter(BandPowerT, RUN_MODE, EXCLUDE_NOISY_SUBJECTS);
if RUN_MODE_N_EXCL_ROWS > 0; fprintf('RUN_MODE=%s: removed %d row(s) of subject(s) {%s} from BandPowerT.\n', RUN_MODE_LABEL, RUN_MODE_N_EXCL_ROWS, num2str(RUN_MODE_EXCL_SUBJ)); end

S11 = load(roi2File, 'ROI2_DEFINITION', 'ROI2_NAMES', 'CONFIRMATORY', 'CONFIRMATORY_METRICS', 'PARAMS11');
ROI2_DEFINITION = S11.ROI2_DEFINITION; ROI2_NAMES = S11.ROI2_NAMES;
CONFIRMATORY = S11.CONFIRMATORY; CONF_METRICS = S11.CONFIRMATORY_METRICS; PARAMS11 = S11.PARAMS11; clear S11;
SM = load(metaFile, 'T'); MetaT = SM.T; clear SM;
[MetaT, RUN_MODE_N_EXCL_ROWS, RUN_MODE_EXCL_SUBJ] = apply_run_mode_filter(MetaT, RUN_MODE, EXCLUDE_NOISY_SUBJECTS);
if RUN_MODE_N_EXCL_ROWS > 0; fprintf('RUN_MODE=%s: removed %d row(s) of subject(s) {%s} from MetaT.\n', RUN_MODE_LABEL, RUN_MODE_N_EXCL_ROWS, num2str(RUN_MODE_EXCL_SUBJ)); end


ALL_BANDS = fieldnames(PARAMS3.BandDefinitionsHz);
nROI = numel(ROI2_NAMES); nBands = numel(ALL_BANDS); nMetrics = numel(METRICS);
BandPowerT.ChannelUC = upper(BandPowerT.Channel);
issues = {};

%% 2. VALIDATION (counts derived from Phase 1, never hard-coded)
fails = 0;
comboT = unique(BandPowerT(:, {'Subject','InterventionNumber','Condition'}), 'rows');
nRec = height(comboT);
if nRec ~= height(MetaT)
    fprintf('  *** FAIL: %d recordings in BandPowerT, Phase 1 has %d.\n', nRec, height(MetaT)); fails = fails + 1;
end
for c = 1:3
    nHere = sum(strcmp(comboT.Condition, CONDITION_ORDER{c}));
    nMeta = sum(strcmp(MetaT.Condition, CONDITION_ORDER{c}));
    if nHere ~= nMeta
        fprintf('  *** FAIL: condition %s has %d recordings here, %d in Phase 1.\n', CONDITION_ORDER{c}, nHere, nMeta); fails = fails + 1;
    end
end
subjNum = local_subject_numeric(comboT.Subject);
uSubj = unique(subjNum);
pairCount = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:numel(uSubj)
    cs = sort(comboT.Condition(subjNum == uSubj(i)));
    if numel(cs) ~= 2 || strcmp(cs{1}, cs{2})
        fprintf('  *** FAIL: subject %g does not have exactly 2 different conditions.\n', uSubj(i)); fails = fails + 1;
        continue;
    end
    key = strjoin(cs', '&');
    if isKey(pairCount, key); pairCount(key) = pairCount(key) + 1; else; pairCount(key) = 1; end
end
allRoiChans = {};
for r = 1:nROI
    labs = upper(ROI2_DEFINITION.(ROI2_NAMES{r}));
    allRoiChans = [allRoiChans, labs]; %#ok<AGROW>
    missing = labs(~ismember(labs, unique(BandPowerT.ChannelUC)));
    if ~isempty(missing)
        fprintf('  *** FAIL: ROI %s channel(s) missing from BandPowerT: %s\n', ROI2_NAMES{r}, strjoin(missing, ', ')); fails = fails + 1;
    end
end
for i = 1:height(CONFIRMATORY)
    if ~ismember(CONFIRMATORY.Band{i}, ALL_BANDS)
        fprintf('  *** FAIL: confirmatory band %s not in PARAMS3.\n', CONFIRMATORY.Band{i}); fails = fails + 1;
    end
end
pairKeys = keys(pairCount);
pairStr = strjoin(cellfun(@(k) sprintf('%s=%d', k, pairCount(k)), pairKeys, 'UniformOutput', false), ', ');
fprintf('  %d recordings, %d subjects, condition pairs: %s\n', nRec, numel(uSubj), pairStr);
if fails > 0
    error('Phase6B:ValidationFailed', '%d validation failure(s). Nothing fitted.', fails);
end

%% 3. ROI AGGREGATION
roiRows = cell(0, 7);
for i = 1:nRec
    recMask = strcmp(BandPowerT.Subject, comboT.Subject{i}) & BandPowerT.InterventionNumber == comboT.InterventionNumber(i);
    for r = 1:nROI
        chans = upper(ROI2_DEFINITION.(ROI2_NAMES{r}));
        chanMask = ismember(BandPowerT.ChannelUC, chans);
        for b = 1:nBands
            rows = recMask & chanMask & strcmp(BandPowerT.Band, ALL_BANDS{b});
            if sum(rows) ~= numel(chans)
                issues{end+1} = sprintf('Sub%s Int%g ROI %s %s: %d/%d channels found.', comboT.Subject{i}, ...
                    comboT.InterventionNumber(i), ROI2_NAMES{r}, ALL_BANDS{b}, sum(rows), numel(chans)); %#ok<SAGROW>
            end
            roiRows(end+1, :) = {comboT.Subject{i}, comboT.InterventionNumber(i), comboT.Condition{i}, ROI2_NAMES{r}, ...
                ALL_BANDS{b}, mean(BandPowerT.AbsolutePower(rows)), mean(BandPowerT.RelativePower(rows))}; %#ok<SAGROW>
        end
    end
end
Roi2T = local_build_table(roiRows, [1 3 4 5], [2 6 7], {'Subject','InterventionNumber','Condition','ROI','Band','AbsolutePower','RelativePower'});
Roi2T.LogAbsolutePower = log(Roi2T.AbsolutePower);
Roi2T.LogAbsolutePower(Roi2T.AbsolutePower <= 0) = NaN;
fprintf('  ROI aggregation: %d rows (%d recordings x %d ROIs x %d bands).\n', height(Roi2T), nRec, nROI, nBands);

%% 4. FIT ALL MODELS (ROI x band x metric)
lmeStore = cell(nROI, nBands, nMetrics);
omniRows = cell(0, 7); emmRows = cell(0, 6); allRows = cell(0, 13);
for r = 1:nROI
    for b = 1:nBands
        for m = 1:nMetrics
            subT = local_prepare(Roi2T, ROI2_NAMES{r}, ALL_BANDS{b}, CONDITION_ORDER);
            try
                lme = fitlme(subT, sprintf('%s ~ Condition + (1|Subject)', METRICS{m}), 'FitMethod', 'REML');
            catch ME
                issues{end+1} = sprintf('fitlme failed: %s %s %s: %s', ROI2_NAMES{r}, ALL_BANDS{b}, METRICS{m}, ME.message); %#ok<SAGROW>
                continue;
            end
            lmeStore{r, b, m} = lme;
            aov = anova(lme, 'DFMethod', 'satterthwaite');
            ri = find(strcmp(aov.Term, 'Condition'), 1);
            omniRows(end+1, :) = {ROI2_NAMES{r}, METRICS{m}, ALL_BANDS{b}, aov.FStat(ri), aov.DF1(ri), aov.DF2(ri), aov.pValue(ri)}; %#ok<SAGROW>
            emm = local_emm(lme);
            emmRows(end+1, :) = {ROI2_NAMES{r}, METRICS{m}, ALL_BANDS{b}, emm(1), emm(2), emm(3)}; %#ok<SAGROW>
            cres = local_contrasts(lme, CONTRASTS, REF_LEVEL);
            adj = local_holm(cres(:, 5), ALPHA);
            for c = 1:3
                allRows(end+1, :) = {ROI2_NAMES{r}, METRICS{m}, ALL_BANDS{b}, CONTRASTS{c, 3}, cres(c, 1), cres(c, 2), ...
                    cres(c, 3), cres(c, 4), cres(c, 5), adj(c), cres(c, 6), cres(c, 7), NaN}; %#ok<SAGROW>
            end
        end
    end
end
OmnibusT = local_build_table(omniRows, [1 2 3], 4:7, {'ROI','Metric','Band','FStat','DF1','DF2','pValue'});
EmmT = local_build_table(emmRows, [1 2 3], 4:6, {'ROI','Metric','Band','EMM_E','EMM_N','EMM_R'});
AllContrastsT = local_build_table(allRows, [1 2 3 4], 5:13, ...
    {'ROI','Metric','Band','Contrast','Estimate','SE','DF','tStat','pValue','pHolm_withinModel','CI_Low','CI_High','OmnibusPassedFDR'});

%% 5. SECONDARY: BH-FDR of the omnibus test within each band x metric (across ROIs)
OmnibusT.pFDR = nan(height(OmnibusT), 1);
OmnibusT.SignificantFDR = false(height(OmnibusT), 1);
for b = 1:nBands
    for m = 1:nMetrics
        idx = strcmp(OmnibusT.Band, ALL_BANDS{b}) & strcmp(OmnibusT.Metric, METRICS{m});
        if ~any(idx); continue; end
        [adjP, rej] = local_fdr_bh(OmnibusT.pValue(idx), ALPHA);
        OmnibusT.pFDR(idx) = adjP; OmnibusT.SignificantFDR(idx) = rej;
    end
end
for k = 1:height(AllContrastsT)
    o = strcmp(OmnibusT.ROI, AllContrastsT.ROI{k}) & strcmp(OmnibusT.Band, AllContrastsT.Band{k}) & strcmp(OmnibusT.Metric, AllContrastsT.Metric{k});
    AllContrastsT.OmnibusPassedFDR(k) = double(any(OmnibusT.SignificantFDR(o)));
end

%% 6. CONFIRMATORY PLANNED CONTRASTS (+ robustness checks)
confRows = cell(0, 19);
confWilcoxon = []; confLOOsubj = []; confEstOrder = [];
for i = 1:height(CONFIRMATORY)
    roiName = CONFIRMATORY.ROI{i}; bandName = CONFIRMATORY.Band{i}; label = CONFIRMATORY.Contrast{i};
    r = find(strcmp(ROI2_NAMES, roiName)); b = find(strcmp(ALL_BANDS, bandName));
    ci = find(strcmp(CONTRASTS(:, 3), label));
    for mm = 1:numel(CONF_METRICS)
        metricName = CONF_METRICS{mm}; m = find(strcmp(METRICS, metricName));
        lme = lmeStore{r, b, m};
        if isempty(lme)
            issues{end+1} = sprintf('Confirmatory %s %s: model missing.', CONFIRMATORY.Hypothesis{i}, metricName); %#ok<SAGROW>
            continue;
        end
        cres = local_contrasts(lme, CONTRASTS(ci, :), REF_LEVEL);   % est se df t p lo hi

        % (a) order-adjusted
        subT = local_prepare(Roi2T, roiName, bandName, CONDITION_ORDER);
        subT.InterventionNumber = categorical(subT.InterventionNumber);
        try
            lmeO = fitlme(subT, sprintf('%s ~ Condition + InterventionNumber + (1|Subject)', metricName), 'FitMethod', 'REML');
            cO = local_contrasts(lmeO, CONTRASTS(ci, :), REF_LEVEL);
            estO = cO(1); pO = cO(5);
        catch
            estO = NaN; pO = NaN;
        end

        % (b) paired subset (subjects who did both conditions of this contrast)
        lvl2 = CONTRASTS{ci, 1}; lvl1 = CONTRASTS{ci, 2};
        [x2, x1] = local_paired(subT, metricName, lvl2, lvl1);
        d = x2 - x1;
        if numel(d) >= 3
            [~, pT] = ttest(d); pW = signrank(d); dz = mean(d) / std(d);
        else
            pT = NaN; pW = NaN; dz = NaN;
        end

        % (c) leave-one-subject-out: worst-case (largest) p
        sNum = local_subject_numeric(cellstr(subT.Subject));
        uS = unique(sNum); pLOO = nan(numel(uS), 1);
        for s = 1:numel(uS)
            jk = subT(sNum ~= uS(s), :);
            jk.Subject = removecats(jk.Subject);
            try
                lJ = fitlme(jk, sprintf('%s ~ Condition + (1|Subject)', metricName), 'FitMethod', 'REML');
                cJ = local_contrasts(lJ, CONTRASTS(ci, :), REF_LEVEL); pLOO(s) = cJ(5);
            catch
            end
        end
        [pLOOmax, iMax] = max(pLOO);

        confRows(end+1, :) = {CONFIRMATORY.Hypothesis{i}, roiName, bandName, metricName, label, ...
            cres(1), cres(2), cres(3), cres(4), cres(5), NaN, cres(6), cres(7), ...
            pO, numel(d), mean(d), dz, pT, pLOOmax}; %#ok<SAGROW>
        confWilcoxon(end+1, 1) = pW; %#ok<SAGROW>
        confLOOsubj(end+1, 1) = uS(iMax); %#ok<SAGROW>
        confEstOrder(end+1, 1) = estO; %#ok<SAGROW>
    end
end
ConfT = local_build_table(confRows, [1 2 3 4 5], 6:19, ...
    {'Hypothesis','ROI','Band','Metric','Contrast','Estimate','SE','DF','tStat','pValue','pHolm_confirmatory', ...
     'CI_Low','CI_High','p_orderAdjusted','N_pairedSubjects','PairedMeanDiff','Cohens_dz','p_pairedT','p_LOO_worst'});
ConfT.Estimate_orderAdjusted = confEstOrder(:);
ConfT.p_Wilcoxon = confWilcoxon(:);
ConfT.LOO_worstSubject = confLOOsubj(:);
% Holm across the WHOLE confirmatory list (all hypotheses x metrics)
[holmP, holmSig] = local_holm(ConfT.pValue, ALPHA);
ConfT.pHolm_confirmatory = holmP;
ConfT.SignificantHolm = holmSig;

%% 7. SAVE
save(fullfile(OUT_DIR, 'roi_v2_statistics.mat'), 'Roi2T', 'OmnibusT', 'AllContrastsT', 'EmmT', 'ConfT', 'PARAMS11', 'issues');
save(fullfile(OUT_DIR, 'roi_v2_bandpower.mat'), 'Roi2T');
writetable(Roi2T, fullfile(OUT_DIR, 'roi_v2_bandpower.csv'));
writetable(ConfT, fullfile(OUT_DIR, 'roi_v2_confirmatory.csv'));
writetable(ConfT, fullfile(OUT_DIR, 'roi_v2_confirmatory.xlsx'));
writetable(OmnibusT, fullfile(OUT_DIR, 'roi_v2_omnibus.csv'));
writetable(OmnibusT, fullfile(OUT_DIR, 'roi_v2_omnibus.xlsx'));
writetable(AllContrastsT, fullfile(OUT_DIR, 'roi_v2_all_contrasts.csv'));
writetable(AllContrastsT, fullfile(OUT_DIR, 'roi_v2_all_contrasts.xlsx'));
writetable(EmmT, fullfile(OUT_DIR, 'roi_v2_emm.csv'));

%% 8. REPORT
L = {};
L{end+1} = sprintf('RUN_MODE = %s | subjects excluded in this run: %s', RUN_MODE_LABEL, local_v3_list(RUN_MODE_EXCL_SUBJ));
L{end+1} = repmat('=', 1, 70);
L{end+1} = 'PHASE 6B -- ROI v2 STATISTICS REPORT';
L{end+1} = repmat('=', 1, 70);
L{end+1} = sprintf('Generated: %s', datestr(now));
L{end+1} = sprintf('Recordings: %d | Subjects: %d | Condition pairs: %s', nRec, numel(uSubj), pairStr);
L{end+1} = 'Model: Metric ~ Condition + (1|Subject), REML, Satterthwaite df. Contrast = first minus second.';
L{end+1} = '';
for r = 1:nROI
    L{end+1} = sprintf('ROI %-15s: %s', ROI2_NAMES{r}, strjoin(ROI2_DEFINITION.(ROI2_NAMES{r}), ' ')); %#ok<SAGROW>
end
L{end+1} = '';
L{end+1} = ['PRESPECIFICATION STATUS: ' PARAMS11.PrespecificationStatus];
L{end+1} = '';
L{end+1} = '--- CONFIRMATORY PLANNED CONTRASTS (Holm across all rows below) ---';
for k = 1:height(ConfT)
    L{end+1} = sprintf('  %s %s %s [%s] %s: est=%.4f (95%% CI %.4f to %.4f) t(%.1f)=%.3f p=%.4f pHolm=%.4f %s', ...
        ConfT.Hypothesis{k}, ConfT.ROI{k}, ConfT.Band{k}, ConfT.Metric{k}, ConfT.Contrast{k}, ConfT.Estimate(k), ...
        ConfT.CI_Low(k), ConfT.CI_High(k), ConfT.DF(k), ConfT.tStat(k), ConfT.pValue(k), ConfT.pHolm_confirmatory(k), ...
        local_yn(ConfT.SignificantHolm(k))); %#ok<SAGROW>
    L{end+1} = sprintf('      robustness: order-adjusted est=%.4f p=%.4f | paired subset n=%d meandiff=%.4f dz=%.2f p_t=%.4f p_Wilcoxon=%.4f | worst leave-one-out p=%.4f (drop subject %g)', ...
        ConfT.Estimate_orderAdjusted(k), ConfT.p_orderAdjusted(k), ConfT.N_pairedSubjects(k), ConfT.PairedMeanDiff(k), ...
        ConfT.Cohens_dz(k), ConfT.p_pairedT(k), ConfT.p_Wilcoxon(k), ConfT.p_LOO_worst(k), ConfT.LOO_worstSubject(k)); %#ok<SAGROW>
end
L{end+1} = '';
L{end+1} = '--- SECONDARY: omnibus Condition test, BH-FDR within band x metric across the 4 ROIs ---';
for k = 1:height(OmnibusT)
    L{end+1} = sprintf('  [%s] %-15s %-25s F(%g,%.1f)=%.3f p=%.4f pFDR=%.4f %s', OmnibusT.Metric{k}, OmnibusT.ROI{k}, ...
        OmnibusT.Band{k}, OmnibusT.DF1(k), OmnibusT.DF2(k), OmnibusT.FStat(k), OmnibusT.pValue(k), OmnibusT.pFDR(k), ...
        local_yn(OmnibusT.SignificantFDR(k))); %#ok<SAGROW>
end
L{end+1} = '';
L{end+1} = '--- SECONDARY: all contrasts (UNGATED; Holm within model; gate = omnibus passed FDR) ---';
for k = 1:height(AllContrastsT)
    gateStr = 'gate=no';
    if AllContrastsT.OmnibusPassedFDR(k) == 1; gateStr = 'gate=yes'; end
    L{end+1} = sprintf('  [%s] %-15s %-25s %-6s est=%.4f (95%% CI %.4f to %.4f) p=%.4f pHolm=%.4f %s', ...
        AllContrastsT.Metric{k}, AllContrastsT.ROI{k}, AllContrastsT.Band{k}, AllContrastsT.Contrast{k}, ...
        AllContrastsT.Estimate(k), AllContrastsT.CI_Low(k), AllContrastsT.CI_High(k), AllContrastsT.pValue(k), ...
        AllContrastsT.pHolm_withinModel(k), gateStr); %#ok<SAGROW>
end
L{end+1} = '';
if isempty(issues)
    L{end+1} = 'Integrity issues: none.';
else
    L{end+1} = '--- INTEGRITY ISSUES ---';
    for k = 1:numel(issues); L{end+1} = ['  - ' issues{k}]; end %#ok<SAGROW>
end
L{end+1} = '';
L{end+1} = 'INTERPRETATION NOTES';
L{end+1} = '  - Only the CONFIRMATORY block carries a confirmatory claim. Everything else is secondary.';
L{end+1} = '  - R is eyes-closed while E and N are eyes-open: R vs E/N differences mostly reflect eye state';
L{end+1} = '    (and possibly drowsiness for theta), not the content of the intervention. E - N is the';
L{end+1} = '    contrast that addresses the Reels-vs-control question.';
L{end+1} = repmat('=', 1, 70);

fid = fopen(fullfile(REPORT_DIR, sprintf('phase6b_roi_v2_statistics_report_%s.txt', RUN_MODE_LABEL)), 'w');
fprintf(fid, '%s\n', L{:});
fclose(fid);
fprintf('%s\n', L{:});
fprintf('\n=== X25_roi_v2_statistics.m complete. ===\n');

%% ======================= LOCAL FUNCTIONS =======================

function subT = local_prepare(Roi2T, roiName, bandName, condOrder)
% Rows for one ROI x band, with Subject/Condition as categoricals.
    subT = Roi2T(strcmp(Roi2T.ROI, roiName) & strcmp(Roi2T.Band, bandName), :);
    subT.Subject = categorical(local_subject_numeric(subT.Subject));
    subT.Condition = categorical(cellstr(subT.Condition), condOrder);
end

function emm = local_emm(lme)
% Estimated marginal means for E, N, R under reference coding (ref = E).
    cn = lme.CoefficientNames; beta = lme.Coefficients.Estimate;
    L = zeros(3, numel(cn));
    L(:, strcmp(cn, '(Intercept)')) = 1;
    L(2, strcmp(cn, 'Condition_N')) = 1;
    L(3, strcmp(cn, 'Condition_R')) = 1;
    emm = L * beta;
end

function res = local_contrasts(lme, contrastRows, refLevel)
% res(k,:) = [est se df t p ciLow ciHigh] for each {lvl2, lvl1, label} row.
    cn = lme.CoefficientNames; beta = lme.Coefficients.Estimate; covB = lme.CoefficientCovariance;
    res = nan(size(contrastRows, 1), 7);
    for k = 1:size(contrastRows, 1)
        C = local_condition_contrast(cn, 'Condition', refLevel, contrastRows{k, 2}, contrastRows{k, 1});
        est = C * beta; se = sqrt(C * covB * C');
        [p, ~, ~, df2] = coefTest(lme, C, 0, 'DFMethod', 'satterthwaite');
        tc = tinv(0.975, df2);
        res(k, :) = [est, se, df2, est / se, p, est - tc * se, est + tc * se];
    end
end

function [x2, x1] = local_paired(subT, metricName, lvl2, lvl1)
% Values of lvl2 and lvl1 for subjects who have BOTH conditions.
    sNum = local_subject_numeric(cellstr(subT.Subject));
    cond = cellstr(subT.Condition);
    uS = unique(sNum); x2 = []; x1 = [];
    for s = 1:numel(uS)
        a = subT.(metricName)(sNum == uS(s) & strcmp(cond, lvl2));
        b = subT.(metricName)(sNum == uS(s) & strcmp(cond, lvl1));
        if numel(a) == 1 && numel(b) == 1 && isfinite(a) && isfinite(b)
            x2(end+1, 1) = a; x1(end+1, 1) = b; %#ok<AGROW>
        end
    end
end

function n = local_subject_numeric(subjVal)
    if iscategorical(subjVal); subjVal = cellstr(subjVal); end
    if iscell(subjVal)
        n = cellfun(@(s) str2double(regexprep(char(s), '\D', '')), subjVal);
    elseif isnumeric(subjVal)
        n = double(subjVal);
    else
        n = str2double(regexprep(char(subjVal), '\D', ''));
    end
end

function [adjP, reject] = local_fdr_bh(pvals, alphaFDR)
    pvals = pvals(:); m = numel(pvals);
    [sp, si] = sort(pvals); adj = nan(m, 1); runMin = 1;
    for i = m:-1:1
        runMin = min(runMin, sp(i) * m / i); adj(i) = runMin;
    end
    adjP = nan(m, 1); adjP(si) = adj; adjP = min(adjP, 1);
    reject = adjP <= alphaFDR;
end

function [adjP, reject] = local_holm(pvals, alphaFWER)
    pvals = pvals(:); m = numel(pvals);
    [sp, si] = sort(pvals); adj = nan(m, 1); runMax = 0;
    for i = 1:m
        runMax = max(runMax, sp(i) * (m - i + 1)); adj(i) = runMax;
    end
    adjP = nan(m, 1); adjP(si) = adj; adjP = min(adjP, 1);
    reject = adjP <= alphaFWER;
end

function cRow = local_condition_contrast(coeffNames, varName, refLevel, level1, level2)
% Row vector for mean(level2) - mean(level1) under reference coding.
    cRow = zeros(1, numel(coeffNames));
    if ~strcmp(level1, refLevel)
        idx1 = find(strcmp(coeffNames, sprintf('%s_%s', varName, level1)));
        if isempty(idx1); error('Phase6B:Contrast', 'Coefficient %s_%s not found.', varName, level1); end
        cRow(idx1) = cRow(idx1) - 1;
    end
    if ~strcmp(level2, refLevel)
        idx2 = find(strcmp(coeffNames, sprintf('%s_%s', varName, level2)));
        if isempty(idx2); error('Phase6B:Contrast', 'Coefficient %s_%s not found.', varName, level2); end
        cRow(idx2) = cRow(idx2) + 1;
    end
end

function s = local_yn(tf)
    if tf; s = 'SIGNIFICANT'; else; s = 'n.s.'; end
end

function T = local_build_table(rows, textColIdx, numColIdx, varNames)
% table() built column by column so text columns stay cellstr.
    nCols = numel(varNames); cols = cell(1, nCols);
    for j = 1:nCols
        if ismember(j, textColIdx)
            cols{j} = rows(:, j);
        elseif ismember(j, numColIdx)
            if isempty(rows); cols{j} = zeros(0, 1); else; cols{j} = cell2mat(rows(:, j)); end
        else
            error('Phase6B:BuildTable', 'Column %d not classified.', j);
        end
    end
    T = table(cols{:}, 'VariableNames', varNames);
end

function s = local_v3_list(v)
    if isempty(v); s = 'none'; else; s = strjoin(cellstr(num2str(v(:))), ', '); end
end
