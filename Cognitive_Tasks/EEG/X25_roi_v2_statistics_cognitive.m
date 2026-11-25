%% ========================================================================
% X25_roi_v2_statistics_cognitive.m
%
% PORTED from the intervention pipeline's X25_roi_v2_statistics.m so the
% cognitive-task results have the same confirmatory-analysis layer as the
% intervention results. What changed in the port (everything else --
% the model, the FDR/Holm scheme, the three robustness checks, the
% confirmatory-vs-secondary separation -- is the intervention X25's own
% logic, copied unchanged):
%   - Runs as THREE SEPARATE, PARALLEL analyses, one per Task (Full /
%     Nback / SART), exactly like X05/X06/X08/X09/X11-X20. Task is NEVER
%     modeled as a factor and is NEVER pooled -- including for the Holm
%     correction across the confirmatory list (Section 6), which is
%     computed independently per Task, not once across all three.
%   - RUN_MODE toggle: this pipeline's own EXCLUSION_CONFIG_cognitive.m
%     (EXCLUDE_NOISY_SUBJECTS = [7,23] -- see that file for the
%     per-subject rationale),
%     applied the SAME way X08_roi_level_statistics_cognitive.m applies it
%     (inline, directly to the channel-level BandPowerT for this Task,
%     before ROI aggregation) rather than through the intervention
%     pipeline's separate apply_run_mode_filter.m helper, which this
%     pipeline does not have (its exclusion logic has always been inlined
%     per-script -- a code-organization difference, not a methodology
%     one).
%   - Input: X04_aggregate_cognitive_subject_level.m's channel-level
%     BandPowerT (Subject, TaskType, BlockNumber, Condition, Channel,
%     Band, AbsolutePower, RelativePower), filtered to this Task and to
%     Condition in {E,N,R} (Baseline excluded) -- identical filtering to
%     X08's own Section 1b. ROI v2 membership/hypotheses come from
%     X24_roi_v2_definitions_cognitive.m's single, shared (non-Task)
%     output.
%   - CONTRASTS / REF_LEVEL: this pipeline's own N-E / R-E / R-N
%     convention (X05/X06/X11-X20), REF_LEVEL='E' -- the OPPOSITE
%     sign/order of the intervention X25's E-N/E-R/N-R. The two
%     confirmatory hypotheses (H1: Occipital Alpha R-E; H2: FrontalMidline
%     Theta N-E) already use this convention's own labels -- see
%     X24_roi_v2_definitions_cognitive.m.
%   - ORDER-ADJUSTED ROBUSTNESS CHECK (Section 6a): uses BlockNumber
%     (1/2 -- this pipeline's within-subject session-order variable) as
%     the additional fixed effect, in place of the intervention X25's
%     InterventionNumber. Same purpose: check whether the confirmatory
%     contrast survives adding session order as a covariate.
%   - local_subject_numeric / local_fdr_bh / local_holm_bonferroni /
%     local_condition_contrast / local_build_table / local_yn: copied
%     verbatim (in convention) from X08_roi_level_statistics_cognitive.m
%     (this pipeline's own existing ROI-level script), NOT re-derived.
%
% MODEL (identical to X05/X06/X08/X09/X11-X20): Metric ~ Condition +
% (1|Subject), REML, Satterthwaite df. Condition reference = E. Contrast
% direction = this pipeline's own N-E/R-E/R-N convention (lvl2 - lvl1).
%
% ROI AGGREGATION (identical rule to X04/X08): ROI AbsolutePower = mean of
% the member channels' AbsolutePower; ROI LogAbsolutePower =
% log(ROI AbsolutePower); ROI RelativePower = mean of the member channels'
% RelativePower.
%
% WHAT IS DIFFERENT FROM X08 (v1 ROI analysis) -- identical in kind to
% what the intervention X25 says about its own v1 ROI script (X08 there)
%   1. ROIs: the small, literature-based ROI v2 set from X24 (4-6
%      electrodes each, no edge electrodes) instead of the 5-region full
%      partition.
%   2. Confirmatory tests are PLANNED CONTRASTS (H1, H2; see X24), tested
%      directly from the LMM -- NOT gated by an omnibus F test and NOT
%      diluted by an FDR family that includes ROIs without a hypothesis.
%      Multiplicity across the short confirmatory list (2 hypotheses x 2
%      metrics = 4 rows, PER TASK) is controlled with Holm (FWER),
%      independently per Task.
%   3. Every confirmatory result comes with three robustness checks:
%        (a) the same contrast with BlockNumber (session order) as an
%            extra fixed effect;
%        (b) a paired-subset check on the subjects who did BOTH
%            conditions of the contrast (paired t, Wilcoxon signed-rank,
%            Cohen's dz);
%        (c) leave-one-subject-out: the largest p-value obtained when any
%            single subject is dropped.
%   4. Secondary analysis: every ROI x band x metric gets an omnibus LMM
%      test (BH-FDR within each band x metric across the 4 ROIs) and all
%      three contrasts are reported UNGATED with Holm within the model,
%      together with a flag saying whether the omnibus passed FDR.
%
% REQUIREMENTS
%   - Statistics and Machine Learning Toolbox (fitlme, anova, coefTest,
%     ttest, signrank).
%   - X04 output: 04_subject_level/phase3_subject_level_bandpower.mat
%   - X24 output: 24_roi_v2_definitions/roi_v2_definitions.mat
%
% OUTPUT (one complete set PER TASK PER RUN_MODE)
%   <ANALYSIS_ROOT>/<Full|Nback|SART>/25_roi_v2_stats/<RUN_MODE>/
%       roi_v2_bandpower.mat / .csv
%       roi_v2_confirmatory.csv / .xlsx
%       roi_v2_omnibus.csv / .xlsx
%       roi_v2_all_contrasts.csv / .xlsx
%       roi_v2_emm.csv
%       roi_v2_statistics.mat
%   <ANALYSIS_ROOT>/09_reports/
%       phase6b_roi_v2_statistics_report_<Full|Nback|SART>_<RunModeLabel>.txt
% ========================================================================

clear; clc;

%% 0. SETTINGS, PATHS, EXCLUSION CONFIG
config;
DATA_ROOT_FOR_SIBLING = DATA_DIR;
[DATA_ROOT_PARENT, ~] = fileparts(DATA_ROOT_FOR_SIBLING);
ANALYSIS_ROOT = OUT_DIR;

SUBJ_DIR   = fullfile(ANALYSIS_ROOT, '04_subject_level');
ROI2_DIR   = fullfile(ANALYSIS_ROOT, '24_roi_v2_definitions');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');
bandPowerFile = fullfile(SUBJ_DIR, 'phase3_subject_level_bandpower.mat');
roi2File      = fullfile(ROI2_DIR, 'roi_v2_definitions.mat');

if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

% ---- Noisy-subject exclusion list (same mechanism as X05/X06/X08/X09/X11-X20) ----
run(fullfile(fileparts(mfilename('fullpath')), 'EXCLUSION_CONFIG_cognitive.m'));

% USER-EDITABLE: run this script ONCE PER VALUE below to get both the
% full-sample and the sensitivity (noisy-subjects-excluded) results, saved
% side by side under separate subfolders.
RUN_MODE = 'sensitivity_excl_noisy';  % or: 'sensitivity_excl_noisy'
if ~ismember(RUN_MODE, RUN_MODES_cognitive)
    error('Phase6B:BadRunMode', 'RUN_MODE must be one of: %s', strjoin(RUN_MODES_cognitive, ', '));
end
RUN_MODE_LABEL = RUN_MODES_DISPLAY_cognitive.(RUN_MODE);
fprintf('RUN_MODE = %s\n', RUN_MODE_LABEL);
if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
    fprintf('Excluding subject(s) {%s} -- see EXCLUSION_CONFIG_cognitive.m for the reasons.\n', ...
        strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','));
end
% ---- END NEW ----

ALPHA = 0.05;
CONDITION_ORDER = {'E', 'N', 'R'};
REF_LEVEL = 'E';
CONTRASTS = {'N','E','N - E'; 'R','E','R - E'; 'R','N','R - N'};   % {lvl2, lvl1, label}; label = lvl2 - lvl1 -- this pipeline's own convention (X05/X06/X11-X20)
METRICS = {'LogAbsolutePower', 'RelativePower'};

for f = {bandPowerFile, roi2File}
    if ~exist(f{1}, 'file')
        error('Phase6B:MissingInput', ['Required input not found: %s\nRun X04_aggregate_cognitive_subject_level.m ' ...
            'and X24_roi_v2_definitions_cognitive.m first.'], f{1});
    end
end
if isempty(ver('stats'))
    error('Phase6B:MissingToolbox', 'Statistics and Machine Learning Toolbox is required.');
end

fprintf('========================================================\n');
fprintf('PHASE 6B -- ROI v2 STATISTICS, COGNITIVE-TASK DATA\n');
fprintf('Running as THREE INDEPENDENT, PARALLEL analyses (Full / Nback / SART).\n');
fprintf('========================================================\n');

%% 1. LOAD (shared across the Task loop)
S4 = load(bandPowerFile, 'BandPowerT', 'PARAMS3');
BandPowerT_ALL = S4.BandPowerT; PARAMS3 = S4.PARAMS3; clear S4;

requiredCols = {'Subject','TaskType','BlockNumber','Condition','Channel','Band','AbsolutePower','RelativePower'};
missingCols = requiredCols(~ismember(requiredCols, BandPowerT_ALL.Properties.VariableNames));
if ~isempty(missingCols)
    error('Phase6B:MissingColumns', 'BandPowerT is missing expected column(s): %s', strjoin(missingCols, ', '));
end

S24 = load(roi2File, 'ROI2_DEFINITION', 'ROI2_NAMES', 'CONFIRMATORY', 'CONFIRMATORY_METRICS', 'PARAMS24');
ROI2_DEFINITION = S24.ROI2_DEFINITION; ROI2_NAMES = S24.ROI2_NAMES;
CONFIRMATORY = S24.CONFIRMATORY; CONF_METRICS = S24.CONFIRMATORY_METRICS; PARAMS24 = S24.PARAMS24; clear S24;

ALL_BANDS = fieldnames(PARAMS3.BandDefinitionsHz);
nROI = numel(ROI2_NAMES); nBands = numel(ALL_BANDS); nMetrics = numel(METRICS);

TASK_TYPES = {'full', 'nback', 'sart'};
TASK_DISPLAY = struct('full', 'Full', 'nback', 'Nback', 'sart', 'SART');

%% ========================================================================
%  TASK LOOP -- everything below is repeated independently per Task
% ========================================================================
for taskI = 1:numel(TASK_TYPES)
    TASK_KEY = TASK_TYPES{taskI};
    TASK_LABEL = TASK_DISPLAY.(TASK_KEY);

    OUT_DIR = fullfile(ANALYSIS_ROOT, TASK_LABEL, '25_roi_v2_stats', RUN_MODE);
    if ~exist(OUT_DIR, 'dir'); mkdir(OUT_DIR); end

    fprintf('\n########################################################\n');
    fprintf('### TASK = %s | RUN_MODE = %s\n', TASK_LABEL, RUN_MODE_LABEL);
    fprintf('########################################################\n\n');

    %% Section 2: filter to this Task, then Condition in {E,N,R}, then RUN_MODE exclusion
    thisTaskT = BandPowerT_ALL(strcmp(BandPowerT_ALL.TaskType, TASK_KEY), :);
    nBaselineExcludedRows = sum(strcmp(thisTaskT.Condition, 'Baseline'));
    BandPowerT = thisTaskT(ismember(thisTaskT.Condition, CONDITION_ORDER), :);
    fprintf('Task=%s: %d total channel-level rows loaded, %d Baseline rows excluded, %d E/N/R rows retained.\n', ...
        TASK_LABEL, height(thisTaskT), nBaselineExcludedRows, height(BandPowerT));

    nExcludedRowsSensitivity = 0;
    if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
        subjNumForExcl = local_subject_numeric(BandPowerT.Subject);
        excludeMaskSensitivity = ismember(subjNumForExcl, EXCLUDE_NOISY_SUBJECTS);
        nExcludedRowsSensitivity = sum(excludeMaskSensitivity);
        BandPowerT = BandPowerT(~excludeMaskSensitivity, :);
        fprintf('Task=%s: RUN_MODE=%s -- excluded %d channel-level row(s) belonging to subject(s) {%s}; %d rows remain.\n', ...
            TASK_LABEL, RUN_MODE_LABEL, nExcludedRowsSensitivity, strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','), height(BandPowerT));
    end

    if isempty(BandPowerT)
        warning('Task=%s: no E/N/R rows found after filtering -- skipping this Task entirely.', TASK_LABEL);
        continue;
    end
    BandPowerT.ChannelUC = upper(BandPowerT.Channel);
    issues = {};

    %% Section 3: PRE-AGGREGATION VALIDATION (structural only -- see X08 header item 4)
    fails = 0;
    comboT = unique(BandPowerT(:, {'Subject','BlockNumber','Condition'}), 'rows');
    nRec = height(comboT);
    subjNum = local_subject_numeric(comboT.Subject);
    uSubj = unique(subjNum);
    pairCount = containers.Map('KeyType', 'char', 'ValueType', 'double');
    condsPerSubjCounts = nan(numel(uSubj), 1);
    for i = 1:numel(uSubj)
        cs = sort(cellstr(comboT.Condition(subjNum == uSubj(i))));
        nDistinct = numel(unique(cs));
        condsPerSubjCounts(i) = nDistinct;
        if numel(cs) ~= 2
            % Reported, not hard-gated -- see X08 header item 4: the exact
            % per-subject recording count is a property of this dataset,
            % not assumed in advance.
            issues{end+1} = sprintf('Subject %g has %d recordings (expected 2) for Task=%s.', uSubj(i), numel(cs), TASK_LABEL); %#ok<AGROW>
            continue;
        end
        if strcmp(cs{1}, cs{2})
            fprintf('  *** FAIL: subject %g has the SAME condition (%s) in both recordings.\n', uSubj(i), cs{1}); fails = fails + 1;
            continue;
        end
        key = strjoin(cs', '&');
        if isKey(pairCount, key); pairCount(key) = pairCount(key) + 1; else; pairCount(key) = 1; end
    end
    % Structural design-assumption check (identical logic/purpose to X08's
    % own check): if every subject in this Task's E/N/R data turns out to
    % have all 3 conditions, the incomplete-block-design justification for
    % an LMM instead of RM-ANOVA no longer applies.
    if ~isempty(condsPerSubjCounts) && all(condsPerSubjCounts >= 3)
        issues{end+1} = sprintf(['Task=%s: every subject appears to have 3+ conditions -- the incomplete-design ' ...
            'justification for using an LMM instead of RM-ANOVA no longer applies. STOP and re-examine.'], TASK_LABEL); %#ok<AGROW>
        fprintf('  *** FAIL: every subject has 3+ distinct conditions -- incomplete-design justification no longer applies.\n');
        fails = fails + 1;
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
        error('Phase6B:ValidationFailed', 'Task=%s: %d validation failure(s). Nothing fitted.', TASK_LABEL, fails);
    end

    %% Section 4: ROI AGGREGATION
    roiRows = cell(0, 7);
    for i = 1:nRec
        recMask = strcmp(BandPowerT.Subject, comboT.Subject{i}) & BandPowerT.BlockNumber == comboT.BlockNumber(i) & strcmp(BandPowerT.Condition, comboT.Condition{i});
        for r = 1:nROI
            chans = upper(ROI2_DEFINITION.(ROI2_NAMES{r}));
            chanMask = ismember(BandPowerT.ChannelUC, chans);
            for b = 1:nBands
                rows = recMask & chanMask & strcmp(BandPowerT.Band, ALL_BANDS{b});
                if sum(rows) ~= numel(chans)
                    issues{end+1} = sprintf('Sub%s Block%g ROI %s %s: %d/%d channels found.', comboT.Subject{i}, ...
                        comboT.BlockNumber(i), ROI2_NAMES{r}, ALL_BANDS{b}, sum(rows), numel(chans)); %#ok<SAGROW>
                end
                roiRows(end+1, :) = {comboT.Subject{i}, comboT.BlockNumber(i), comboT.Condition{i}, ROI2_NAMES{r}, ...
                    ALL_BANDS{b}, mean(BandPowerT.AbsolutePower(rows)), mean(BandPowerT.RelativePower(rows))}; %#ok<SAGROW>
            end
        end
    end
    Roi2T = local_build_table(roiRows, [1 3 4 5], [2 6 7], {'Subject','BlockNumber','Condition','ROI','Band','AbsolutePower','RelativePower'});
    Roi2T.LogAbsolutePower = log(Roi2T.AbsolutePower);
    Roi2T.LogAbsolutePower(Roi2T.AbsolutePower <= 0) = NaN;
    fprintf('  ROI aggregation complete (Task=%s): %d rows (%d recordings x %d ROIs x %d bands).\n', TASK_LABEL, height(Roi2T), nRec, nROI, nBands);

    %% Section 5: FIT ALL MODELS (ROI x band x metric)
    lmeStore = cell(nROI, nBands, nMetrics);
    omniRows = cell(0, 7); emmRows = cell(0, 6); allRows = cell(0, 13);
    for r = 1:nROI
        for b = 1:nBands
            for m = 1:nMetrics
                subT = local_prepare(Roi2T, ROI2_NAMES{r}, ALL_BANDS{b}, CONDITION_ORDER);
                try
                    lme = fitlme(subT, sprintf('%s ~ Condition + (1|Subject)', METRICS{m}), 'FitMethod', 'REML');
                catch ME
                    issues{end+1} = sprintf('fitlme failed: %s %s %s (Task=%s): %s', ROI2_NAMES{r}, ALL_BANDS{b}, METRICS{m}, TASK_LABEL, ME.message); %#ok<SAGROW>
                    continue;
                end
                lmeStore{r, b, m} = lme;
                aov = anova(lme, 'DFMethod', 'satterthwaite');
                ri = find(strcmp(aov.Term, 'Condition'), 1);
                omniRows(end+1, :) = {ROI2_NAMES{r}, METRICS{m}, ALL_BANDS{b}, aov.FStat(ri), aov.DF1(ri), aov.DF2(ri), aov.pValue(ri)}; %#ok<SAGROW>
                emm = local_emm(lme);
                emmRows(end+1, :) = {ROI2_NAMES{r}, METRICS{m}, ALL_BANDS{b}, emm(1), emm(2), emm(3)}; %#ok<SAGROW>
                cres = local_contrasts(lme, CONTRASTS, REF_LEVEL);
                adj = local_holm_bonferroni(cres(:, 5), ALPHA);
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

    %% Section 6: SECONDARY -- BH-FDR of the omnibus test within each band x metric (across ROIs)
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

    %% Section 7: CONFIRMATORY PLANNED CONTRASTS (+ robustness checks)
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
                issues{end+1} = sprintf('Confirmatory %s %s (Task=%s): model missing.', CONFIRMATORY.Hypothesis{i}, metricName, TASK_LABEL); %#ok<SAGROW>
                continue;
            end
            cres = local_contrasts(lme, CONTRASTS(ci, :), REF_LEVEL);   % est se df t p lo hi

            % (a) order-adjusted (BlockNumber as extra fixed effect)
            subT = local_prepare(Roi2T, roiName, bandName, CONDITION_ORDER);
            subT.BlockNumber = categorical(subT.BlockNumber);
            try
                lmeO = fitlme(subT, sprintf('%s ~ Condition + BlockNumber + (1|Subject)', metricName), 'FitMethod', 'REML');
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
    % Holm across THIS TASK's own confirmatory list only (never pooled across Tasks)
    [holmP, holmSig] = local_holm_bonferroni(ConfT.pValue, ALPHA);
    ConfT.pHolm_confirmatory = holmP;
    ConfT.SignificantHolm = holmSig;

    %% Section 8: SAVE
    save(fullfile(OUT_DIR, 'roi_v2_statistics.mat'), 'Roi2T', 'OmnibusT', 'AllContrastsT', 'EmmT', 'ConfT', 'PARAMS24', 'issues');
    save(fullfile(OUT_DIR, 'roi_v2_bandpower.mat'), 'Roi2T');
    writetable(Roi2T, fullfile(OUT_DIR, 'roi_v2_bandpower.csv'));
    writetable(ConfT, fullfile(OUT_DIR, 'roi_v2_confirmatory.csv'));
    writetable(ConfT, fullfile(OUT_DIR, 'roi_v2_confirmatory.xlsx'));
    writetable(OmnibusT, fullfile(OUT_DIR, 'roi_v2_omnibus.csv'));
    writetable(OmnibusT, fullfile(OUT_DIR, 'roi_v2_omnibus.xlsx'));
    writetable(AllContrastsT, fullfile(OUT_DIR, 'roi_v2_all_contrasts.csv'));
    writetable(AllContrastsT, fullfile(OUT_DIR, 'roi_v2_all_contrasts.xlsx'));
    writetable(EmmT, fullfile(OUT_DIR, 'roi_v2_emm.csv'));

    %% Section 9: REPORT
    L = {};
    L{end+1} = sprintf('Task=%s | RUN_MODE = %s | subjects excluded in this run: %s', TASK_LABEL, RUN_MODE_LABEL, local_list(EXCLUDE_NOISY_SUBJECTS, RUN_MODE));
    L{end+1} = repmat('=', 1, 70);
    L{end+1} = sprintf('PHASE 6B -- ROI v2 STATISTICS REPORT (Task=%s)', TASK_LABEL);
    L{end+1} = repmat('=', 1, 70);
    L{end+1} = sprintf('Generated: %s', datestr(now));
    L{end+1} = sprintf('Recordings: %d | Subjects: %d | Condition pairs: %s', nRec, numel(uSubj), pairStr);
    L{end+1} = 'Model: Metric ~ Condition + (1|Subject), REML, Satterthwaite df. Contrast = this pipeline''s own N-E/R-E/R-N convention (lvl2 - lvl1).';
    L{end+1} = '';
    for r = 1:nROI
        L{end+1} = sprintf('ROI %-15s: %s', ROI2_NAMES{r}, strjoin(ROI2_DEFINITION.(ROI2_NAMES{r}), ' ')); %#ok<SAGROW>
    end
    L{end+1} = '';
    L{end+1} = ['PRESPECIFICATION STATUS: ' PARAMS24.PrespecificationStatus];
    L{end+1} = '';
    L{end+1} = '--- CONFIRMATORY PLANNED CONTRASTS (Holm across all rows below, this Task only) ---';
    for k = 1:height(ConfT)
        L{end+1} = sprintf('  %s %s %s [%s] %s: est=%.4f (95%% CI %.4f to %.4f) t(%.1f)=%.3f p=%.4f pHolm=%.4f %s', ...
            ConfT.Hypothesis{k}, ConfT.ROI{k}, ConfT.Band{k}, ConfT.Metric{k}, ConfT.Contrast{k}, ConfT.Estimate(k), ...
            ConfT.CI_Low(k), ConfT.CI_High(k), ConfT.DF(k), ConfT.tStat(k), ConfT.pValue(k), ConfT.pHolm_confirmatory(k), ...
            local_yn(ConfT.SignificantHolm(k))); %#ok<SAGROW>
        L{end+1} = sprintf('      robustness: order-adjusted (BlockNumber) est=%.4f p=%.4f | paired subset n=%d meandiff=%.4f dz=%.2f p_t=%.4f p_Wilcoxon=%.4f | worst leave-one-out p=%.4f (drop subject %g)', ...
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
    L{end+1} = '  - This Task''s confirmatory result is independent of the other two Tasks'' -- never pooled, never compared directly here.';
    L{end+1} = repmat('=', 1, 70);

    fid = fopen(fullfile(REPORT_DIR, sprintf('phase6b_roi_v2_statistics_report_%s_%s.txt', TASK_LABEL, RUN_MODE_LABEL)), 'w');
    fprintf(fid, '%s\n', L{:});
    fclose(fid);
    fprintf('%s\n', L{:});
    fprintf('\nSaved output to: %s\n', OUT_DIR);

end % TASK_TYPES loop

fprintf('\n=== X25_roi_v2_statistics_cognitive.m complete for all three Tasks. ===\n');

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
% Identical implementation to X05/X06/X08/X09/X11-X20's own copy.
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
% Identical implementation to X05/X06/X08/X09/X11-X20's own copy.
    pvals = pvals(:); m = numel(pvals);
    [sp, si] = sort(pvals); adj = nan(m, 1); runMin = 1;
    for i = m:-1:1
        runMin = min(runMin, sp(i) * m / i); adj(i) = runMin;
    end
    adjP = nan(m, 1); adjP(si) = adj; adjP = min(adjP, 1);
    reject = adjP <= alphaFDR;
end

function [adjP, reject] = local_holm_bonferroni(pvals, alphaFWER)
% Identical implementation to X05/X06/X08/X09/X11-X20's own copy.
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
% Identical implementation to X05/X06/X08/X09/X11-X20's own copy.
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
% table() built column by column so text columns stay cellstr. Identical
% implementation to X05/X06/X08/X09/X11-X20's own copy.
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

function s = local_list(v, runMode)
% Printable excluded-subject list for the report header -- 'none' for
% full_sample regardless of the vector's contents (sensitivity-only list).
    if ~strcmp(runMode, 'sensitivity_excl_noisy') || isempty(v)
        s = 'none';
    else
        s = strjoin(cellstr(num2str(v(:))), ', ');
    end
end
