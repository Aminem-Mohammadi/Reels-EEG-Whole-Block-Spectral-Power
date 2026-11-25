%% ========================================================================
% X18B_electrode_all_contrasts_ungated_cognitive.m
%
% COMPANION SCRIPT TO X18_position_pos1_pos2_cognitive.m -- NOT a
% modification of X18, NOT a new numbered phase. A purely additive
% companion, same pattern and reason as X11B/X12B/X17B: closing a
% schema gap needed to build a full-scalp (ALL 62 electrodes) topographic
% figure for X18's own recording-position (block-order) family.
%
% WHY THIS SCRIPT EXISTS
%   X18's own phase18_pairwise_*_marginal_electrode.csv files only contain
%   rows for electrodes whose own term-family FDR already survived
%   BH-FDR. Not enough for a full-scalp topomap.
%   UNLIKE X11/X12/X17, X18 does NOT save its own per-electrode,
%   Position-labeled source table to disk anywhere -- verified by reading
%   X18_position_pos1_pos2_cognitive.m's every save() call directly: it
%   saves phase18_omnibus_electrode.mat (OmnibusElectrodeT -- F-stats
%   only, ungated, all electrodes) and the (gated) pairwise tables, but
%   never the labeled BandPowerT-with-Position table itself. That table is
%   built in-memory each run from TWO read-only sources X18 itself already
%   reads: X04's raw electrode BandPowerT, and X01's own audit table
%   (audit_table_cognitive.csv), joined by Subject+BlockNumber to derive
%   Position (First/Second) from TaskOrderSequenceRaw, with any
%   missing/ambiguous/unparseable match EXCLUDED (never guessed).
%   X18's situation needs more work than X11B/X12B/X17B's "just refit from
%   an already-saved table", since this Position-derivation + audit-join
%   logic has to be re-run -- but, exactly like X05B's own precedent for
%   X05's un-saved filtered input, it is STILL possible without running
%   X18 or any earlier script: every input this script needs (X04's raw
%   output, X01's audit table) is already on disk, read-only.
%   This script changes NOTHING about X18's own models, FDR, or Holm
%   logic, and does not touch X05/X05B/X09/X09B/X11/X11B/X12/X12B/X17/X17B
%   or their companions. Electrode level ONLY. It does NOT run X18 itself.
%
% WHAT THIS SCRIPT DOES
%   For each RUN_MODE in RUN_MODES_cognitive x each Task in {Nback, SART}
%   (4 independent passes, never pooled -- Full is categorically out of
%   scope, identical restriction to X18 itself: TASK_TYPES =
%   {'nback','sart'}):
%     1. Re-derives X18's own Task+RUN_MODE-filtered, Position-labeled
%        electrode BandPowerT, step for step, from X04's raw electrode
%        BandPowerT + X01's audit table -- IDENTICAL logic to X18's own
%        Sections 1-3b, reproduced verbatim below (filter TaskType,
%        exclude Baseline, join audit table on Subject+BlockNumber within
%        this TaskType, derive Position via local_position_from_seq,
%        exclude any un-Position-labeled row, then (RUN_MODE=
%        'sensitivity_excl_noisy' only) exclude EXCLUDE_NOISY_SUBJECTS).
%     2. Loads that Task+RUN_MODE's ALREADY-SAVED
%        phase18_omnibus_electrode.mat (OmnibusElectrodeT -- X18's own
%        ungated, all-electrode F-test table), for cross-checking.
%        REQUIRES X18 to already have been run for that Task+RUN_MODE;
%        errors out if missing.
%     3. Re-fits the SAME Electrode x Band x Metric Condition*Position
%        models X18 already fit (formula "<Metric> ~ Condition*Position +
%        (1|Subject)", REML, Condition ref='E', Position ref='First',
%        identical to X18's own local_fit_condition_position_omnibus,
%        reproduced verbatim below), for ALL 62 electrodes x 9 bands x 2
%        metrics (1116 models).
%     4. CROSS-CHECKS every successful refit's THREE omnibus F/DF1/DF2/p
%        (Condition, Position, Condition:Position) against X18's own
%        already-saved OmnibusElectrodeT rows for that exact
%        Electrode/Metric/Band/Term, within floating-point tolerance. Any
%        mismatch beyond tolerance is a hard STOP -- same safeguard as
%        X09B/X05B/X11B/X12B/X17B.
%     5. From each successful refit, computes ALL THREE pairwise-contrast
%        families, ungated, for EVERY electrode, using X18's OWN
%        contrast-construction functions reproduced VERBATIM:
%          - ConditionMarginal: N-E/R-E/R-N, averaged over Position,
%            Holm-3 (via local_condition_position_contrast w1=w2=0.5)
%          - PositionMarginal: Second-First, averaged over Condition,
%            single comparison, no correction
%            (local_position_marginal_contrast)
%          - InteractionSimple: N-E/R-E/R-N evaluated separately AT
%            Position=First and Position=Second (6 contrasts), Holm-6
%            (via local_condition_position_contrast w1=w2=0 or 1)
%     6. Saves ElectrodeAllContrastsT_X18:
%          - appended into the SAME phase18_omnibus_electrode.mat X18
%            already wrote (via save(...,'-append'));
%          - electrode_all_contrasts_ungated.csv, written to the SAME
%            18_position_pos1_pos2/<RUN_MODE> folder X18 already writes to.
%
% WHAT THIS SCRIPT DOES NOT DO
%   - Does not change X18's FDR, Holm, or any inferential conclusion.
%   - Does not touch Global-level output.
%   - Does not run X18 itself, and does not attempt Full (same categorical
%     exclusion as X18 itself -- logged, not silently skipped).
%   - Does not change how Position is derived or how ambiguous/missing
%     audit rows are handled -- reproduces X18's own exclusion logic
%     exactly, including its own integrity-issue counting.
%
% REQUIREMENTS
%   - MATLAB with Statistics and Machine Learning Toolbox.
%   - X18_position_pos1_pos2_cognitive.m must already have been run for
%     EVERY (Task, RUN_MODE) combination this script attempts.
%
% INPUT (read-only, not modified):
%   <ANALYSIS_ROOT>/04_subject_level/phase3_subject_level_bandpower.mat
%     (BandPowerT, PARAMS3 -- shared, Task-invariant source)
%   <ANALYSIS_ROOT>/01_audit/audit_table_cognitive.csv (shared, read-only)
%   <ANALYSIS_ROOT>/<Nback|SART>/18_position_pos1_pos2/<RUN_MODE>/
%     phase18_omnibus_electrode.mat (OmnibusElectrodeT -- produced by X18,
%     read-only)
%
% OUTPUT (added, per Task per RUN_MODE; nothing pre-existing is overwritten):
%   <ANALYSIS_ROOT>/<Nback|SART>/18_position_pos1_pos2/<RUN_MODE>/
%       electrode_all_contrasts_ungated.csv       (new)
%       phase18_omnibus_electrode.mat              (ElectrodeAllContrastsT_X18 appended)
%   <ANALYSIS_ROOT>/09_reports/
%       phase18B_electrode_all_contrasts_report_<Nback|SART>_<RunModeLabel>.txt
%
% The report, per Task per RUN_MODE, includes the cross-check section
% confirming the refit omnibus results match X18's own saved
% OmnibusElectrodeT exactly, and the Position-derivation counts
% (missing-audit-match / ambiguous / unparseable), which should match
% X18's own report for that Task/RUN_MODE.
% ========================================================================

clear; clc;

%% ------------------------------------------------------------------
%  0. SETTINGS, PATHS, EXISTENCE CHECKS
% ------------------------------------------------------------------

config;
DATA_ROOT_FOR_SIBLING = DATA_DIR;
[DATA_ROOT_PARENT, ~] = fileparts(DATA_ROOT_FOR_SIBLING);
ANALYSIS_ROOT = OUT_DIR;

SUBJ_DIR   = fullfile(ANALYSIS_ROOT, '04_subject_level');
AUDIT_DIR  = fullfile(ANALYSIS_ROOT, '01_audit');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');

bandPowerFile = fullfile(SUBJ_DIR, 'phase3_subject_level_bandpower.mat');
auditFile     = fullfile(AUDIT_DIR, 'audit_table_cognitive.csv');

if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

run(fullfile(fileparts(mfilename('fullpath')), 'EXCLUSION_CONFIG_cognitive.m'));

fprintf('========================================================\n');
fprintf('X18B -- ELECTRODE-LEVEL ALL-CONTRASTS, UNGATED (every electrode), COGNITIVE-TASK DATA\n');
fprintf('Companion to X18. Re-derives X18''s own filtered, Position-labeled electrode input\n');
fprintf('(from X04''s raw output + X01''s audit table, identical logic to X18''s own Sections 1-3b),\n');
fprintf('re-fits the SAME electrode-level Condition*Position models X18 already fit, per Task\n');
fprintf('per RUN_MODE, purely to export all three contrast families (ConditionMarginal/\n');
fprintf('PositionMarginal/InteractionSimple) for EVERY electrode (ungated, descriptive).\n');
fprintf('========================================================\n');
fprintf('ANALYSIS_ROOT   : %s\n', ANALYSIS_ROOT);
fprintf('--------------------------------------------------------\n\n');

if ~exist(bandPowerFile, 'file')
    error('Phase18B:MissingInput', 'X04 band-power file not found: %s', bandPowerFile);
end
if ~exist(auditFile, 'file')
    error('Phase18B:MissingInput', 'X01 audit table not found: %s', auditFile);
end
if isempty(ver('stats'))
    error('Phase18B:MissingToolbox', ...
        'Statistics and Machine Learning Toolbox not detected. This script requires fitlme, anova(LinearMixedModel), and coefTest.');
end

%% ------------------------------------------------------------------
%  1. SHARED (Task- and RUN_MODE-invariant) REFERENCE DATA
% ------------------------------------------------------------------

S4 = load(bandPowerFile, 'BandPowerT', 'PARAMS3');
BandPowerT_ALL = S4.BandPowerT;
PARAMS3 = S4.PARAMS3;
clear S4;

requiredElecCols = {'Subject','TaskType','BlockNumber','Condition','Channel','Band', ...
    'LogAbsolutePower','RelativePower'};
missingE = requiredElecCols(~ismember(requiredElecCols, BandPowerT_ALL.Properties.VariableNames));
if ~isempty(missingE)
    error('Phase18B:MissingColumns', 'BandPowerT is missing expected column(s): %s', strjoin(missingE, ', '));
end

ALL_BANDS = fieldnames(PARAMS3.BandDefinitionsHz);
nBands = numel(ALL_BANDS);
METRICS = {'LogAbsolutePower', 'RelativePower'};
nMetrics = numel(METRICS);
CONDITION_ORDER = {'E', 'N', 'R'};
REF_LEVEL = CONDITION_ORDER{1};
POSITION_ORDER = {'First', 'Second'};
POSITION_REF = POSITION_ORDER{1};
POSITION_SECOND = POSITION_ORDER{2};
% Reproduced VERBATIM from X18's own CONTRASTS declaration (N-E, R-E,
% R-N) -- NOT relabeled. See disclosure above.
CONTRASTS = {'N','E','N - E'; 'R','E','R - E'; 'R','N','R - N'};

TASK_TYPES = {'nback', 'sart'};   % Full is OUT OF SCOPE -- identical restriction to X18 itself
TASK_DISPLAY = struct('nback', 'Nback', 'sart', 'SART');

FIT_TOL_F = 1e-4;     % relative tolerance for cross-checking refit FStat vs X18's saved FStat
FIT_TOL_P = 1e-6;     % absolute tolerance for cross-checking refit pValue vs X18's saved pValue

%% ------------------------------------------------------------------
%  2. LOAD X01 AUDIT TABLE (read-only). Reproduced VERBATIM from X18's
%     own Section 2.
% ------------------------------------------------------------------
AuditT = readtable(auditFile);
requiredAuditCols = {'Subject','TaskType','BlockNumber','TaskOrderSequenceRaw','LoadStatus','ParseOK'};
missingAudit = requiredAuditCols(~ismember(requiredAuditCols, AuditT.Properties.VariableNames));
if ~isempty(missingAudit)
    error('Phase18B:MissingColumns', ...
        'audit_table_cognitive.csv is missing expected column(s): %s.', strjoin(missingAudit, ', '));
end
AuditT.Subject               = local_robust_cellstr(AuditT.Subject);
AuditT.TaskType               = local_robust_cellstr(AuditT.TaskType);
AuditT.TaskOrderSequenceRaw   = local_robust_cellstr(AuditT.TaskOrderSequenceRaw);
AuditT.LoadStatus             = local_robust_cellstr(AuditT.LoadStatus);
AuditT.ParseOK                = local_robust_cellstr(AuditT.ParseOK);
AuditOKMask = strcmpi(AuditT.LoadStatus, 'OK') & (strcmpi(AuditT.ParseOK, 'true') | strcmpi(AuditT.ParseOK, '1'));
AuditTOK = AuditT(AuditOKMask, :);
fprintf('Loaded audit table: %d total rows, %d with LoadStatus=OK and ParseOK=true.\n\n', height(AuditT), height(AuditTOK));

grandCrossCheckFails = 0;

%% ========================================================================
%  MAIN LOOP -- RUN_MODE x TASK, each a fully independent pass
% ========================================================================
for rmI = 1:numel(RUN_MODES_cognitive)
    RUN_MODE = RUN_MODES_cognitive{rmI};
    RUN_MODE_LABEL = RUN_MODES_DISPLAY_cognitive.(RUN_MODE);

    for taskI = 1:numel(TASK_TYPES)
        TASK_KEY = TASK_TYPES{taskI};
        TASK_LABEL = TASK_DISPLAY.(TASK_KEY);

        STATS_DIR = fullfile(ANALYSIS_ROOT, TASK_LABEL, '18_position_pos1_pos2', RUN_MODE);
        omnibusFile = fullfile(STATS_DIR, 'phase18_omnibus_electrode.mat');

        fprintf('\n########################################################\n');
        fprintf('### TASK = %s | RUN_MODE = %s\n', TASK_LABEL, RUN_MODE_LABEL);
        fprintf('########################################################\n\n');

        if ~exist(omnibusFile, 'file')
            error('Phase18B:MissingX18Output', ...
                ['X18''s saved output not found for Task=%s, RUN_MODE=%s:\n  %s\n' ...
                 'Run X18_position_pos1_pos2_cognitive.m first with RUN_MODE=''%s'' ' ...
                 'before running this companion script for that Task/RUN_MODE.'], ...
                TASK_LABEL, RUN_MODE_LABEL, omnibusFile, RUN_MODE);
        end
        So = load(omnibusFile, 'OmnibusElectrodeT');
        if ~isfield(So, 'OmnibusElectrodeT')
            error('Phase18B:MissingVariables', '%s does not contain OmnibusElectrodeT as expected from X18.', omnibusFile);
        end
        OmnibusElectrodeT_X18 = So.OmnibusElectrodeT;
        clear So;

        %% --- Section 3a: filter X04 electrode table to this Task, then E/N/R
        %      (VERBATIM logic from X18's own Section 3a) ---
        thisTaskElec = BandPowerT_ALL(strcmp(BandPowerT_ALL.TaskType, TASK_KEY), :);
        nBaselineExcludedElec = sum(strcmp(thisTaskElec.Condition, 'Baseline'));
        BandPowerT = thisTaskElec(ismember(thisTaskElec.Condition, CONDITION_ORDER), :);
        fprintf('Task=%s: Electrode rows %d total, %d Baseline excluded, %d E/N/R retained.\n', ...
            TASK_LABEL, height(thisTaskElec), nBaselineExcludedElec, height(BandPowerT));

        if isempty(BandPowerT)
            warning('Task=%s: no E/N/R rows found after filtering -- skipping this Task entirely.', TASK_LABEL);
            continue;
        end

        %% --- Section 3b: derive Position from X01's audit table
        %      (VERBATIM logic from X18's own Section 3b) ---
        auditThisTask = AuditTOK(strcmp(AuditTOK.TaskType, TASK_KEY), :);
        recCombos = unique(BandPowerT(:, {'Subject','BlockNumber'}), 'rows');
        nRecCombos = height(recCombos);
        posLabelByCombo = cell(nRecCombos, 1);
        nAmbiguous = 0; nMissingAudit = 0; nUnparseable = 0;
        for i = 1:nRecCombos
            subjStr = recCombos.Subject{i};
            blockNum = recCombos.BlockNumber(i);
            matchMask = strcmp(auditThisTask.Subject, subjStr) & auditThisTask.BlockNumber == blockNum;
            nMatch = sum(matchMask);
            if nMatch == 0
                posLabelByCombo{i} = '';
                nMissingAudit = nMissingAudit + 1;
                continue;
            elseif nMatch > 1
                posLabelByCombo{i} = '';
                nAmbiguous = nAmbiguous + 1;
                continue;
            end
            seqRaw = auditThisTask.TaskOrderSequenceRaw{matchMask};
            seqU = upper(strtrim(seqRaw));
            posLabel = local_position_from_seq(TASK_KEY, seqU);
            if isempty(posLabel)
                nUnparseable = nUnparseable + 1;
            end
            posLabelByCombo{i} = posLabel;
        end
        fprintf('Position derivation (Task=%s): %d recordings, %d missing-audit-match, %d ambiguous-match, %d unparseable SEQ -- all EXCLUDED, never guessed.\n', ...
            TASK_LABEL, nRecCombos, nMissingAudit, nAmbiguous, nUnparseable);

        BandPowerT.Position = repmat({''}, height(BandPowerT), 1);
        for i = 1:nRecCombos
            subjStr = recCombos.Subject{i};
            blockNum = recCombos.BlockNumber(i);
            eMask = strcmp(BandPowerT.Subject, subjStr) & BandPowerT.BlockNumber == blockNum;
            BandPowerT.Position(eMask) = {posLabelByCombo{i}};
        end
        nExcludedElec = sum(strcmp(BandPowerT.Position, ''));
        BandPowerT = BandPowerT(~strcmp(BandPowerT.Position, ''), :);
        fprintf('After excluding un-Position-labeled rows: Electrode %d rows retained (%d excluded).\n', ...
            height(BandPowerT), nExcludedElec);

        if isempty(BandPowerT)
            warning('Task=%s, RUN_MODE=%s: no rows remain after Position derivation -- skipping.', TASK_LABEL, RUN_MODE_LABEL);
            continue;
        end

        %% --- sensitivity exclusion (RUN_MODE='sensitivity_excl_noisy' only) ---
        nExcludedRowsSensitivityElec = 0;
        if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
            subjNumForExclElec = local_subject_numeric(BandPowerT.Subject);
            excludeMaskElec = ismember(subjNumForExclElec, EXCLUDE_NOISY_SUBJECTS);
            nExcludedRowsSensitivityElec = sum(excludeMaskElec);
            BandPowerT = BandPowerT(~excludeMaskElec, :);
            fprintf('Task=%s: RUN_MODE=%s -- excluded %d Electrode-level row(s) belonging to subject(s) {%s}; %d rows remain.\n', ...
                TASK_LABEL, RUN_MODE_LABEL, nExcludedRowsSensitivityElec, strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','), height(BandPowerT));
        end
        if isempty(BandPowerT)
            warning('Task=%s, RUN_MODE=%s: no rows remain after sensitivity exclusion -- skipping.', TASK_LABEL, RUN_MODE_LABEL);
            continue;
        end

        ELECTRODES = unique(BandPowerT.Channel);
        nElec = numel(ELECTRODES);

        %% ----------------------------------------------------------------
        %  REFIT all Electrode x Band x Metric Condition*Position models,
        %  cross-check against X18's own saved omnibus result (all 3
        %  terms), and extract all 3 contrast families (ungated).
        % ----------------------------------------------------------------
        allRows = cell(0, 13);
        crossCheckFails = {};
        nRefitOK = 0; nRefitFailed = 0; nCrossCheckOK = 0;
        TERM_KEYS = {'Condition', 'Position', 'Condition:Position'};

        for e = 1:nElec
            elecName = ELECTRODES{e};
            for m = 1:nMetrics
                metricName = METRICS{m};
                for b = 1:nBands
                    bandName = ALL_BANDS{b};
                    subT = BandPowerT(strcmp(BandPowerT.Channel, elecName) & strcmp(BandPowerT.Band, bandName), :);
                    fr = local_fit_condition_position_omnibus(subT, metricName, CONDITION_ORDER, POSITION_ORDER);
                    if ~fr.ok
                        nRefitFailed = nRefitFailed + 1;
                        continue;
                    end

                    mismatchHere = false;
                    for tI = 1:3
                        x18Row = OmnibusElectrodeT_X18(strcmp(OmnibusElectrodeT_X18.Electrode, elecName) & strcmp(OmnibusElectrodeT_X18.Metric, metricName) & ...
                            strcmp(OmnibusElectrodeT_X18.Band, bandName) & strcmp(OmnibusElectrodeT_X18.Term, TERM_KEYS{tI}), :);
                        if height(x18Row) ~= 1
                            mismatchHere = true;
                            crossCheckFails{end+1} = sprintf('%s/%s/%s/%s: expected exactly 1 matching row in X18''s OmnibusElectrodeT, found %d.', ...
                                elecName, metricName, bandName, TERM_KEYS{tI}, height(x18Row)); %#ok<AGROW>
                            continue;
                        end
                        fStat = fr.rows(tI).FStat; pVal = fr.rows(tI).pValue;
                        fMismatch = abs(fStat - x18Row.FStat(1)) > FIT_TOL_F * max(1, abs(x18Row.FStat(1)));
                        pMismatch = abs(pVal - x18Row.pValue(1)) > FIT_TOL_P;
                        if fMismatch || pMismatch
                            mismatchHere = true;
                            crossCheckFails{end+1} = sprintf('%s/%s/%s/%s: refit F=%.6g p=%.6g vs X18-saved F=%.6g p=%.6g (mismatch beyond tolerance).', ...
                                elecName, metricName, bandName, TERM_KEYS{tI}, fStat, pVal, x18Row.FStat(1), x18Row.pValue(1)); %#ok<AGROW>
                        end
                    end
                    if mismatchHere
                        continue;
                    end
                    nCrossCheckOK = nCrossCheckOK + 1;

                    cn = fr.lme.CoefficientNames;
                    beta = fr.lme.Coefficients.Estimate;
                    covB = fr.lme.CoefficientCovariance;

                    % --- ConditionMarginal: N-E/R-E/R-N averaged over Position, Holm-3 ---
                    cp3 = nan(3, 1); cdat3 = cell(3, 8);
                    for c = 1:3
                        lvl2 = CONTRASTS{c,1}; lvl1 = CONTRASTS{c,2}; label = CONTRASTS{c,3};
                        C = local_condition_position_contrast(cn, REF_LEVEL, lvl1, lvl2, POSITION_SECOND, 0.5, 0.5);
                        [est, se, tstat, p, df2c, ciLow, ciHigh] = local_apply_contrast(fr.lme, C, beta, covB);
                        cp3(c) = p;
                        cdat3(c, :) = {label, est, se, df2c, tstat, p, ciLow, ciHigh};
                    end
                    adj3 = local_holm_bonferroni(cp3, 0.05);
                    for c = 1:3
                        allRows(end+1, :) = [{elecName, metricName, bandName, 'ConditionMarginal'}, cdat3(c, 1:6), {adj3(c)}, cdat3(c, 7:8)]; %#ok<AGROW>
                    end

                    % --- PositionMarginal: Second-First averaged over Condition, no correction ---
                    Cpos = local_position_marginal_contrast(cn, CONDITION_ORDER, POSITION_SECOND);
                    [est, se, tstat, p, df2c, ciLow, ciHigh] = local_apply_contrast(fr.lme, Cpos, beta, covB);
                    allRows(end+1, :) = [{elecName, metricName, bandName, 'PositionMarginal'}, {'Second - First', est, se, df2c, tstat, p, p, ciLow, ciHigh}]; %#ok<AGROW>

                    % --- InteractionSimple: N-E/R-E/R-N at Position=First and Position=Second, Holm-6 ---
                    cp6 = nan(6, 1); cdat6 = cell(6, 9);
                    rowIdx6 = 0;
                    for p6 = 1:numel(POSITION_ORDER)
                        posLvl = POSITION_ORDER{p6};
                        w = 0; if strcmp(posLvl, POSITION_SECOND); w = 1; end
                        for c = 1:3
                            lvl2 = CONTRASTS{c,1}; lvl1 = CONTRASTS{c,2}; label = CONTRASTS{c,3};
                            rowIdx6 = rowIdx6 + 1;
                            C = local_condition_position_contrast(cn, REF_LEVEL, lvl1, lvl2, POSITION_SECOND, w, w);
                            [est2, se2, tstat2, p2, df2c2, ciLow2, ciHigh2] = local_apply_contrast(fr.lme, C, beta, covB);
                            cp6(rowIdx6) = p2;
                            cdat6(rowIdx6, :) = {sprintf('%s @ Position=%s', label, posLvl), est2, se2, df2c2, tstat2, p2, ciLow2, ciHigh2, NaN};
                        end
                    end
                    adj6 = local_holm_bonferroni(cp6, 0.05);
                    for r = 1:6
                        allRows(end+1, :) = [{elecName, metricName, bandName, 'InteractionSimple'}, cdat6(r, 1:6), {adj6(r)}, cdat6(r, 7:8)]; %#ok<AGROW>
                    end

                    nRefitOK = nRefitOK + 1;
                end
            end
            if mod(e, 10) == 0 || e == nElec
                fprintf('  ... %d/%d electrodes processed (Task=%s, RUN_MODE=%s; %d models exported OK, %d cross-check fails so far)\n', ...
                    e, nElec, TASK_LABEL, RUN_MODE_LABEL, nRefitOK, numel(crossCheckFails));
            end
        end

        ElectrodeAllContrastsT_X18 = local_build_table(allRows, [1 2 3 4 5], [6 7 8 9 10 11 12 13], ...
            {'Electrode','Metric','Band','ContrastFamily','Contrast','Estimate','SE','DF','tStat','pValue','pHolm_withinFamily','CI_Low','CI_High'});

        fprintf('\nRefit summary (Task=%s, RUN_MODE=%s): %d models refit and exported OK (x10 rows each), %d refit/skip failures, %d cross-check mismatches.\n', ...
            TASK_LABEL, RUN_MODE_LABEL, nRefitOK, nRefitFailed, numel(crossCheckFails));

        if ~isempty(crossCheckFails)
            grandCrossCheckFails = grandCrossCheckFails + numel(crossCheckFails);
            fprintf('*** CROSS-CHECK MISMATCHES (refit result does not match X18''s own saved result) ***\n');
            for i = 1:numel(crossCheckFails)
                fprintf('  - %s\n', crossCheckFails{i});
            end
            error('Phase18B:CrossCheckFailed', ...
                ['Task=%s, RUN_MODE=%s: %d refit(s) did not match X18''s own saved OmnibusElectrodeT result ' ...
                 'within tolerance (F: rel %.1e, p: abs %.1e). Stopping before saving any output for this ' ...
                 'Task/RUN_MODE.'], TASK_LABEL, RUN_MODE_LABEL, numel(crossCheckFails), FIT_TOL_F, FIT_TOL_P);
        else
            fprintf('ALL %d refit models cross-checked EXACTLY against X18''s own saved omnibus result (within tolerance, all 3 terms). No mismatches.\n', nCrossCheckOK);
        end

        %% ----------------------------------------------------------------
        %  SAVE -- append-only into X18's own .mat; new .csv
        % ----------------------------------------------------------------
        save(omnibusFile, 'ElectrodeAllContrastsT_X18', '-append');
        writetable(ElectrodeAllContrastsT_X18, fullfile(STATS_DIR, 'electrode_all_contrasts_ungated.csv'));
        fprintf('\nSaved electrode_all_contrasts_ungated.csv (%d rows) and appended ElectrodeAllContrastsT_X18 into %s.\n', ...
            height(ElectrodeAllContrastsT_X18), omnibusFile);

        %% ----------------------------------------------------------------
        %  REPORT
        % ----------------------------------------------------------------
        reportLines = {};
        reportLines{end+1} = repmat('=', 1, 60);
        reportLines{end+1} = sprintf('X18B -- ELECTRODE-LEVEL ALL-CONTRASTS UNGATED REPORT (Task=%s, RUN_MODE=%s)', TASK_LABEL, RUN_MODE_LABEL);
        reportLines{end+1} = repmat('=', 1, 60);
        reportLines{end+1} = '';
        reportLines{end+1} = 'Companion to X18_position_pos1_pos2_cognitive.m -- does not modify X18 and does not';
        reportLines{end+1} = 'change any FDR/Holm/inferential result. Re-derives X18''s own filtered, Position-';
        reportLines{end+1} = 'labeled electrode input (from X04''s raw output + X01''s audit table, same logic X18';
        reportLines{end+1} = 'itself uses) and re-fits the SAME electrode-level Condition*Position models X18';
        reportLines{end+1} = 'already fit, purely to export all 3 contrast families (ConditionMarginal/';
        reportLines{end+1} = 'PositionMarginal/InteractionSimple) for EVERY electrode (ungated, descriptive --';
        reportLines{end+1} = 'needed for a full-scalp topomap). pHolm_withinFamily is NOT a substitute for X18''s';
        reportLines{end+1} = 'own per-term-type BH-FDR-gated pairwise tables.';
        reportLines{end+1} = '';
        reportLines{end+1} = sprintf('Position derivation: %d recordings, %d missing-audit-match, %d ambiguous-match, %d unparseable SEQ (should match X18''s own report for this Task).', ...
            nRecCombos, nMissingAudit, nAmbiguous, nUnparseable);
        reportLines{end+1} = sprintf('Models refit and exported: %d / %d expected (%d electrodes x %d bands x %d metrics; 10 rows per model)', ...
            nRefitOK, nElec*nBands*nMetrics, nElec, nBands, nMetrics);
        reportLines{end+1} = sprintf('Refit/skip failures: %d', nRefitFailed);
        reportLines{end+1} = sprintf('Cross-check against X18''s own saved omnibus result (all 3 terms): %d checked, %d mismatches.', nCrossCheckOK, numel(crossCheckFails));
        if isempty(crossCheckFails)
            reportLines{end+1} = 'Every refit model reproduced X18''s own F-statistics and p-values (all 3 terms) within tolerance.';
        end
        reportLines{end+1} = '';
        reportLines{end+1} = 'Output:';
        reportLines{end+1} = ['  ' fullfile(STATS_DIR, 'electrode_all_contrasts_ungated.csv')];
        reportLines{end+1} = ['  ' omnibusFile ' (ElectrodeAllContrastsT_X18 appended)'];
        reportLines{end+1} = repmat('=', 1, 60);

        reportText = strjoin(reportLines, '\n');
        fprintf('%s\n', reportText);

        fid1 = fopen(fullfile(REPORT_DIR, sprintf('phase18B_electrode_all_contrasts_report_%s_%s.txt', TASK_LABEL, RUN_MODE_LABEL)), 'w');
        fprintf(fid1, '%s\n', reportText);
        fclose(fid1);
    end
end

fprintf('\n=== X18B complete for all RUN_MODEs x both Tasks (Nback, SART -- Full out of scope). ===\n');
if grandCrossCheckFails == 0
    fprintf('No Task/RUN_MODE combination had any cross-check mismatch.\n');
end

%% ======================= LOCAL FUNCTIONS (all after top-level code) =======================
%  local_robust_cellstr, local_position_from_seq,
%  local_fit_condition_position_omnibus, local_find_interaction_name,
%  local_condition_position_contrast, local_position_marginal_contrast,
%  local_apply_contrast, local_condition_contrast are ALL reproduced
%  VERBATIM from X18_position_pos1_pos2_cognitive.m -- this companion
%  script deliberately does not re-derive this contrast arithmetic or the
%  Position-derivation logic independently.

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
            error('Phase18B:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name1);
        end
        cRow(idx1) = cRow(idx1) - 1;
    end
    if ~strcmp(level2, refLevel)
        name2 = sprintf('%s_%s', varName, level2);
        idx2 = find(strcmp(coeffNames, name2));
        if isempty(idx2)
            error('Phase18B:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name2);
        end
        cRow(idx2) = cRow(idx2) + 1;
    end
end

function T = local_build_table(rows, textColIdx, numColIdx, varNames)
    nCols = numel(varNames);
    if size(rows, 2) ~= nCols
        error('Phase18B:BuildTableColumnMismatch', ...
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
            error('Phase18B:BuildTableColumnMismatch', 'Column %d not classified as text or numeric.', j);
        end
    end
    T = table(cols{:}, 'VariableNames', varNames);
end

function c = local_robust_cellstr(col)
% Reproduced VERBATIM from X18_position_pos1_pos2_cognitive.m.
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
% Reproduced VERBATIM from X18_position_pos1_pos2_cognitive.m.
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
% Reproduced VERBATIM from X18_position_pos1_pos2_cognitive.m.
    cand1 = sprintf('Condition_%s:Position_%s', condLevel, posLevel);
    cand2 = sprintf('Position_%s:Condition_%s', posLevel, condLevel);
    if any(strcmp(coeffNames, cand1))
        nm = cand1;
    elseif any(strcmp(coeffNames, cand2))
        nm = cand2;
    else
        error('Phase18B:InteractionNameNotFound', ...
            'Neither "%s" nor "%s" found in CoefficientNames.', cand1, cand2);
    end
end

function cRow = local_condition_position_contrast(coeffNames, refLevel, level1, level2, posSecondLabel, w1, w2)
% Reproduced VERBATIM from X18_position_pos1_pos2_cognitive.m.
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
% Reproduced VERBATIM from X18_position_pos1_pos2_cognitive.m.
    cRow = zeros(1, numel(coeffNames));
    idxPos = find(strcmp(coeffNames, sprintf('Position_%s', posSecondLabel)));
    if isempty(idxPos)
        error('Phase18B:ContrastNameNotFound', 'Coefficient name "Position_%s" not found in model.', posSecondLabel);
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
% Reproduced VERBATIM from X18_position_pos1_pos2_cognitive.m.
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
        fr.errMsg = sprintf('anova(lme) missing one or more expected terms. Terms found: %s', strjoin(aov.Term, ', '));
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
% Reproduced VERBATIM from X18_position_pos1_pos2_cognitive.m.
    est = C * beta;
    se = sqrt(C * covB * C');
    [p, ~, ~, df2c] = coefTest(lme, C, 0, 'DFMethod', 'satterthwaite');
    tstat = est / se;
    tcrit = tinv(0.975, df2c);
    ciLow = est - tcrit * se;
    ciHigh = est + tcrit * se;
end
