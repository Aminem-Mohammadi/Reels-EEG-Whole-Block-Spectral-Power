%% ========================================================================
% X00_noise_audit_cognitive.m
%
% Purpose: pre-statistics data-quality audit. Computes an objective,
% pre-specified composite noise score per subject, from four metrics, to
% produce a candidate "noisy subject" exclusion list for the sensitivity
% re-run (see EXCLUSION_CONFIG_cognitive.m). This script does not exclude
% anyone itself and does not fit any inferential model; it only produces
% the candidate list, which is reviewed before use.
%
% Method: four metrics, each computed from the N-back + SART recordings
% only (never Full, since Full's extra non-task time would bias
% amplitude/spectral estimates for reasons unrelated to signal quality):
%   M1: window-rejection rate (100 - mean PercentValid), from X04's
%       phase3_recording_summary.mat / RecordingQCT.
%   M2: high-frequency power ratio, Power(20-45Hz) / Power(1-45Hz), from
%       X04's GlobalBandPowerT -- the classic residual muscle/EMG
%       contamination signature that can survive ICA/ICLabel cleaning.
%   M3: |z-score| of mean log(broadband 1-45Hz power), two-sided (both
%       abnormally high power, e.g. residual movement/muscle artifact,
%       and abnormally low power, e.g. excessive attenuation or a
%       gain/calibration problem, are data-quality concerns).
%   M4: NBadChannels + NBadICs, averaged per subject from
%       Preprocessing_QC.xlsx (supplementary/corroborating only).
% Composite score = unweighted mean of {z(M1), z(M2), |z(M3)|, z(M4)}
% (3-metric mean of {z(M1), z(M2), |z(M3)|} if Preprocessing_QC.xlsx is
% not found). Threshold is pre-specified (composite z > 2.0) and is not
% tuned after seeing who it flags. FORCED_EXCLUSION_SUBJECTS below allows
% forcing specific subjects onto the candidate list for an independent,
% non-audit reason, disclosed in FORCED_EXCLUSION_RATIONALE; empty by
% default. If
% the z > 2.0 rule (plus any forced subjects) flags more than ~3-4
% subjects, the
% script prints a warning rather than auto-tightening the threshold,
% since adjusting the threshold after seeing the names would reintroduce
% the exact circularity this design is meant to avoid; the ranked list is
% reviewed and the exclusion decision made manually.
%
% This script does not fit any LMM, does not touch FDR/Holm, does not
% compute any inferential statistic, and does not exclude anyone from any
% existing output -- it only produces the candidate list above.
%
% Inputs (re-run X03 and X04 first if any of these is missing):
%   <ANALYSIS_ROOT>/04_subject_level/phase3_recording_summary.mat  (RecordingQCT)
%   <ANALYSIS_ROOT>/04_subject_level/phase3_subject_level_bandpower.mat (GlobalBandPowerT, PARAMS3)
%   Preprocessing_QC.xlsx (path set below; supplementary, optional)
%
% Outputs:
%   <ANALYSIS_ROOT>/00_noise_audit/
%       noise_audit_subject_scores.xlsx / .mat
%       noise_audit_parameters.mat
%   <ANALYSIS_ROOT>/09_reports/
%       noise_audit_report.txt
%
% RUN_MODE: not applicable; this script runs once, upstream of the
% full_sample / sensitivity_excl_noisy split it feeds.
% ========================================================================

clear; clc;

%% ------------------------------------------------------------------
%  0. USER-EDITABLE SETTINGS
% ------------------------------------------------------------------

config;
DATA_ROOT_FOR_SIBLING = DATA_DIR;
[DATA_ROOT_PARENT, ~] = fileparts(DATA_ROOT_FOR_SIBLING);
ANALYSIS_ROOT = OUT_DIR;

SUBJ_DIR    = fullfile(ANALYSIS_ROOT, '04_subject_level');
AUDIT_DIR   = fullfile(ANALYSIS_ROOT, '00_noise_audit');
REPORT_DIR  = fullfile(ANALYSIS_ROOT, '09_reports');
if ~exist(AUDIT_DIR, 'dir');  mkdir(AUDIT_DIR);  end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

% EDIT THIS PATH if Preprocessing_QC.xlsx lives somewhere else on your
% machine. This input is supplementary only (M4); the script still runs
% (with a 3-metric composite, clearly disclosed) if it is not found here.
PREPROCESSING_QC_FILE = fullfile(DATA_ROOT_PARENT, 'Preprocessing_QC.xlsx');

FORCED_EXCLUSION_SUBJECTS = [];
FORCED_EXCLUSION_RATIONALE = struct();

COMPOSITE_Z_THRESHOLD = 2.0;
EXCLUSION_COUNT_MONITORING_CEILING = 4;   % Advisory only -- see header.

TASK_TYPES_FOR_AUDIT = {'nback', 'sart'};   % Full deliberately excluded -- see header.

fprintf('========================================================\n');
fprintf('X00 -- PRE-STATISTICS NOISE AUDIT (COGNITIVE-TASK DATA)\n');
fprintf('========================================================\n');
fprintf('ANALYSIS_ROOT : %s\n', ANALYSIS_ROOT);
fprintf('OUTPUT        : %s\n', AUDIT_DIR);
fprintf('--------------------------------------------------------\n\n');

%% ------------------------------------------------------------------
%  1. LOAD X04 OUTPUTS (source of truth; nothing recomputed from PSD)
% ------------------------------------------------------------------

qcFile = fullfile(SUBJ_DIR, 'phase3_recording_summary.mat');
bpFile = fullfile(SUBJ_DIR, 'phase3_subject_level_bandpower.mat');

if ~exist(qcFile, 'file')
    error('X00Noise:MissingInput', 'Required input file not found: %s\nRun X04_aggregate_cognitive_subject_level.m first.', qcFile);
end
if ~exist(bpFile, 'file')
    error('X00Noise:MissingInput', 'Required input file not found: %s\nRun X04_aggregate_cognitive_subject_level.m first.', bpFile);
end

S1 = load(qcFile, 'RecordingQCT');
RecordingQCT_ALL = S1.RecordingQCT;
clear S1;

S2 = load(bpFile, 'GlobalBandPowerT', 'PARAMS3');
GlobalBandPowerT_ALL = S2.GlobalBandPowerT;
PARAMS3 = S2.PARAMS3;
clear S2;

REQUIRED_BANDS_M2 = {'Delta', 'Theta', 'Alpha', 'Beta', 'Beta_high_exploratory', 'Gamma_low_exploratory'};
definedBands = fieldnames(PARAMS3.BandDefinitionsHz);
missingBands = REQUIRED_BANDS_M2(~ismember(REQUIRED_BANDS_M2, definedBands));
if ~isempty(missingBands)
    error('X00Noise:MissingBands', ...
        ['This script needs Beta_high_exploratory (added alongside Gamma) to compute the ' ...
         'high-frequency-ratio metric (M2). Missing band(s): %s.\n' ...
         'Re-run X03_compute_cognitive_psd.m and X04_aggregate_cognitive_subject_level.m ' ...
         'first (they now include Beta_low_exploratory/Beta_high_exploratory), then re-run this script.'], ...
        strjoin(missingBands, ', '));
end

%% ------------------------------------------------------------------
%  2. LOAD Preprocessing_QC.xlsx (supplementary; graceful degrade)
% ------------------------------------------------------------------

haveQCFile = exist(PREPROCESSING_QC_FILE, 'file') == 2;
if haveQCFile
    try
        PreprocQCT = readtable(PREPROCESSING_QC_FILE);
        requiredQCCols = {'Subject', 'NBadChannels', 'NBadICs'};
        missingQCCols = requiredQCCols(~ismember(requiredQCCols, PreprocQCT.Properties.VariableNames));
        if ~isempty(missingQCCols)
            warning('X00Noise:QCFileColumnsMissing', ...
                'Preprocessing_QC.xlsx is missing column(s): %s -- proceeding WITHOUT M4 (3-metric composite).', ...
                strjoin(missingQCCols, ', '));
            haveQCFile = false;
        end
    catch ME
        warning('X00Noise:QCFileReadFailed', 'Could not read %s (%s) -- proceeding WITHOUT M4 (3-metric composite).', ...
            PREPROCESSING_QC_FILE, ME.message);
        haveQCFile = false;
    end
else
    warning('X00Noise:QCFileNotFound', ...
        'Preprocessing_QC.xlsx not found at %s -- proceeding WITHOUT M4 (3-metric composite). Edit PREPROCESSING_QC_FILE at the top of this script if it lives elsewhere.', ...
        PREPROCESSING_QC_FILE);
end

if haveQCFile
    fprintf('Preprocessing_QC.xlsx loaded: %d rows. M4 (NBadChannels+NBadICs) WILL be included.\n\n', height(PreprocQCT));
else
    fprintf('Preprocessing_QC.xlsx NOT used. Composite score will be based on M1+M2+M3 only (disclosed in output).\n\n');
end

%% ------------------------------------------------------------------
%  3. RESTRICT TO N-back + SART (never Full -- see header)
% ------------------------------------------------------------------

qcMask = ismember(RecordingQCT_ALL.TaskType, TASK_TYPES_FOR_AUDIT);
RecordingQCT = RecordingQCT_ALL(qcMask, :);

gbpMask = ismember(GlobalBandPowerT_ALL.TaskType, TASK_TYPES_FOR_AUDIT);
GlobalBandPowerT = GlobalBandPowerT_ALL(gbpMask, :);

fprintf('Restricted to TaskType in {%s}: %d RecordingQCT rows, %d GlobalBandPowerT rows.\n\n', ...
    strjoin(TASK_TYPES_FOR_AUDIT, ','), height(RecordingQCT), height(GlobalBandPowerT));

subjNum_qc  = local_subject_numeric(RecordingQCT.Subject);
subjNum_gbp = local_subject_numeric(GlobalBandPowerT.Subject);
allSubjects = unique([subjNum_qc; subjNum_gbp]);
nSubj = numel(allSubjects);

integrityIssues = {};

%% ------------------------------------------------------------------
%  4. M1: window-rejection rate, per subject
% ------------------------------------------------------------------

M1_RejectionRatePct = nan(nSubj, 1);
for si = 1:nSubj
    idx = subjNum_qc == allSubjects(si);
    if ~any(idx)
        integrityIssues{end+1} = sprintf('Subject %g: no RecordingQCT rows found for Nback/SART.', allSubjects(si)); %#ok<AGROW>
        continue;
    end
    M1_RejectionRatePct(si) = 100 - mean(RecordingQCT.PercentValid(idx), 'omitnan');
end

%% ------------------------------------------------------------------
%  5. M2: high-frequency power ratio, per subject
%     Power(20-45Hz) / Power(1-45Hz) = (Beta_high_exploratory+Gamma_low_exploratory)
%                                       / (Delta+Theta+Alpha+Beta+Gamma_low_exploratory)
% ------------------------------------------------------------------

% Build a per-row ratio, then average per subject (each row = one
% Subject x BlockNumber x Condition recording, so this weights every
% recording equally, matching the "one row per recording" convention
% used everywhere else in this pipeline).
uniqueBandRows = unique(GlobalBandPowerT(:, {'Subject','TaskType','BlockNumber','Condition'}), 'rows');
nRowsU = height(uniqueBandRows);
rowHighFreq = nan(nRowsU, 1);
rowTotal    = nan(nRowsU, 1);
rowSubjNum  = local_subject_numeric(uniqueBandRows.Subject);

for r = 1:nRowsU
    m = strcmp(GlobalBandPowerT.Subject, uniqueBandRows.Subject{r}) & ...
        strcmp(GlobalBandPowerT.TaskType, uniqueBandRows.TaskType{r}) & ...
        GlobalBandPowerT.BlockNumber == uniqueBandRows.BlockNumber(r) & ...
        strcmp(GlobalBandPowerT.Condition, uniqueBandRows.Condition{r});
    thisRows = GlobalBandPowerT(m, :);

    getBand = @(bn) thisRows.AbsolutePower(strcmp(thisRows.Band, bn));
    delta = getBand('Delta'); theta = getBand('Theta'); alpha = getBand('Alpha');
    beta  = getBand('Beta');  gammaLow = getBand('Gamma_low_exploratory'); betaHigh = getBand('Beta_high_exploratory');

    if numel(delta)~=1 || numel(theta)~=1 || numel(alpha)~=1 || numel(beta)~=1 || numel(gammaLow)~=1 || numel(betaHigh)~=1
        integrityIssues{end+1} = sprintf('Subject %s Task%s Block%g Cond%s: could not find exactly one row per required band for M2 -- row skipped.', ...
            uniqueBandRows.Subject{r}, uniqueBandRows.TaskType{r}, uniqueBandRows.BlockNumber(r), uniqueBandRows.Condition{r}); %#ok<AGROW>
        continue;
    end

    rowTotal(r)    = delta + theta + alpha + beta + gammaLow;
    rowHighFreq(r) = betaHigh + gammaLow;
end
rowRatio = rowHighFreq ./ rowTotal;

M2_HighFreqRatio = nan(nSubj, 1);
for si = 1:nSubj
    idx = rowSubjNum == allSubjects(si);
    if any(idx)
        M2_HighFreqRatio(si) = mean(rowRatio(idx), 'omitnan');
    else
        integrityIssues{end+1} = sprintf('Subject %g: no GlobalBandPowerT rows found for M2.', allSubjects(si)); %#ok<AGROW>
    end
end

%% ------------------------------------------------------------------
%  6. M3: |z-score| of mean log(broadband 1-45Hz power), per subject
% ------------------------------------------------------------------

rowLogTotal = log(rowTotal);
M3_LogBroadbandPower = nan(nSubj, 1);
for si = 1:nSubj
    idx = rowSubjNum == allSubjects(si);
    if any(idx)
        M3_LogBroadbandPower(si) = mean(rowLogTotal(idx), 'omitnan');
    end
end

%% ------------------------------------------------------------------
%  7. M4: NBadChannels + NBadICs, per subject (supplementary, optional)
% ------------------------------------------------------------------

M4_QCBadCount = nan(nSubj, 1);
if haveQCFile
    qcSubjNum = double(PreprocQCT.Subject);
    for si = 1:nSubj
        idx = qcSubjNum == allSubjects(si);
        if any(idx)
            M4_QCBadCount(si) = mean(PreprocQCT.NBadChannels(idx) + PreprocQCT.NBadICs(idx), 'omitnan');
        else
            integrityIssues{end+1} = sprintf('Subject %g: no Preprocessing_QC.xlsx rows found for M4.', allSubjects(si)); %#ok<AGROW>
        end
    end
end

%% ------------------------------------------------------------------
%  8. Z-SCORE EACH METRIC ACROSS SUBJECTS; COMPOSITE
% ------------------------------------------------------------------

z1 = local_zscore(M1_RejectionRatePct);
z2 = local_zscore(M2_HighFreqRatio);
z3 = abs(local_zscore(M3_LogBroadbandPower));
if haveQCFile
    z4 = local_zscore(M4_QCBadCount);
    CompositeZ = mean([z1, z2, z3, z4], 2, 'omitnan');
    nMetricsUsed = 4;
else
    z4 = nan(nSubj, 1);
    CompositeZ = mean([z1, z2, z3], 2, 'omitnan');
    nMetricsUsed = 3;
end

%% ------------------------------------------------------------------
%  9. FLAGGING (pre-specified threshold + forced-exclusion list)
% ------------------------------------------------------------------

FlaggedByScore = CompositeZ > COMPOSITE_Z_THRESHOLD;
ForcedExclusion = ismember(allSubjects, FORCED_EXCLUSION_SUBJECTS);
FinalExclusion = FlaggedByScore | ForcedExclusion;

nFlaggedTotal = sum(FinalExclusion);
nFlaggedByScoreOnly = sum(FlaggedByScore & ~ForcedExclusion);

%% ------------------------------------------------------------------
%  10. ASSEMBLE + SAVE
% ------------------------------------------------------------------

ScoreT = table(allSubjects, M1_RejectionRatePct, M2_HighFreqRatio, M3_LogBroadbandPower, M4_QCBadCount, ...
    z1, z2, z3, z4, CompositeZ, FlaggedByScore, ForcedExclusion, FinalExclusion, ...
    'VariableNames', {'Subject', 'M1_RejectionRatePct', 'M2_HighFreqRatio', 'M3_LogBroadbandPower', 'M4_QCBadCount', ...
    'z_M1', 'z_M2', 'z_M3_abs', 'z_M4', 'CompositeZ', 'FlaggedByScore', 'ForcedExclusion', 'FinalExclusion'});
ScoreT = sortrows(ScoreT, 'CompositeZ', 'descend');

fprintf('--- Noise audit scores, sorted by CompositeZ (descending) ---\n');
disp(ScoreT);

writetable(ScoreT, fullfile(AUDIT_DIR, 'noise_audit_subject_scores.xlsx'));
save(fullfile(AUDIT_DIR, 'noise_audit_subject_scores.mat'), 'ScoreT');

PARAMS_NOISE = struct();
PARAMS_NOISE.ScriptName = 'X00_noise_audit_cognitive.m';
PARAMS_NOISE.GeneratedOn = datestr(now); %#ok<TNOW1,DATST>
PARAMS_NOISE.Purpose = 'Pre-specified, audit-derived composite noise score per subject, computed BEFORE any inferential model, to produce a disclosed candidate exclusion list for a later sensitivity re-run.';
PARAMS_NOISE.TasksUsed = TASK_TYPES_FOR_AUDIT;
PARAMS_NOISE.MetricsUsed = nMetricsUsed;
PARAMS_NOISE.M1Definition = '100 - mean(PercentValid) across Nback+SART recordings (window-rejection rate).';
PARAMS_NOISE.M2Definition = 'mean over recordings of Power(20-45Hz)/Power(1-45Hz) = (Beta_high_exploratory+Gamma_low_exploratory)/(Delta+Theta+Alpha+Beta+Gamma_low_exploratory).';
PARAMS_NOISE.M3Definition = '|z-score| of mean log(broadband 1-45Hz power); two-sided by design (both high and low are quality concerns).';
PARAMS_NOISE.M4Definition = 'mean(NBadChannels+NBadICs) from Preprocessing_QC.xlsx, supplementary only.';
PARAMS_NOISE.M4Available = haveQCFile;
PARAMS_NOISE.CompositeFormula = 'Unweighted mean of {z(M1), z(M2), |z(M3)|, z(M4)} (equal weights, disclosed); 3-metric mean if M4 unavailable.';
PARAMS_NOISE.CompositeZThreshold = COMPOSITE_Z_THRESHOLD;
PARAMS_NOISE.ThresholdPreSpecified = true;
PARAMS_NOISE.ForcedExclusionSubjects = FORCED_EXCLUSION_SUBJECTS;
PARAMS_NOISE.ForcedExclusionRationale = FORCED_EXCLUSION_RATIONALE;
PARAMS_NOISE.ExclusionCountMonitoringCeiling = EXCLUSION_COUNT_MONITORING_CEILING;
PARAMS_NOISE.ExclusionCountCeilingIsAdvisoryNotAutoTuned = true;
save(fullfile(AUDIT_DIR, 'noise_audit_parameters.mat'), 'PARAMS_NOISE');

%% ------------------------------------------------------------------
%  11. REPORT
% ------------------------------------------------------------------

reportLines = {};
reportLines{end+1} = repmat('=', 1, 60);
reportLines{end+1} = 'X00 -- PRE-STATISTICS NOISE AUDIT (COGNITIVE-TASK DATA)';
reportLines{end+1} = repmat('=', 1, 60);
reportLines{end+1} = '';
reportLines{end+1} = sprintf('Metrics used: %d (%s)', nMetricsUsed, iif(haveQCFile, 'M1+M2+M3+M4', 'M1+M2+M3 -- Preprocessing_QC.xlsx not used, see warning above'));
reportLines{end+1} = sprintf('Computed from TaskType in {%s} only (never Full).', strjoin(TASK_TYPES_FOR_AUDIT, ','));
reportLines{end+1} = sprintf('Composite z-score threshold (pre-specified): %.1f', COMPOSITE_Z_THRESHOLD);
reportLines{end+1} = '';
reportLines{end+1} = 'Ranked scores (all subjects):';
for r = 1:height(ScoreT)
    reportLines{end+1} = sprintf('  Subject %-4g CompositeZ=%7.3f | M1=%.2f%% z=%.2f | M2=%.4f z=%.2f | M3(log)=%.3f |z|=%.2f | M4=%s z=%s | Flag=%d Forced=%d', ...
        ScoreT.Subject(r), ScoreT.CompositeZ(r), ScoreT.M1_RejectionRatePct(r), ScoreT.z_M1(r), ...
        ScoreT.M2_HighFreqRatio(r), ScoreT.z_M2(r), ScoreT.M3_LogBroadbandPower(r), ScoreT.z_M3_abs(r), ...
        local_numstr(ScoreT.M4_QCBadCount(r)), local_numstr(ScoreT.z_M4(r)), ...
        ScoreT.FlaggedByScore(r), ScoreT.ForcedExclusion(r)); %#ok<AGROW>
end
reportLines{end+1} = '';
reportLines{end+1} = sprintf('FINAL CANDIDATE EXCLUSION LIST (%d subjects): %s', nFlaggedTotal, ...
    strjoin(cellstr(num2str(ScoreT.Subject(ScoreT.FinalExclusion), '%g')), ', '));
reportLines{end+1} = sprintf('  - Flagged by composite score (z>%.1f): %d subject(s)', COMPOSITE_Z_THRESHOLD, nFlaggedByScoreOnly);
if isempty(FORCED_EXCLUSION_SUBJECTS)
    reportLines{end+1} = '  - Force-included for an independent, non-audit reason: (none)';
else
    reportLines{end+1} = sprintf('  - Force-included for an independent, non-audit reason: %s', strjoin(cellstr(num2str(FORCED_EXCLUSION_SUBJECTS(:), '%g')), ', '));
end
reportLines{end+1} = '';
if nFlaggedTotal > EXCLUSION_COUNT_MONITORING_CEILING
    reportLines{end+1} = sprintf(['*** MONITORING CEILING EXCEEDED: %d subjects flagged, ceiling was %d. ' ...
        'Threshold was NOT auto-adjusted (see header rationale) -- please review the ranked list above and decide manually. ***'], ...
        nFlaggedTotal, EXCLUSION_COUNT_MONITORING_CEILING);
else
    reportLines{end+1} = sprintf('Total flagged (%d) is within the %d-subject monitoring ceiling.', nFlaggedTotal, EXCLUSION_COUNT_MONITORING_CEILING);
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
reportLines{end+1} = 'Output files written to:';
reportLines{end+1} = ['  ' AUDIT_DIR];
reportLines{end+1} = repmat('=', 1, 60);

reportText = strjoin(reportLines, '\n');
fprintf('%s\n', reportText);
fid = fopen(fullfile(REPORT_DIR, 'noise_audit_report.txt'), 'w');
fprintf(fid, '%s\n', reportText);
fclose(fid);

fprintf('\n=== X00 complete. Review noise_audit_subject_scores.xlsx before this exclusion list is used anywhere downstream. ===\n');

%% ======================= LOCAL FUNCTIONS =======================

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

function z = local_zscore(x)
    mu = mean(x, 'omitnan');
    sigma = std(x, 0, 'omitnan');
    if sigma == 0 || isnan(sigma)
        z = zeros(size(x));
    else
        z = (x - mu) / sigma;
    end
end

function s = local_numstr(v)
    if isnan(v)
        s = 'NA';
    else
        s = sprintf('%.2f', v);
    end
end

function out = iif(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end
