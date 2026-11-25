%% ========================================================================
% X00_noise_audit_intervention.m
%
% Composite noise score per participant for the intervention recordings.
% Quality indicators per recording: percentage of rejected windows, high-frequency power
% ratio P(20-45 Hz)/P(1-45 Hz), absolute standardized deviation of log broadband power, and
% the number of quality problems in the preprocessing log (bad channels + bad ICs).
% Indicators are z-scored within condition, averaged per participant and z-scored across
% participants (CompositeZ). Threshold z > 2.00; 1.75-2.00 is listed as borderline.
% Inputs : 02_metadata, 04_subject_level.
% Outputs: 00_noise_audit/ (noise_audit_subject_scores.xlsx, noise_audit_recording_scores.xlsx),
%          09_reports/noise_audit_report.txt.
% ========================================================================

clear; clc;

%% 0. SETTINGS
config;
ANALYSIS_ROOT = OUT_DIR;
SUBJ_DIR   = fullfile(ANALYSIS_ROOT, '04_subject_level');
META_FILE  = fullfile(ANALYSIS_ROOT, '02_metadata', 'intervention_metadata.mat');
AUDIT_DIR  = fullfile(ANALYSIS_ROOT, '00_noise_audit');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');
if ~exist(AUDIT_DIR, 'dir');  mkdir(AUDIT_DIR);  end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

COMPOSITE_Z_THRESHOLD = 2.0;
BORDERLINE_Z = 1.75;                 % listed only, never flagged automatically
EXCLUSION_COUNT_MONITORING_CEILING = 4;
CONDITION_ORDER = {'E', 'N', 'R'};

fprintf('========================================================\n');
fprintf('X00 -- PRE-STATISTICS NOISE AUDIT (INTERVENTION DATA)\n');
fprintf('========================================================\n');

%% 1. LOAD X04 OUTPUTS
qcFile = fullfile(SUBJ_DIR, 'phase3_recording_summary.mat');
bpFile = fullfile(SUBJ_DIR, 'phase3_subject_level_bandpower.mat');
if ~exist(qcFile, 'file') || ~exist(bpFile, 'file')
    error('X00Noise:MissingInput', 'X04 outputs not found in %s. Run X03 and X04 first.', SUBJ_DIR);
end
S1 = load(qcFile, 'RecordingQCT');            RecordingQCT = S1.RecordingQCT; clear S1;
S2 = load(bpFile, 'GlobalBandPowerT', 'PARAMS3');
GlobalBandPowerT = S2.GlobalBandPowerT; PARAMS3 = S2.PARAMS3; clear S2;

REQUIRED_BANDS_M2 = {'Delta', 'Theta', 'Alpha', 'Beta', 'Beta_high_exploratory', 'Gamma_low_exploratory'};
missingBands = REQUIRED_BANDS_M2(~ismember(REQUIRED_BANDS_M2, fieldnames(PARAMS3.BandDefinitionsHz)));
if ~isempty(missingBands)
    error('X00Noise:MissingBands', ['Missing band(s): %s. M2 needs Beta_high_exploratory, which X03 ' ...
        'adds. Re-run X03 and X04 first.'], strjoin(missingBands, ', '));
end

integrityIssues = {};

%% 2. ONE ROW PER RECORDING: M1, M2, M3
recT = unique(GlobalBandPowerT(:, {'Subject', 'InterventionNumber', 'Condition'}), 'rows');
nRec = height(recT);
recSubj = local_subject_numeric(recT.Subject);
M1 = nan(nRec, 1); M2 = nan(nRec, 1); M3 = nan(nRec, 1); M4 = nan(nRec, 1);
qcSubj = local_subject_numeric(RecordingQCT.Subject);
for r = 1:nRec
    % M1
    q = qcSubj == recSubj(r) & RecordingQCT.InterventionNumber == recT.InterventionNumber(r) & ...
        strcmp(RecordingQCT.Condition, recT.Condition{r});
    if sum(q) == 1
        M1(r) = 100 - RecordingQCT.PercentValid(q);
    else
        integrityIssues{end+1} = sprintf('Sub%g Int%g %s: %d RecordingQCT rows (expected 1) -- M1 missing.', ...
            recSubj(r), recT.InterventionNumber(r), recT.Condition{r}, sum(q)); %#ok<AGROW>
    end
    % M2, M3
    m = strcmp(GlobalBandPowerT.Subject, recT.Subject{r}) & ...
        GlobalBandPowerT.InterventionNumber == recT.InterventionNumber(r) & ...
        strcmp(GlobalBandPowerT.Condition, recT.Condition{r});
    rowsHere = GlobalBandPowerT(m, :);
    getBand = @(bn) rowsHere.AbsolutePower(strcmp(rowsHere.Band, bn));
    d = getBand('Delta'); t = getBand('Theta'); a = getBand('Alpha'); b = getBand('Beta');
    gl = getBand('Gamma_low_exploratory'); bh = getBand('Beta_high_exploratory');
    if all([numel(d) numel(t) numel(a) numel(b) numel(gl) numel(bh)] == 1)
        total = d + t + a + b + gl;
        M2(r) = (bh + gl) / total;
        M3(r) = log(total);
    else
        integrityIssues{end+1} = sprintf('Sub%g Int%g %s: band rows missing -- M2/M3 missing.', ...
            recSubj(r), recT.InterventionNumber(r), recT.Condition{r}); %#ok<AGROW>
    end
end

%% 3. M4 FROM PHASE 1 METADATA (optional)
haveM4 = false;
if exist(META_FILE, 'file')
    SM = load(META_FILE, 'T'); MetaT = SM.T; clear SM;
    if all(ismember({'QC_NBadChannels', 'QC_NBadICs'}, MetaT.Properties.VariableNames))
        metaSubj = local_subject_numeric(MetaT.Subject);
        nbc = local_to_double(MetaT.QC_NBadChannels);
        nbi = local_to_double(MetaT.QC_NBadICs);
        metaInt = local_to_double(MetaT.InterventionNumber);
        for r = 1:nRec
            k = metaSubj == recSubj(r) & metaInt == recT.InterventionNumber(r);
            if sum(k) == 1
                M4(r) = nbc(k) + nbi(k);
            end
        end
        haveM4 = sum(isfinite(M4)) >= 0.8 * nRec;
        if ~haveM4
            warning('X00Noise:M4Sparse', 'NBadChannels/NBadICs available for only %d/%d recordings -- M4 NOT used.', sum(isfinite(M4)), nRec);
        end
    else
        warning('X00Noise:M4Missing', 'QC_NBadChannels / QC_NBadICs not in metadata T -- M4 NOT used.');
    end
else
    warning('X00Noise:NoMeta', 'Metadata file not found (%s) -- M4 NOT used.', META_FILE);
end

%% 4. RECORDING-LEVEL Z, WITHIN CONDITION (M1-M3); OVERALL (M4)
zr1 = nan(nRec, 1); zr2 = nan(nRec, 1); zr3 = nan(nRec, 1);
for c = 1:numel(CONDITION_ORDER)
    k = strcmp(recT.Condition, CONDITION_ORDER{c});
    zr1(k) = local_zscore(M1(k));
    zr2(k) = local_zscore(M2(k));
    zr3(k) = local_zscore(M3(k));
end
zr4 = local_zscore(M4);

RecordingScoreT = table(recSubj, recT.InterventionNumber, recT.Condition, M1, M2, M3, M4, zr1, zr2, zr3, zr4, ...
    'VariableNames', {'Subject', 'InterventionNumber', 'Condition', 'M1_RejectionRatePct', 'M2_HighFreqRatio', ...
    'M3_LogBroadbandPower', 'M4_QCBadCount', 'zWithinCond_M1', 'zWithinCond_M2', 'zWithinCond_M3', 'z_M4'});
RecordingScoreT = sortrows(RecordingScoreT, {'Subject', 'InterventionNumber'});

%% 5. SUBJECT LEVEL: mean of the subject's recordings, then z across subjects
allSubjects = unique(recSubj);
nSubj = numel(allSubjects);
s1 = nan(nSubj, 1); s2 = nan(nSubj, 1); s3 = nan(nSubj, 1); s4 = nan(nSubj, 1);
condPair = cell(nSubj, 1); worstRec = cell(nSubj, 1);
for i = 1:nSubj
    k = recSubj == allSubjects(i);
    s1(i) = mean(zr1(k), 'omitnan'); s2(i) = mean(zr2(k), 'omitnan');
    s3(i) = mean(zr3(k), 'omitnan'); s4(i) = mean(zr4(k), 'omitnan');
    condPair{i} = strjoin(sort(recT.Condition(k))', '&');
    [~, w] = max(M1(k)); ints = recT.InterventionNumber(k);
    if ~isempty(w); worstRec{i} = sprintf('Int%g (%.1f%% rejected)', ints(w), max(M1(k))); else; worstRec{i} = ''; end
end
z1 = local_zscore(s1);
z2 = local_zscore(s2);
z3 = abs(local_zscore(s3));   % two-sided by design
if haveM4
    z4 = local_zscore(s4);
    CompositeZ = mean([z1, z2, z3, z4], 2, 'omitnan');
    nMetricsUsed = 4;
else
    z4 = nan(nSubj, 1);
    CompositeZ = mean([z1, z2, z3], 2, 'omitnan');
    nMetricsUsed = 3;
end
FlaggedByScore = CompositeZ > COMPOSITE_Z_THRESHOLD;
Borderline = CompositeZ > BORDERLINE_Z & ~FlaggedByScore;

ScoreT = table(allSubjects, condPair, z1, z2, z3, z4, CompositeZ, FlaggedByScore, Borderline, worstRec, ...
    'VariableNames', {'Subject', 'ConditionPair', 'z_M1', 'z_M2', 'z_M3_abs', 'z_M4', 'CompositeZ', ...
    'FlaggedByScore', 'Borderline_1p75', 'WorstRecordingM1'});
ScoreT = sortrows(ScoreT, 'CompositeZ', 'descend');
disp(ScoreT);

writetable(ScoreT, fullfile(AUDIT_DIR, 'noise_audit_subject_scores.xlsx'));
writetable(RecordingScoreT, fullfile(AUDIT_DIR, 'noise_audit_recording_scores.xlsx'));
save(fullfile(AUDIT_DIR, 'noise_audit_subject_scores.mat'), 'ScoreT', 'RecordingScoreT');

PARAMS_NOISE = struct();
PARAMS_NOISE.ScriptName = 'X00_noise_audit_intervention.m';
PARAMS_NOISE.GeneratedOn = datestr(now); %#ok<TNOW1,DATST>
PARAMS_NOISE.MetricsUsed = nMetricsUsed;
PARAMS_NOISE.Standardisation = 'recording value -> z within condition (M1-M3) or over all recordings (M4) -> mean per subject -> z across subjects; |z| for M3';
PARAMS_NOISE.M2Definition = '(Beta_high_exploratory + Gamma_low_exploratory) / (Delta+Theta+Alpha+Beta+Gamma_low_exploratory) = P(20-45)/P(1-45)';
PARAMS_NOISE.CompositeFormula = 'unweighted mean of {z(M1), z(M2), |z(M3)|, z(M4)}; 3-metric mean if M4 unavailable';
PARAMS_NOISE.CompositeZThreshold = COMPOSITE_Z_THRESHOLD;
PARAMS_NOISE.BorderlineZListedOnly = BORDERLINE_Z;
PARAMS_NOISE.ForcedExclusionSubjects = [];
PARAMS_NOISE.ExclusionCountMonitoringCeiling = EXCLUSION_COUNT_MONITORING_CEILING;
save(fullfile(AUDIT_DIR, 'noise_audit_parameters.mat'), 'PARAMS_NOISE');

%% 6. REPORT
L = {};
L{end+1} = repmat('=', 1, 70);
L{end+1} = 'X00 -- PRE-STATISTICS NOISE AUDIT (INTERVENTION DATA)';
L{end+1} = repmat('=', 1, 70);
L{end+1} = sprintf('Metrics used: %d (%s)', nMetricsUsed, local_iif(haveM4, 'M1+M2+M3+M4', 'M1+M2+M3; M4 not available'));
L{end+1} = 'Standardisation: z WITHIN condition per recording -> subject mean -> z across subjects (see header).';
L{end+1} = sprintf('Pre-specified threshold: CompositeZ > %.2f (borderline %.2f-%.2f listed, not flagged).', ...
    COMPOSITE_Z_THRESHOLD, BORDERLINE_Z, COMPOSITE_Z_THRESHOLD);
L{end+1} = '';
L{end+1} = 'Ranked subjects:';
for r = 1:height(ScoreT)
    L{end+1} = sprintf('  Sub %-3g [%s] CompositeZ=%6.2f | zM1=%5.2f zM2=%5.2f |zM3|=%5.2f zM4=%s | worst M1: %s%s', ...
        ScoreT.Subject(r), ScoreT.ConditionPair{r}, ScoreT.CompositeZ(r), ScoreT.z_M1(r), ScoreT.z_M2(r), ...
        ScoreT.z_M3_abs(r), local_numstr(ScoreT.z_M4(r)), ScoreT.WorstRecordingM1{r}, ...
        local_iif(ScoreT.FlaggedByScore(r), '  <-- FLAGGED', local_iif(ScoreT.Borderline_1p75(r), '  <-- borderline', ''))); %#ok<SAGROW>
end
flagged = ScoreT.Subject(ScoreT.FlaggedByScore);
border  = ScoreT.Subject(ScoreT.Borderline_1p75);
L{end+1} = '';
L{end+1} = sprintf('FLAGGED (z > %.2f): %s', COMPOSITE_Z_THRESHOLD, local_list(flagged));
L{end+1} = sprintf('Borderline (%.2f < z <= %.2f, NOT flagged): %s', BORDERLINE_Z, COMPOSITE_Z_THRESHOLD, local_list(border));
if numel(flagged) > EXCLUSION_COUNT_MONITORING_CEILING
    L{end+1} = sprintf('*** %d subjects flagged, above the monitoring ceiling of %d. Threshold NOT changed automatically -- review by hand. ***', ...
        numel(flagged), EXCLUSION_COUNT_MONITORING_CEILING);
end
L{end+1} = '';
L{end+1} = 'Suggested line for EXCLUSION_CONFIG_intervention.m (flagged subjects + any subject you keep';
L{end+1} = 'for an independent, documented reason, e.g. 7 for its Int1 boundary loss):';
L{end+1} = sprintf('  EXCLUDE_NOISY_SUBJECTS = [%s];', strjoin(cellstr(num2str(flagged(:))), ', '));
L{end+1} = '';
if isempty(integrityIssues)
    L{end+1} = 'Integrity issues: none.';
else
    L{end+1} = '--- Integrity issues ---';
    for i = 1:numel(integrityIssues); L{end+1} = ['  - ' integrityIssues{i}]; end %#ok<SAGROW>
end
L{end+1} = repmat('=', 1, 70);
reportText = strjoin(L, '\n');
fprintf('%s\n', reportText);
fid = fopen(fullfile(REPORT_DIR, 'noise_audit_report.txt'), 'w');
fprintf(fid, '%s\n', reportText);
fclose(fid);
fprintf('\n=== X00 complete. Edit EXCLUSION_CONFIG_intervention.m by hand if needed, BEFORE running X05+. ===\n');

%% ======================= LOCAL FUNCTIONS =======================
function n = local_subject_numeric(subjVal)
    if iscategorical(subjVal) || isstring(subjVal); subjVal = cellstr(subjVal); end
    if iscell(subjVal)
        n = cellfun(@(s) str2double(regexprep(char(s), '\D', '')), subjVal);
    elseif isnumeric(subjVal)
        n = double(subjVal);
    else
        n = str2double(regexprep(char(subjVal), '\D', ''));
    end
end

function x = local_to_double(v)
    if isnumeric(v) || islogical(v)
        x = double(v);
    elseif iscategorical(v) || isstring(v)
        x = str2double(cellstr(v));
    elseif iscell(v)
        x = cellfun(@(e) local_one(e), v);
    else
        x = str2double(v);
    end
end

function y = local_one(e)
    if isnumeric(e) && isscalar(e); y = double(e);
    elseif isempty(e); y = NaN;
    else; y = str2double(char(e));
    end
end

function z = local_zscore(x)
    mu = mean(x, 'omitnan'); sigma = std(x, 0, 'omitnan');
    if sigma == 0 || isnan(sigma)
        z = zeros(size(x)); z(isnan(x)) = NaN;
    else
        z = (x - mu) / sigma;
    end
end

function s = local_numstr(v)
    if isnan(v); s = '  NA'; else; s = sprintf('%5.2f', v); end
end

function s = local_list(v)
    if isempty(v); s = 'none'; else; s = strjoin(cellstr(num2str(v(:))), ', '); end
end

function out = local_iif(cond, a, b)
    if cond; out = a; else; out = b; end
end
