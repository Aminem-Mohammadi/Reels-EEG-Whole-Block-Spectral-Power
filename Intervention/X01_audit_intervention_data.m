%% ========================================================================
% X01_audit_intervention_data.m
%
% Data inventory of the preprocessed intervention recordings (read-only).
% Finds every *.set file under DATA_DIR, parses subject / intervention number / condition
% from the file name (sub<N>_in<K>_<COND>_pp.set; COND i/I = Reels (E), N = neutral images,
% R = rest), loads each recording header with EEGLAB and checks sampling rate, channels,
% duration, events and condition-pair coverage.
% Inputs : DATA_DIR (*.set, Preprocessing_QC*.xlsx).
% Outputs: 01_audit/ (audit_table_intervention.csv and related tables), 09_reports/.
% ========================================================================

clear; clc;

%% ------------------------------------------------------------------

% Data folder: set DATA_DIR in config.m.
config;
DATA_ROOT = DATA_DIR;

% Reference expectations, used only to flag deviations in the report.
% These are NOT enforced/assumed anywhere in the processing below.
EXPECTED_SAMPLING_RATE = 512;   % Hz, as set in preprocessing
EXPECTED_NUM_CHANNELS  = 62;    % scalp channels after A1/A2 removal,
                                % per the 22-step preprocessing pipeline
EXPECTED_N_SUBJECTS    = 30;
EXPECTED_INTERVENTIONS_PER_SUBJECT = 2;   % Int1 and Int2 (design uses
                                           % only 2 of 3 conditions per
                                           % subject; confirmed by the
                                           % allocation table below)
DURATION_OUTLIER_FRACTION = 0.10;  % flag files whose duration deviates
                                    % from the group median by >10%

ALLOC_CROSS_VALIDATION_TRUSTED = false;

[DATA_ROOT_PARENT, ~] = fileparts(DATA_ROOT);
ANALYSIS_ROOT = OUT_DIR;
AUDIT_DIR  = fullfile(ANALYSIS_ROOT, '01_audit');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');

if ~exist(AUDIT_DIR, 'dir');  mkdir(AUDIT_DIR);  end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

fprintf('========================================================\n');
fprintf('PHASE 0 -- INTERVENTION DATA AUDIT (READ-ONLY)\n');
fprintf('========================================================\n');
fprintf('DATA_ROOT     : %s\n', DATA_ROOT);
fprintf('ANALYSIS_ROOT : %s\n', ANALYSIS_ROOT);
fprintf('--------------------------------------------------------\n\n');

if exist('pop_loadset', 'file') ~= 2
    error(['EEGLAB (pop_loadset) not found on the MATLAB path. ' ...
           'Start EEGLAB first, e.g.:  eeglab(''nogui'');  or  eeglab;']);
end

if ~exist(DATA_ROOT, 'dir')
    error('DATA_ROOT does not exist: %s', DATA_ROOT);
end

if ~ALLOC_CROSS_VALIDATION_TRUSTED
    fprintf(['NOTE: ALLOC_CROSS_VALIDATION_TRUSTED = false.\n' ...
             'Condition letters parsed from filenames will NOT be checked\n' ...
             'against the roadmap allocation table in Section 1, because the\n' ...
             'actual subject recording order is not yet confirmed to match\n' ...
             'that table''s row order. Filename-derived conditions are being\n' ...
             'reported as PROVISIONAL until a trustworthy mapping is supplied.\n\n']);
end

%% ------------------------------------------------------------------
%  1. PLANNED COUNTERBALANCING TABLE (used only for an optional cross-check)
%     Transcribed directly from "Experimental_Design_30_Subjects_3.pdf"
%     (columns: Subject, Gender, T1, Int1, T2, Int2, T3).
%
%     Subject IDs are stored as CHAR (not numeric double), matching the
%     digit string extracted from filenames by parse_filename() below,
%     so that IDs are compared as text.
% ------------------------------------------------------------------

alloc_raw = {
    '1',  'Female','SN','R','NS','N','SN';
    '2',  'Female','NS','R','SN','N','NS';
    '3',  'Female','NS','R','NS','N','SN';
    '4',  'Male',  'SN','R','SN','N','NS';
    '5',  'Male',  'SN','R','NS','N','NS';
    '6',  'Male',  'NS','N','SN','R','SN';
    '7',  'Male',  'SN','N','NS','R','NS';
    '8',  'Male',  'SN','N','SN','R','NS';
    '9',  'Female','NS','N','NS','R','SN';
    '10', 'Female','NS','N','SN','R','NS';
    '11', 'Female','SN','N','NS','E','SN';
    '12', 'Female','NS','N','SN','E','NS';
    '13', 'Female','SN','N','SN','E','NS';
    '14', 'Male',  'NS','N','NS','E','SN';
    '15', 'Male',  'SN','N','NS','E','NS';
    '16', 'Male',  'NS','E','SN','N','SN';
    '17', 'Male',  'SN','E','NS','N','NS';
    '18', 'Male',  'NS','E','NS','N','SN';
    '19', 'Female','SN','E','SN','N','NS';
    '20', 'Female','NS','E','SN','N','NS';
    '21', 'Female','SN','E','NS','R','SN';
    '22', 'Female','NS','E','SN','R','NS';
    '23', 'Female','SN','E','SN','R','NS';
    '24', 'Male',  'NS','E','NS','R','SN';
    '25', 'Male',  'SN','E','NS','R','NS';
    '26', 'Male',  'NS','R','SN','E','SN';
    '27', 'Male',  'SN','R','NS','E','NS';
    '28', 'Male',  'NS','R','NS','E','SN';
    '29', 'Female','SN','R','SN','E','NS';
    '30', 'Female','NS','R','SN','E','NS';
};

ALLOC = cell2table(alloc_raw, 'VariableNames', ...
    {'Subject','Gender','T1','Int1','T2','Int2','T3'});
ALLOC = normalize_text_columns(ALLOC);  % guarantees ALLOC.Subject{i} etc. always work

for i = 1:height(ALLOC)
    c1 = ALLOC.Int1{i}; c2 = ALLOC.Int2{i};
    if strcmp(c1, c2)
        warning('ALLOCATION TABLE: Subject %s has identical Int1/Int2 condition (%s). Check transcription.', ...
            ALLOC.Subject{i}, c1);
    end
    if ~all(ismember({c1, c2}, {'E','N','R'}))
        warning('ALLOCATION TABLE: Subject %s has an unrecognized condition code (%s / %s).', ...
            ALLOC.Subject{i}, c1, c2);
    end
end

%% ------------------------------------------------------------------
%  2. DISCOVER .set FILES (recursive)
% ------------------------------------------------------------------

set_list = dir(fullfile(DATA_ROOT, '**', '*.set'));
set_list = set_list(~[set_list.isdir]);

fprintf('Found %d .set file(s) under DATA_ROOT.\n\n', numel(set_list));

if isempty(set_list)
    error('No .set files found under DATA_ROOT. Check the path.');
end

% Also (informationally only) locate QC Excel files, one per subject
% folder, as mentioned in the project description. Not parsed here.
qc_list = dir(fullfile(DATA_ROOT, '**', 'Preprocessing_QC*.xlsx'));
qc_list = qc_list(~[qc_list.isdir]);
fprintf('Found %d Preprocessing_QC*.xlsx file(s) (logged only, not parsed in Phase 0).\n\n', numel(qc_list));

%% ------------------------------------------------------------------
%  3. FILENAME PARSING
% ------------------------------------------------------------------
%
% i.e.  sub<N>_in<K>_<COND>_pp.set   -- case of every letter is NOT
% guaranteed consistent (seen: N, R, i, r -- both cases occur even
% within the same dataset). The regex below is fully case-insensitive
% and tolerant of "in"/"int" and optional leading zeros.

% NOTE: the functions parse_filename(), map_condition_letter(), and
% get_expected_condition() that implement this parsing (including the
% exact regex pattern) are defined as LOCAL FUNCTIONS at the very end
% of this file -- MATLAB requires all local functions in a script to
% follow all top-level script code, so they cannot live here inline.

%% ------------------------------------------------------------------
%  4. MAIN LOOP -- LOAD EACH FILE, EXTRACT METADATA, CROSS-VALIDATE
% ------------------------------------------------------------------

n = numel(set_list);
rows = cell(n, 1);
chanlabel_store = cell(n, 1);
event_detail_store = cell(n, 1);

for i = 1:n

    fpath = fullfile(set_list(i).folder, set_list(i).name);
    fprintf('[%3d/%3d] %s ... ', i, n, set_list(i).name);

    row = struct();
    row.FileName = set_list(i).name;
    row.FolderPath = set_list(i).folder;
    row.FullPath = fpath;

    % ---- filename parsing ----
    [subj, intnum, condRaw, parseOK] = parse_filename(set_list(i).name);
    row.ParseOK = parseOK;
    row.Subject = subj;
    row.InterventionNumber = intnum;
    row.ConditionLetterRaw = condRaw;

    if parseOK
        [condCandidate, needsConfirm] = map_condition_letter(condRaw);
    else
        condCandidate = ''; needsConfirm = true;
    end
    row.ConditionCandidate = condCandidate;
    row.ConditionNeedsConfirmation = needsConfirm;

    % ---- cross-validate against allocation table (only if trusted) ----
    if ALLOC_CROSS_VALIDATION_TRUSTED && parseOK && ~isempty(subj) && ismember(intnum, [1 2])
        expectedCond = get_expected_condition(ALLOC, subj, intnum);
    else
        expectedCond = '';
    end
    row.ConditionExpectedFromAllocation = expectedCond;

    if ~ALLOC_CROSS_VALIDATION_TRUSTED
        % Cross-check disabled (see Section 0): no MATCH/MISMATCH is reported.
        row.ConditionMatch = 'NOT_CROSS_CHECKED';
    elseif isempty(expectedCond) || isempty(condCandidate)
        row.ConditionMatch = 'N/A';
    elseif strcmp(condCandidate, expectedCond)
        row.ConditionMatch = 'MATCH';
    else
        row.ConditionMatch = 'MISMATCH';
    end

    % ---- load the dataset (metadata + data, read-only) ----
    notes = {};
    try
        EEG = pop_loadset('filename', set_list(i).name, ...
                           'filepath', set_list(i).folder, ...
                           'loadmode', 'all');

        row.LoadStatus = 'OK';
        row.LoadErrorMessage = '';

        row.SamplingRateHz = EEG.srate;
        row.NumChannels    = EEG.nbchan;
        row.NumSamples     = EEG.pnts;
        row.DurationSec    = EEG.pnts / EEG.srate;
        row.DurationMin    = row.DurationSec / 60;
        row.SetName        = EEG.setname;

        % channel labels
        if isfield(EEG, 'chanlocs') && ~isempty(EEG.chanlocs) && isfield(EEG.chanlocs, 'labels')
            labels = {EEG.chanlocs.labels};
        else
            labels = {};
        end
        chanlabel_store{i} = labels;
        row.NumChannelLabels = numel(labels);

        hasXYZ = false;
        if ~isempty(EEG.chanlocs) && isfield(EEG.chanlocs, 'X') && ~isempty(EEG.chanlocs(1).X)
            hasXYZ = true;
        end
        row.ChanlocsHaveXYZ = hasXYZ;

        % events
        if isfield(EEG, 'event') && ~isempty(EEG.event)
            row.NumEvents = numel(EEG.event);
            if isfield(EEG.event, 'type')
                etypes = {EEG.event.type};
                etypes = cellfun(@event_type_to_str, etypes, 'UniformOutput', false);
                [u, ~, ic] = unique(etypes);
                cnts = accumarray(ic, 1);
                summary_parts = arrayfun(@(k) sprintf('%s:%d', u{k}, cnts(k)), 1:numel(u), 'UniformOutput', false);
                row.EventTypesSummary = strjoin(summary_parts, '; ');
            else
                row.EventTypesSummary = '';
            end
        else
            row.NumEvents = 0;
            row.EventTypesSummary = '';
        end
        event_detail_store{i} = row.EventTypesSummary;

        % ICA
        row.HasICAWeights = ~isempty(EEG.icaweights);
        if row.HasICAWeights
            row.NumICs = size(EEG.icaweights, 1);
        else
            row.NumICs = 0;
        end

        % NaN / Inf
        if ~isempty(EEG.data)
            row.NaNCount = sum(isnan(EEG.data(:)));
            row.InfCount = sum(isinf(EEG.data(:)));
        else
            row.NaNCount = NaN;
            row.InfCount = NaN;
            notes{end+1} = 'EMPTY_DATA';
        end

        % ---- flags vs. reference expectations ----
        if EEG.srate ~= EXPECTED_SAMPLING_RATE
            notes{end+1} = sprintf('SRATE_UNEXPECTED(%g)', EEG.srate);
        end
        if EEG.nbchan ~= EXPECTED_NUM_CHANNELS
            notes{end+1} = sprintf('NCHAN_UNEXPECTED(%d)', EEG.nbchan);
        end
        if row.NaNCount > 0
            notes{end+1} = sprintf('CONTAINS_NAN(%d)', row.NaNCount);
        end
        if row.InfCount > 0
            notes{end+1} = sprintf('CONTAINS_INF(%d)', row.InfCount);
        end
        if strcmp(row.ConditionMatch, 'MISMATCH')
            notes{end+1} = sprintf('CONDITION_MISMATCH(file=%s,expected=%s)', condCandidate, expectedCond);
        end
        if ~parseOK
            notes{end+1} = 'UNPARSEABLE_FILENAME';
        end
        if needsConfirm && parseOK
            notes{end+1} = sprintf('CONDITION_LETTER_NEEDS_CONFIRMATION(raw=%s)', condRaw);
        end

        fprintf('OK  (Fs=%g, chans=%d, dur=%.1f s, cond[file/alloc]=%s/%s)\n', ...
            EEG.srate, EEG.nbchan, row.DurationSec, condCandidate, expectedCond);

    catch ME
        row.LoadStatus = 'LOAD_ERROR';
        row.LoadErrorMessage = ME.message;
        row.SamplingRateHz = NaN;
        row.NumChannels = NaN;
        row.NumSamples = NaN;
        row.DurationSec = NaN;
        row.DurationMin = NaN;
        row.SetName = '';
        row.NumChannelLabels = NaN;
        row.ChanlocsHaveXYZ = false;
        row.NumEvents = NaN;
        row.EventTypesSummary = '';
        row.HasICAWeights = false;
        row.NumICs = NaN;
        row.NaNCount = NaN;
        row.InfCount = NaN;
        notes{end+1} = 'LOAD_ERROR';
        chanlabel_store{i} = {};
        event_detail_store{i} = '';
        fprintf('LOAD_ERROR: %s\n', ME.message);
    end

    row.Notes = strjoin(notes, ' | ');
    if isempty(notes) || (numel(notes) == 1 && strcmp(notes{1}, ''))
        row.Status = 'OK';
    elseif any(strcmp(notes, 'LOAD_ERROR')) || any(strcmp(notes, 'UNPARSEABLE_FILENAME'))
        row.Status = 'ERROR';
    else
        row.Status = 'REVIEW_NEEDED';
    end

    rows{i} = row;
end

%% ------------------------------------------------------------------

durations = cellfun(@(r) r.DurationSec, rows);
valid_dur = durations(~isnan(durations));
if ~isempty(valid_dur)
    med_dur = median(valid_dur);
    for i = 1:n
        d = rows{i}.DurationSec;
        if ~isnan(d) && med_dur > 0 && abs(d - med_dur) / med_dur > DURATION_OUTLIER_FRACTION
            rows{i}.Notes = append_note(rows{i}.Notes, sprintf('DURATION_OUTLIER(%.1fs vs median %.1fs)', d, med_dur));
            if strcmp(rows{i}.Status, 'OK')
                rows{i}.Status = 'REVIEW_NEEDED';
            end
        end
    end
else
    med_dur = NaN;
end

%% ------------------------------------------------------------------
%  6. DUPLICATE DETECTION (same Subject + InterventionNumber, >1 file)
% ------------------------------------------------------------------

subj_all = cellfun(@(r) r.Subject, rows, 'UniformOutput', false);  % cell array of char
int_all  = cellfun(@(r) r.InterventionNumber, rows);               % numeric (1/2/NaN)

dup_rows = cell(0, 4);  % pre-sized (0 rows, 4 cols) so cell2table works even if empty
valid_idx = find(~cellfun(@isempty, subj_all) & ~isnan(int_all));
combo_keys = arrayfun(@(i) sprintf('%s_%d', subj_all{i}, int_all(i)), valid_idx, 'UniformOutput', false);
[u_keys, ~, ic] = unique(combo_keys);
for k = 1:numel(u_keys)
    members = valid_idx(ic == k);
    if numel(members) > 1
        for m = 1:numel(members)
            r = rows{members(m)};
            dup_rows(end+1, :) = {r.Subject, r.InterventionNumber, r.FileName, r.FullPath}; %#ok<AGROW>
            rows{members(m)}.Notes = append_note(rows{members(m)}.Notes, 'DUPLICATE_RECORDING');
            rows{members(m)}.Status = 'REVIEW_NEEDED';
        end
    end
end

%% ------------------------------------------------------------------

missing_rows = cell(0, 3);  % pre-sized so cell2table works even if empty
if ALLOC_CROSS_VALIDATION_TRUSTED
    for s = 1:height(ALLOC)
        subjID = ALLOC.Subject{s};
        for intnum = 1:2
            expectedCond = get_expected_condition(ALLOC, subjID, intnum);
            found = any(strcmp(subj_all, subjID) & int_all == intnum);
            if ~found
                missing_rows(end+1, :) = {subjID, intnum, expectedCond}; %#ok<AGROW>
            end
        end
    end
end

%% ------------------------------------------------------------------
%  8. CONDITION MISMATCH SUB-TABLE
% ------------------------------------------------------------------

mismatch_rows = cell(0, 6);  % pre-sized so cell2table works even if empty
for i = 1:n
    if strcmp(rows{i}.ConditionMatch, 'MISMATCH')
        r = rows{i};
        mismatch_rows(end+1, :) = {r.Subject, r.InterventionNumber, r.FileName, ...
            r.ConditionLetterRaw, r.ConditionCandidate, r.ConditionExpectedFromAllocation}; %#ok<AGROW>
    end
end

%% ------------------------------------------------------------------
%  9. CHANNEL LABEL CONSISTENCY CHECK
% ------------------------------------------------------------------

ok_idx = find(strcmp(cellfun(@(r) r.LoadStatus, rows, 'UniformOutput', false), 'OK'));
label_mismatch_files = {};
ref_labels = {};
if ~isempty(ok_idx)
    ref_labels = chanlabel_store{ok_idx(1)};
    for k = 1:numel(ok_idx)
        i = ok_idx(k);
        lbl = chanlabel_store{i};
        if ~isequal(sort(lbl), sort(ref_labels))
            label_mismatch_files{end+1} = rows{i}.FileName; %#ok<AGROW>
            rows{i}.Notes = append_note(rows{i}.Notes, 'CHANNEL_LABEL_SET_DIFFERS_FROM_REFERENCE');
            if strcmp(rows{i}.Status, 'OK')
                rows{i}.Status = 'REVIEW_NEEDED';
            end
        end
    end
end

%% ------------------------------------------------------------------
%  10. ASSEMBLE MASTER TABLE AND SAVE
% ------------------------------------------------------------------

T = struct2table([rows{:}]);

% Reorder columns for readability
preferredOrder = {'Subject','InterventionNumber','ConditionLetterRaw', ...
    'ConditionCandidate','ConditionNeedsConfirmation', ...
    'ConditionExpectedFromAllocation','ConditionMatch', ...
    'FileName','FolderPath','LoadStatus','LoadErrorMessage', ...
    'DurationSec','DurationMin','SamplingRateHz','NumChannels', ...
    'NumSamples','NumChannelLabels','ChanlocsHaveXYZ','NumEvents', ...
    'EventTypesSummary','HasICAWeights','NumICs','NaNCount','InfCount', ...
    'SetName','ParseOK','FullPath','Notes','Status'};
existingCols = preferredOrder(ismember(preferredOrder, T.Properties.VariableNames));
T = T(:, existingCols);

% Safety: struct2table/cell2table can silently store a text column as a
% char MATRIX (instead of a cell array of char) whenever every value in
% that column happens to have the same length in a given run (e.g. if
% every file loads OK, LoadStatus would otherwise become a char matrix).
% That breaks the T.Col{i}-style indexing used throughout this script.
% Force every text column to a cellstr so indexing is always uniform,
% regardless of what values happen to occur in any particular run.
T = normalize_text_columns(T);

writetable(T, fullfile(AUDIT_DIR, 'audit_table_intervention.csv'));

MissingT = cell2table(missing_rows, 'VariableNames', {'Subject','InterventionNumber','ExpectedCondition'});
MissingT = normalize_text_columns(MissingT);
writetable(MissingT, fullfile(AUDIT_DIR, 'missing_recordings.csv'));

DupT = cell2table(dup_rows, 'VariableNames', {'Subject','InterventionNumber','FileName','FullPath'});
DupT = normalize_text_columns(DupT);
writetable(DupT, fullfile(AUDIT_DIR, 'duplicate_recordings.csv'));

MismatchT = cell2table(mismatch_rows, 'VariableNames', ...
    {'Subject','InterventionNumber','FileName','ConditionLetterRaw','ConditionCandidate','ConditionExpectedFromAllocation'});
MismatchT = normalize_text_columns(MismatchT);
writetable(MismatchT, fullfile(AUDIT_DIR, 'condition_mismatches.csv'));

save(fullfile(AUDIT_DIR, 'audit_table_intervention.mat'), ...
    'T', 'MissingT', 'DupT', 'MismatchT', 'ALLOC', 'chanlabel_store', ...
    'ref_labels', 'label_mismatch_files', 'med_dur');

%% ------------------------------------------------------------------
%  11. HUMAN-READABLE REPORT
% ------------------------------------------------------------------

rpt = fullfile(REPORT_DIR, 'phase0_audit_report_intervention.txt');
fid = fopen(rpt, 'w');

fprintf(fid, '========================================================\n');
fprintf(fid, 'PHASE 0 -- INTERVENTION DATA AUDIT REPORT\n');
fprintf(fid, 'Generated: %s\n', datestr(now));
fprintf(fid, '========================================================\n\n');

fprintf(fid, 'DATA_ROOT: %s\n\n', DATA_ROOT);

fprintf(fid, '--- ALLOCATION-TABLE CROSS-VALIDATION STATUS ---\n');
if ALLOC_CROSS_VALIDATION_TRUSTED
    fprintf(fid, '  ENABLED. Filename-parsed conditions were checked against the\n');
    fprintf(fid, '  hard-coded roadmap table and any disagreement is listed below.\n');
else
    fprintf(fid, '  DISABLED. The roadmap allocation table (Section 1 of the script)\n');
    fprintf(fid, '  was the PLANNED counterbalancing order, but actual data collection\n');
    fprintf(fid, '  may not have followed that same subject order, so "Subject N" in\n');
    fprintf(fid, '  the table is not yet confirmed to be the same participant as\n');
    fprintf(fid, '  "subN" in the filenames. Every condition below is therefore the\n');
    fprintf(fid, '  PROVISIONAL, filename-derived value only -- NOT cross-checked\n');
    fprintf(fid, '  against any independent ground truth. "Missing expected recordings"\n');
    fprintf(fid, '  and "Condition mismatches" cannot be computed while this is\n');
    fprintf(fid, '  disabled, and are not reported as if they were.\n');
end
fprintf(fid, '\n');

fprintf(fid, '--- OVERVIEW ---\n');
fprintf(fid, 'Total .set files found            : %d\n', n);
fprintf(fid, 'Expected subjects (per allocation) : %d\n', EXPECTED_N_SUBJECTS);
fprintf(fid, 'Expected recordings/subject        : %d\n', EXPECTED_INTERVENTIONS_PER_SUBJECT);
fprintf(fid, 'Expected total recordings          : %d\n', EXPECTED_N_SUBJECTS * EXPECTED_INTERVENTIONS_PER_SUBJECT);
fprintf(fid, 'Files successfully loaded          : %d\n', sum(strcmp(T.LoadStatus, 'OK')));
fprintf(fid, 'Files with LOAD_ERROR              : %d\n', sum(strcmp(T.LoadStatus, 'LOAD_ERROR')));
fprintf(fid, 'Files with unparseable filenames   : %d\n', sum(~T.ParseOK));
if ALLOC_CROSS_VALIDATION_TRUSTED
    fprintf(fid, 'Missing expected recordings        : %d\n', height(MissingT));
    fprintf(fid, 'Condition mismatches (file vs alloc): %d\n', height(MismatchT));
else
    fprintf(fid, 'Missing expected recordings        : N/A (cross-validation disabled)\n');
    fprintf(fid, 'Condition mismatches (file vs alloc): N/A (cross-validation disabled)\n');
end
fprintf(fid, 'Duplicate recordings (file rows)   : %d\n', height(DupT));
fprintf(fid, 'Files flagged CONTAINS_NAN         : %d\n', sum(contains(T.Notes, 'CONTAINS_NAN')));
fprintf(fid, 'Files flagged CONTAINS_INF         : %d\n', sum(contains(T.Notes, 'CONTAINS_INF')));
fprintf(fid, 'Files with channel-label-set mismatch vs reference: %d\n', numel(label_mismatch_files));
fprintf(fid, '\n');

fprintf(fid, '--- SAMPLING RATE ---\n');
srates = T.SamplingRateHz(~isnan(T.SamplingRateHz));
if ~isempty(srates)
    u_sr = unique(srates);
    for k = 1:numel(u_sr)
        fprintf(fid, '  Fs = %g Hz : %d file(s)%s\n', u_sr(k), sum(srates == u_sr(k)), ...
            iif(u_sr(k) ~= EXPECTED_SAMPLING_RATE, '   <-- UNEXPECTED', ''));
    end
end
fprintf(fid, '\n');

fprintf(fid, '--- CHANNEL COUNT ---\n');
nchans = T.NumChannels(~isnan(T.NumChannels));
if ~isempty(nchans)
    u_nc = unique(nchans);
    for k = 1:numel(u_nc)
        fprintf(fid, '  %d channels : %d file(s)%s\n', u_nc(k), sum(nchans == u_nc(k)), ...
            iif(u_nc(k) ~= EXPECTED_NUM_CHANNELS, '   <-- UNEXPECTED', ''));
    end
end
fprintf(fid, '\n');

fprintf(fid, '--- DURATION ---\n');
if ~isempty(valid_dur)
    fprintf(fid, '  Median duration : %.2f s (%.2f min)\n', med_dur, med_dur/60);
    fprintf(fid, '  Min duration    : %.2f s\n', min(valid_dur));
    fprintf(fid, '  Max duration    : %.2f s\n', max(valid_dur));
    fprintf(fid, '  Files flagged as duration outliers (>%.0f%% from median): %d\n', ...
        DURATION_OUTLIER_FRACTION*100, sum(contains(T.Notes, 'DURATION_OUTLIER')));
end
fprintf(fid, '\n');

if ALLOC_CROSS_VALIDATION_TRUSTED
    fprintf(fid, '--- CONDITION COUNTS (per ConditionExpectedFromAllocation, authoritative) ---\n');
    condsFound = T.ConditionExpectedFromAllocation(strcmp(T.LoadStatus,'OK'));
else
    fprintf(fid, '--- CONDITION COUNTS (per ConditionCandidate, PROVISIONAL / filename-derived only) ---\n');
    condsFound = T.ConditionCandidate(strcmp(T.LoadStatus,'OK'));
end
for c = {'E','N','R'}
    fprintf(fid, '  %s : %d recording(s)\n', c{1}, sum(strcmp(condsFound, c{1})));
end
fprintf(fid, '\n');

fprintf(fid, '--- INTERVENTION POSITION COUNTS ---\n');
for k = 1:2
    fprintf(fid, '  Int%d : %d recording(s)\n', k, sum(T.InterventionNumber == k));
end
fprintf(fid, '\n');

fprintf(fid, '--- MISSING EXPECTED RECORDINGS ---\n');
if ~ALLOC_CROSS_VALIDATION_TRUSTED
    fprintf(fid, '  N/A -- cross-validation disabled (see status note above).\n');
elseif isempty(MissingT)
    fprintf(fid, '  None.\n');
else
    for i = 1:height(MissingT)
        fprintf(fid, '  Sub%s, Int%d (expected condition: %s)\n', ...
            MissingT.Subject{i}, MissingT.InterventionNumber(i), MissingT.ExpectedCondition{i});
    end
end
fprintf(fid, '\n');

fprintf(fid, '--- DUPLICATE RECORDINGS ---\n');
if isempty(DupT)
    fprintf(fid, '  None.\n');
else
    for i = 1:height(DupT)
        fprintf(fid, '  Sub%s, Int%d : %s\n', DupT.Subject{i}, DupT.InterventionNumber(i), DupT.FileName{i});
    end
end
fprintf(fid, '\n');

fprintf(fid, '--- CONDITION MISMATCHES (filename-derived vs. allocation-table) ---\n');
if ~ALLOC_CROSS_VALIDATION_TRUSTED
    fprintf(fid, '  N/A -- cross-validation disabled (see status note above).\n');
elseif isempty(MismatchT)
    fprintf(fid, '  None. All parsed condition letters agree with the allocation table.\n');
else
    fprintf(fid, '  ACTION NEEDED -- do not resolve automatically:\n');
    for i = 1:height(MismatchT)
        fprintf(fid, '  Sub%s, Int%d, file=%s : raw="%s" -> candidate=%s, but allocation table expects %s\n', ...
            MismatchT.Subject{i}, MismatchT.InterventionNumber(i), MismatchT.FileName{i}, ...
            MismatchT.ConditionLetterRaw{i}, MismatchT.ConditionCandidate{i}, MismatchT.ConditionExpectedFromAllocation{i});
    end
end
fprintf(fid, '\n');

fprintf(fid, '--- FILES NEEDING CONDITION-LETTER CONFIRMATION (e.g. raw "I"/"i") ---\n');
confirmIdx = find(contains(T.Notes, 'CONDITION_LETTER_NEEDS_CONFIRMATION'));
if isempty(confirmIdx)
    fprintf(fid, '  None.\n');
else
    for i = confirmIdx'
        fprintf(fid, '  %s : raw="%s" -> assumed %s (matches allocation: %s)\n', ...
            T.FileName{i}, T.ConditionLetterRaw{i}, T.ConditionCandidate{i}, T.ConditionMatch{i});
    end
end
fprintf(fid, '\n');

fprintf(fid, '--- UNPARSEABLE FILENAMES ---\n');
badIdx = find(~T.ParseOK);
if isempty(badIdx)
    fprintf(fid, '  None.\n');
else
    for i = badIdx'
        fprintf(fid, '  %s\n', T.FileName{i});
    end
end
fprintf(fid, '\n');

fprintf(fid, '--- LOAD ERRORS ---\n');
errIdx = find(strcmp(T.LoadStatus, 'LOAD_ERROR'));
if isempty(errIdx)
    fprintf(fid, '  None.\n');
else
    for i = errIdx'
        fprintf(fid, '  %s : %s\n', T.FileName{i}, T.LoadErrorMessage{i});
    end
end
fprintf(fid, '\n');

fprintf(fid, '--- CHANNEL LABEL SET MISMATCHES (vs. first successfully loaded file) ---\n');
if isempty(ok_idx)
    fprintf(fid, '  N/A -- no file loaded successfully, so no reference labels are available.\n');
elseif isempty(label_mismatch_files)
    fprintf(fid, '  Reference file: %s\n', T.FileName{ok_idx(1)});
    fprintf(fid, '  None -- all files share the same channel-label set.\n');
else
    fprintf(fid, '  Reference file: %s\n', T.FileName{ok_idx(1)});
    for i = 1:numel(label_mismatch_files)
        fprintf(fid, '  %s\n', label_mismatch_files{i});
    end
end
fprintf(fid, '\n');

fprintf(fid, '--- ROWS FLAGGED REVIEW_NEEDED OR ERROR ---\n');
flagIdx = find(~strcmp(T.Status, 'OK'));
if isempty(flagIdx)
    fprintf(fid, '  None -- every file passed all checks.\n');
else
    for i = flagIdx'
        fprintf(fid, '  [%s] %s : %s\n', T.Status{i}, T.FileName{i}, T.Notes{i});
    end
end
fprintf(fid, '\n');

fprintf(fid, '--- QC EXCEL FILES FOUND (logged only, not parsed) ---\n');
fprintf(fid, '  %d Preprocessing_QC*.xlsx file(s) found.\n', numel(qc_list));
for i = 1:numel(qc_list)
    fprintf(fid, '  %s\n', fullfile(qc_list(i).folder, qc_list(i).name));
end

fclose(fid);

fprintf('\n========================================================\n');
fprintf('AUDIT COMPLETE.\n');
fprintf('Master table : %s\n', fullfile(AUDIT_DIR, 'audit_table_intervention.csv'));
fprintf('Report       : %s\n', rpt);
fprintf('========================================================\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%  (must appear here, after all top-level script code, per MATLAB rules
%  for scripts that mix executable statements with local functions)
% ========================================================================

function [subj, intnum, condRaw, ok] = parse_filename(fname)
    % Expected convention: sub<N>_in<K>_<COND>_pp.set
    % Fully case-insensitive; tolerant of "in"/"int" and leading zeros.
    %
    % subj is returned as CHAR (e.g. '16'), not a numeric double, so it
    % is compared as text with the allocation table (Section 1).
    FNAME_PATTERN = '^[Ss]ub0*(\d+)_[Ii][Nn][Tt]?0*(\d+)_([A-Za-z]+)_pp$';
    [~, base, ~] = fileparts(fname);
    tok = regexp(base, FNAME_PATTERN, 'tokens', 'once');
    if isempty(tok)
        subj = ''; intnum = NaN; condRaw = ''; ok = false;
    else
        subj    = tok{1};              % CHAR, e.g. '16'
        intnum  = str2double(tok{2});  % numeric is fine here (only 1 or 2)
        condRaw = tok{3};
        ok = true;
    end
end

function [condCandidate, needsConfirm] = map_condition_letter(condRaw)
    c = upper(condRaw);
    switch c
        case 'E'
            condCandidate = 'E'; needsConfirm = false;
        case 'I'
            condCandidate = 'E'; needsConfirm = true;
        case 'N'
            condCandidate = 'N'; needsConfirm = false;
        case 'R'
            condCandidate = 'R'; needsConfirm = false;
        otherwise
            condCandidate = ''; needsConfirm = true;
    end
end

function expected = get_expected_condition(ALLOC, subj, intnum)
    % subj is CHAR (e.g. '16'); matched by string comparison, not
    % numeric equality, against ALLOC.Subject (also char/cellstr).
    expected = '';
    if isempty(subj)
        return;
    end
    idx = find(strcmp(ALLOC.Subject, subj), 1);
    if isempty(idx)
        return;
    end
    if intnum == 1
        expected = ALLOC.Int1{idx};
    elseif intnum == 2
        expected = ALLOC.Int2{idx};
    end
end

function out = iif(cond, a, b)
    % Inline conditional string, used only for report text.
    if cond
        out = a;
    else
        out = b;
    end
end

function Tout = normalize_text_columns(Tin)
    % Ensures every text-like table variable is a cellstr (not a plain
    % char matrix), so T.Col{i}-style curly-brace indexing is always
    % valid downstream, regardless of whether struct2table/cell2table
    % happened to collapse a same-length text column into a char matrix.
    Tout = Tin;
    vn = Tout.Properties.VariableNames;
    for k = 1:numel(vn)
        col = Tout.(vn{k});
        if ischar(col)
            Tout.(vn{k}) = cellstr(col);
        end
    end
end

function s = event_type_to_str(x)
    % EEG.event.type may be numeric or char depending on EEGLAB version /
    % how the artificial marker channel was converted -- handle both
    % rather than assuming, since num2str() errors on char input.
    if ischar(x)
        s = x;
    elseif isstring(x)
        s = char(x);
    elseif isnumeric(x) || islogical(x)
        s = num2str(x);
    else
        s = 'UNKNOWN_TYPE';
    end
end

function out = append_note(existing, newNote)
    if isempty(strtrim(existing))
        out = newNote;
    else
        out = [strtrim(existing) ' | ' newNote];
    end
end
