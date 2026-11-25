%% ========================================================================
% X02_build_intervention_metadata.m
%
% Builds the recording-level metadata table (one row per recording: Subject,
% InterventionNumber, Condition, file paths, duration, events, boundary information) and
% attaches the preprocessing QC values (bad channels, bad ICs) read from
% Preprocessing_QC_<n>.xlsx.
% Inputs : DATA_DIR.
% Outputs: 02_metadata/intervention_metadata.mat/.xlsx, qc_excel_inspection.csv, 09_reports/.
% ========================================================================

clear; clc;

%% ------------------------------------------------------------------

% Data folder: set DATA_DIR in config.m.
config;
DATA_ROOT = DATA_DIR;

% Reference expectations, used ONLY to flag deviations for review.
EXPECTED_SAMPLING_RATE = 512;
EXPECTED_NUM_CHANNELS  = 62;
EXPECTED_N_SUBJECTS    = 30;
EXPECTED_INTERVENTIONS_PER_SUBJECT = 2;
DURATION_OUTLIER_FRACTION = 0.10;
NOMINAL_DURATION_SEC = 600;  % nominal 10-min intervention; reference only,
                             % NOT used to exclude/alter any recording

[DATA_ROOT_PARENT, ~] = fileparts(DATA_ROOT);
ANALYSIS_ROOT = OUT_DIR;
META_DIR   = fullfile(ANALYSIS_ROOT, '02_metadata');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');
if ~exist(META_DIR, 'dir');   mkdir(META_DIR);   end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

fprintf('========================================================\n');
fprintf('PHASE 1 -- BUILD INTERVENTION METADATA (READ-ONLY)\n');
fprintf('========================================================\n');
fprintf('DATA_ROOT   : %s\n', DATA_ROOT);
fprintf('OUTPUT      : %s\n', META_DIR);
fprintf('--------------------------------------------------------\n\n');

if exist('pop_loadset', 'file') ~= 2
    error(['EEGLAB (pop_loadset) not found on the MATLAB path. ' ...
           'Start EEGLAB first, e.g.:  eeglab(''nogui'');  or  eeglab;']);
end
if ~exist(DATA_ROOT, 'dir')
    error('DATA_ROOT does not exist: %s', DATA_ROOT);
end

%% ------------------------------------------------------------------

QC_FIELD_MAP = struct( ...
    'Subject_raw',         {{'Subject'}}, ...
    'Task_raw',            {{'Task'}}, ...
    'TaskNumber_raw',      {{'TaskNumber'}}, ...
    'Dataset_raw',         {{'Dataset'}}, ...
    'InitialChannels',     {{'InitialChannels'}}, ...
    'BadChannels',         {{'BadChannels'}}, ...
    'NBadChannels',        {{'NBadChannels'}}, ...
    'NBadICs',             {{'NBadICs'}}, ...
    'Brain',               {{'Brain'}}, ...
    'Muscle',              {{'Muscle'}}, ...
    'Eye',                 {{'Eye'}}, ...
    'Heart',               {{'Heart'}}, ...
    'LineNoise',           {{'LineNoise'}}, ...
    'ChannelNoise',        {{'ChannelNoise'}}, ...
    'Other',               {{'Other'}}, ...
    'FinalChannels',       {{'FinalChannels'}}, ...
    'Duration',            {{'Duration', 'DurationSec', 'RecordingDuration'}}, ...
    'SamplingRate',        {{'SamplingRate', 'Fs', 'SampleRate'}}, ...
    'PreprocessingStatus', {{'Status', 'PreprocessingStatus'}} ...
);
QC_CANON_FIELDS = fieldnames(QC_FIELD_MAP);

%% ------------------------------------------------------------------
%  2. DISCOVER .set FILES (recursive)
% ------------------------------------------------------------------

set_list = dir(fullfile(DATA_ROOT, '**', '*.set'));
set_list = set_list(~[set_list.isdir]);
fprintf('Found %d .set file(s) under DATA_ROOT.\n\n', numel(set_list));
if isempty(set_list)
    error('No .set files found under DATA_ROOT. Check the path.');
end

%% ------------------------------------------------------------------
%  3. MAIN LOOP -- PARSE, LOAD METADATA, LOCATE/INSPECT QC, VALIDATE
% ------------------------------------------------------------------

n = numel(set_list);
rows = cell(n, 1);
chanlabel_store = cell(n, 1);
qc_inspect_log = cell(0, 10);  % pre-sized (0 rows, 10 cols) so cell2table works even if empty

for i = 1:n

    fpath = fullfile(set_list(i).folder, set_list(i).name);
    fprintf('[%3d/%3d] %s ... ', i, n, set_list(i).name);

    row = struct();
    row.FileName   = set_list(i).name;
    row.FolderPath = set_list(i).folder;
    row.FullPath   = fpath;

    % ---- filename parsing ----
    [subj, intnum, condRaw, parseOK] = parse_filename(set_list(i).name);
    row.ParseOK = parseOK;
    row.Subject = subj;                       % char, e.g. '16'
    row.InterventionNumber = intnum;          % numeric 1/2/NaN
    row.ConditionLetterRaw = condRaw;          % preserved as-is
    row.Condition = normalize_condition(condRaw);   % 'E'/'N'/'R'/''

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

        row.HasICAWeights = ~isempty(EEG.icaweights);
        if row.HasICAWeights
            row.NumICs = size(EEG.icaweights, 1);
        else
            row.NumICs = 0;
        end

        if ~isempty(EEG.data)
            row.NaNCount = sum(isnan(EEG.data(:)));
            row.InfCount = sum(isinf(EEG.data(:)));
        else
            row.NaNCount = NaN;
            row.InfCount = NaN;
            notes{end+1} = 'EMPTY_DATA'; %#ok<AGROW>
        end

        [s_subj, s_int, s_cond, s_ok] = parse_setname(EEG.setname);
        row.SetNameParsed = s_ok;
        mismatchParts = {};
        if s_ok
            if ~isempty(s_subj) && ~isempty(subj) && ~strcmp(s_subj, subj)
                mismatchParts{end+1} = sprintf('subject(file=%s,setname=%s)', subj, s_subj); %#ok<AGROW>
            end
            if ~isnan(s_int) && ~isnan(intnum) && s_int ~= intnum
                mismatchParts{end+1} = sprintf('intervention(file=%d,setname=%g)', intnum, s_int); %#ok<AGROW>
            end
            if ~isempty(s_cond)
                s_cond_norm = normalize_condition(s_cond);
                if ~isempty(s_cond_norm) && ~isempty(row.Condition) && ~strcmp(s_cond_norm, row.Condition)
                    mismatchParts{end+1} = sprintf('condition(file=%s,setname=%s)', row.Condition, s_cond_norm); %#ok<AGROW>
                end
            end
        end
        row.SetNameFilenameMismatch = ~isempty(mismatchParts);
        row.SetNameMismatchDetail = strjoin(mismatchParts, '; ');
        if row.SetNameFilenameMismatch
            notes{end+1} = ['SETNAME_MISMATCH(' row.SetNameMismatchDetail ')']; %#ok<AGROW>
        end

        if EEG.srate ~= EXPECTED_SAMPLING_RATE
            notes{end+1} = sprintf('SRATE_UNEXPECTED(%g)', EEG.srate); %#ok<AGROW>
        end
        if EEG.nbchan ~= EXPECTED_NUM_CHANNELS
            notes{end+1} = sprintf('NCHAN_UNEXPECTED(%d)', EEG.nbchan); %#ok<AGROW>
        end
        if row.NaNCount > 0
            notes{end+1} = sprintf('CONTAINS_NAN(%d)', row.NaNCount); %#ok<AGROW>
        end
        if row.InfCount > 0
            notes{end+1} = sprintf('CONTAINS_INF(%d)', row.InfCount); %#ok<AGROW>
        end

        fprintf('OK  (Fs=%g, chans=%d, dur=%.1f s, cond=%s)\n', ...
            EEG.srate, EEG.nbchan, row.DurationSec, row.Condition);

    catch ME
        row.LoadStatus = 'LOAD_ERROR';
        row.LoadErrorMessage = ME.message;
        row.SamplingRateHz = NaN; row.NumChannels = NaN; row.NumSamples = NaN;
        row.DurationSec = NaN; row.DurationMin = NaN; row.SetName = '';
        row.NumChannelLabels = NaN; row.ChanlocsHaveXYZ = false;
        row.NumEvents = NaN; row.EventTypesSummary = '';
        row.HasICAWeights = false; row.NumICs = NaN;
        row.NaNCount = NaN; row.InfCount = NaN;
        row.SetNameParsed = false; row.SetNameFilenameMismatch = false;
        row.SetNameMismatchDetail = '';
        notes{end+1} = 'LOAD_ERROR'; %#ok<AGROW>
        chanlabel_store{i} = {};
        fprintf('LOAD_ERROR: %s\n', ME.message);
    end

    % ---- locate + inspect the corresponding QC Excel file ----
    [qcPath, qcFileName, qcHowFound] = find_qc_file(set_list(i).folder, intnum);
    if isempty(qcPath)
        row.QCFileExists = false;
        row.QCFileName = '';
        row.QCFilePath = '';
        qcResult = empty_qc_result(QC_CANON_FIELDS);
        qcResult.Notes = 'No matching Preprocessing_QC file found for this intervention.';
    else
        row.QCFileExists = true;
        row.QCFileName = qcFileName;
        row.QCFilePath = qcPath;
        try
            qcResult = inspect_qc_excel(qcPath, subj, intnum, QC_FIELD_MAP);
        catch ME2
            qcResult = empty_qc_result(QC_CANON_FIELDS);
            qcResult.Notes = ['Unhandled error inspecting QC file: ' ME2.message];
        end
        if ~strcmp(qcHowFound, 'exact')
            qcResult.Notes = strtrim([qcResult.Notes ' | Matched by fuzzy filename search, not the exact expected name.']);
        end
    end

    % flatten QC canonical fields onto the row with a QC_ prefix
    for k = 1:numel(QC_CANON_FIELDS)
        row.(['QC_' QC_CANON_FIELDS{k}]) = qcResult.Mapped.(QC_CANON_FIELDS{k});
    end
    row.QCSheetName        = qcResult.ChosenSheet;
    row.QCRowCount         = qcResult.RowCount;
    row.QCRowMatched       = qcResult.RowMatched;
    row.QCAmbiguous        = qcResult.Ambiguous;
    row.QCReadOK           = qcResult.ReadOK;
    row.QCUnmappedColumns  = strjoin(qcResult.UnmappedColumns, '; ');
    row.QCNotes            = qcResult.Notes;

    if row.QCFileExists
        qc_inspect_log(end+1, :) = { ...
            row.FileName, qcPath, strjoin(qcResult.SheetNames, '; '), ...
            qcResult.ChosenSheet, strjoin(qcResult.ColumnNames, '; '), ...
            qcResult.RowCount, qcResult.RowMatched, qcResult.Ambiguous, ...
            strjoin(qcResult.UnmappedColumns, '; '), qcResult.Notes}; %#ok<AGROW>
    end

    % ---- objective validation flags (never used to exclude anything) ----
    row.LoadOK = strcmp(row.LoadStatus, 'OK');
    row.MissingQC = ~row.QCFileExists;
    row.Check_Fs512          = row.LoadOK && (row.SamplingRateHz == EXPECTED_SAMPLING_RATE);
    row.Check_NumChannels62  = row.LoadOK && (row.NumChannels == EXPECTED_NUM_CHANNELS);
    row.Check_NoNaNInf       = row.LoadOK && (row.NaNCount == 0) && (row.InfCount == 0);
    row.Check_HasICA         = row.LoadOK && row.HasICAWeights;
    row.Check_ParseOK        = row.ParseOK;
    row.Check_ConditionValid = ismember(row.Condition, {'E', 'N', 'R'});
    row.Check_InterventionValid = ismember(row.InterventionNumber, [1 2]);
    row.MetadataOK = row.LoadOK && row.Check_Fs512 && row.Check_NumChannels62 && ...
        row.Check_NoNaNInf && row.Check_HasICA && row.Check_ParseOK && ...
        row.Check_ConditionValid && row.Check_InterventionValid;

    row.Notes = strjoin(notes, ' | ');
    rows{i} = row;
end

%% ------------------------------------------------------------------
%  4. DURATION OUTLIER FLAGS (relative to this dataset's own median)
% ------------------------------------------------------------------

durations = cellfun(@(r) r.DurationSec, rows);
valid_dur = durations(~isnan(durations));
if ~isempty(valid_dur)
    med_dur = median(valid_dur);
else
    med_dur = NaN;
end
for i = 1:n
    d = rows{i}.DurationSec;
    if ~isnan(d) && ~isnan(med_dur) && med_dur > 0
        rows{i}.DurationDeviationFromMedianPct = 100 * (d - med_dur) / med_dur;
        rows{i}.DurationOutlier = abs(d - med_dur) / med_dur > DURATION_OUTLIER_FRACTION;
    else
        rows{i}.DurationDeviationFromMedianPct = NaN;
        rows{i}.DurationOutlier = false;
    end
    if ~isnan(d)
        rows{i}.DurationDeviationFromNominal600Sec = d - NOMINAL_DURATION_SEC;
    else
        rows{i}.DurationDeviationFromNominal600Sec = NaN;
    end
end

%% ------------------------------------------------------------------
%  5. CHANNEL-LABEL CONSISTENCY CHECK (vs. first successfully loaded file)
% ------------------------------------------------------------------

ok_idx = find(strcmp(cellfun(@(r) r.LoadStatus, rows, 'UniformOutput', false), 'OK'));
ref_labels = {};
if ~isempty(ok_idx)
    ref_labels = chanlabel_store{ok_idx(1)};
end
for i = 1:n
    if strcmp(rows{i}.LoadStatus, 'OK') && ~isempty(ref_labels)
        rows{i}.ChannelLabelsConsistentWithReference = isequal(sort(chanlabel_store{i}), sort(ref_labels));
    else
        rows{i}.ChannelLabelsConsistentWithReference = false;
    end
end

%% ------------------------------------------------------------------
%  6. MISSING RECORDINGS (based ONLY on observed files -- no roadmap)
%     Also sets, per observed row, whether that subject's OTHER
%     expected intervention is missing (MissingExpectedIntervention).
% ------------------------------------------------------------------

subj_all = cellfun(@(r) r.Subject, rows, 'UniformOutput', false);
int_all  = cellfun(@(r) r.InterventionNumber, rows);

validSubjIdx = ~cellfun(@isempty, subj_all);
observedSubjects = unique(subj_all(validSubjIdx));
[~, sortOrd] = sort(str2double(observedSubjects));
observedSubjects = observedSubjects(sortOrd);

missing_rows = cell(0, 2);
for si = 1:numel(observedSubjects)
    sid = observedSubjects{si};
    for k = 1:2
        found = any(strcmp(subj_all, sid) & int_all == k);
        if ~found
            missing_rows(end+1, :) = {sid, k}; %#ok<AGROW>
        end
    end
end
fullyAbsentSubjects = {};
for sidNum = 1:EXPECTED_N_SUBJECTS
    sidStr = num2str(sidNum);
    if ~any(strcmp(observedSubjects, sidStr))
        fullyAbsentSubjects{end+1} = sidStr; %#ok<AGROW>
        missing_rows(end+1, :) = {sidStr, 1}; %#ok<AGROW>
        missing_rows(end+1, :) = {sidStr, 2}; %#ok<AGROW>
    end
end

numericObserved = str2double(observedSubjects);
outOfRangeSubjects = observedSubjects(numericObserved < 1 | numericObserved > EXPECTED_N_SUBJECTS | isnan(numericObserved));

for i = 1:n
    sid = rows{i}.Subject;
    if isempty(sid)
        rows{i}.MissingExpectedIntervention = false;
        continue;
    end
    countForSubj = sum(strcmp(subj_all, sid));
    rows{i}.MissingExpectedIntervention = countForSubj < EXPECTED_INTERVENTIONS_PER_SUBJECT;
end

%% ------------------------------------------------------------------
%  7. ASSEMBLE MASTER TABLE AND SAVE
% ------------------------------------------------------------------

T = struct2table([rows{:}]);
T = normalize_text_columns(T);  % guards against char-matrix collapse (see note in 01_audit script)

writetable(T, fullfile(META_DIR, 'intervention_metadata.xlsx'));

MissingT = cell2table(missing_rows, 'VariableNames', {'Subject', 'InterventionNumber'});
MissingT = normalize_text_columns(MissingT);
writetable(MissingT, fullfile(META_DIR, 'missing_recordings_phase1.csv'));

if ~isempty(qc_inspect_log)
    QCInspectT = cell2table(qc_inspect_log, 'VariableNames', ...
        {'SetFileName', 'QCFilePath', 'SheetNamesFound', 'ChosenSheet', ...
         'ColumnNamesFound', 'RowCount', 'RowMatched', 'Ambiguous', ...
         'UnmappedColumns', 'Notes'});
    QCInspectT = normalize_text_columns(QCInspectT);
else
    QCInspectT = cell2table(cell(0, 10), 'VariableNames', ...
        {'SetFileName', 'QCFilePath', 'SheetNamesFound', 'ChosenSheet', ...
         'ColumnNamesFound', 'RowCount', 'RowMatched', 'Ambiguous', ...
         'UnmappedColumns', 'Notes'});
end
writetable(QCInspectT, fullfile(META_DIR, 'qc_excel_inspection.csv'));

save(fullfile(META_DIR, 'intervention_metadata.mat'), ...
    'T', 'MissingT', 'QCInspectT', 'QC_FIELD_MAP', 'chanlabel_store', ...
    'ref_labels', 'med_dur', 'fullyAbsentSubjects', 'outOfRangeSubjects');

%% ------------------------------------------------------------------
%  8. QC-EXCEL INSPECTION REPORT (separate text report)
% ------------------------------------------------------------------

qcRpt = fullfile(REPORT_DIR, 'phase1_qc_excel_inspection_report.txt');
fid2 = fopen(qcRpt, 'w');
fprintf(fid2, '========================================================\n');
fprintf(fid2, 'PHASE 1 -- QC EXCEL INSPECTION REPORT\n');
fprintf(fid2, 'Generated: %s\n', datestr(now));
fprintf(fid2, '========================================================\n\n');
fprintf(fid2, 'Known QC field candidates searched for (mapped ONLY if actually present):\n');
for k = 1:numel(QC_CANON_FIELDS)
    fprintf(fid2, '  %-22s <- %s\n', QC_CANON_FIELDS{k}, strjoin(QC_FIELD_MAP.(QC_CANON_FIELDS{k}), ' / '));
end
fprintf(fid2, '\n');

nWithQC = sum(cellfun(@(r) r.QCFileExists, rows));
fprintf(fid2, 'Recordings with a QC file found : %d / %d\n', nWithQC, n);
fprintf(fid2, 'Recordings with NO QC file      : %d / %d\n\n', n - nWithQC, n);

fprintf(fid2, '--- PER-FILE QC INSPECTION DETAIL ---\n');
for i = 1:n
    r = rows{i};
    fprintf(fid2, '%s :\n', r.FileName);
    if ~r.QCFileExists
        fprintf(fid2, '    No QC file found.\n\n');
        continue;
    end
    fprintf(fid2, '    QC file       : %s\n', r.QCFileName);
    fprintf(fid2, '    Sheet chosen  : %s\n', r.QCSheetName);
    fprintf(fid2, '    Row matched   : %d (ambiguous=%d, readOK=%d)\n', r.QCRowMatched, r.QCAmbiguous, r.QCReadOK);
    if ~isempty(r.QCUnmappedColumns)
        fprintf(fid2, '    Unmapped columns found (not auto-extracted -- review manually):\n');
        unmapped = strsplit(r.QCUnmappedColumns, '; ');
        for u = 1:numel(unmapped)
            if ~isempty(unmapped{u})
                fprintf(fid2, '      QC field not automatically mapped: "%s"\n', unmapped{u});
            end
        end
    end
    if ~isempty(r.QCNotes)
        fprintf(fid2, '    Notes: %s\n', r.QCNotes);
    end
    fprintf(fid2, '\n');
end
fclose(fid2);

%% ------------------------------------------------------------------
%  9. PHASE 1 METADATA REPORT (text file + Command Window)
% ------------------------------------------------------------------

rpt = fullfile(REPORT_DIR, 'phase1_metadata_report.txt');
fid = fopen(rpt, 'w');

report_lines = {};
report_lines{end+1} = '========================================================';
report_lines{end+1} = 'PHASE 1 -- INTERVENTION METADATA REPORT';
report_lines{end+1} = '========================================================';
report_lines{end+1} = '';
report_lines{end+1} = sprintf('Subjects found        : %d', numel(observedSubjects));
report_lines{end+1} = sprintf('Recordings found      : %d', n);
report_lines{end+1} = sprintf('Expected recordings   : %d (%d subjects x %d interventions)', ...
    EXPECTED_N_SUBJECTS * EXPECTED_INTERVENTIONS_PER_SUBJECT, EXPECTED_N_SUBJECTS, EXPECTED_INTERVENTIONS_PER_SUBJECT);
report_lines{end+1} = '';
report_lines{end+1} = 'Condition counts (filename-derived, normalized):';
for c = {'E', 'N', 'R'}
    report_lines{end+1} = sprintf('  %s : %d', c{1}, sum(strcmp(T.Condition, c{1}))); %#ok<AGROW>
end
report_lines{end+1} = '';
report_lines{end+1} = 'Intervention position:';
report_lines{end+1} = sprintf('  Int 1 : %d', sum(T.InterventionNumber == 1));
report_lines{end+1} = sprintf('  Int 2 : %d', sum(T.InterventionNumber == 2));
report_lines{end+1} = '';
report_lines{end+1} = sprintf('Missing recordings (%d total):', height(MissingT));
if isempty(MissingT)
    report_lines{end+1} = '  None.';
else
    for i = 1:height(MissingT)
        report_lines{end+1} = sprintf('  Sub%s, Int%d', MissingT.Subject{i}, MissingT.InterventionNumber(i)); %#ok<AGROW>
    end
end
if ~isempty(fullyAbsentSubjects)
    report_lines{end+1} = sprintf('  (Note: subject(s) %s have NO recordings at all in this dataset.)', strjoin(fullyAbsentSubjects, ', '));
end
if ~isempty(outOfRangeSubjects)
    report_lines{end+1} = sprintf('  (Note: subject ID(s) %s fall outside the expected 1-%d range -- review for a labeling issue.)', ...
        strjoin(outOfRangeSubjects, ', '), EXPECTED_N_SUBJECTS);
end
report_lines{end+1} = '';
report_lines{end+1} = 'QC Excel files:';
report_lines{end+1} = sprintf('  Found   : %d', nWithQC);
report_lines{end+1} = sprintf('  Missing : %d', n - nWithQC);
nRowMatched = sum(cellfun(@(r) r.QCFileExists && r.QCRowMatched, rows));
report_lines{end+1} = sprintf('  QC files with a successfully matched data row : %d', nRowMatched);
nAmbiguous = sum(cellfun(@(r) r.QCFileExists && r.QCAmbiguous, rows));
report_lines{end+1} = sprintf('  QC files that were ambiguous (multiple rows, could not disambiguate) : %d', nAmbiguous);
report_lines{end+1} = '';
report_lines{end+1} = 'Duration:';
if ~isnan(med_dur)
    report_lines{end+1} = sprintf('  Median   : %.2f s (%.2f min)', med_dur, med_dur/60);
    report_lines{end+1} = sprintf('  Min      : %.2f s', min(valid_dur));
    report_lines{end+1} = sprintf('  Max      : %.2f s', max(valid_dur));
end
report_lines{end+1} = sprintf('  Outliers (>%.0f%% from median) : %d', DURATION_OUTLIER_FRACTION*100, sum(T.DurationOutlier));
report_lines{end+1} = '';
report_lines{end+1} = 'Metadata problems:';
report_lines{end+1} = sprintf('  Load errors                       : %d', sum(~T.LoadOK));
report_lines{end+1} = sprintf('  Unparseable filenames              : %d', sum(~T.ParseOK));
report_lines{end+1} = sprintf('  Channel-label-set inconsistencies  : %d', sum(strcmp(T.LoadStatus,'OK') & ~T.ChannelLabelsConsistentWithReference));
report_lines{end+1} = sprintf('  Sampling-rate != %d                 : %d', EXPECTED_SAMPLING_RATE, sum(T.LoadOK & ~T.Check_Fs512));
report_lines{end+1} = sprintf('  Channel count != %d                 : %d', EXPECTED_NUM_CHANNELS, sum(T.LoadOK & ~T.Check_NumChannels62));
report_lines{end+1} = sprintf('  Recordings overall MetadataOK=false : %d', sum(~T.MetadataOK));
report_lines{end+1} = '';
report_lines{end+1} = sprintf('SetName/Filename inconsistencies : %d', sum(T.SetNameFilenameMismatch));
if any(T.SetNameFilenameMismatch)
    idxMismatch = find(T.SetNameFilenameMismatch);
    for i = idxMismatch'
        report_lines{end+1} = sprintf('  %s : SetName="%s" -> %s', T.FileName{i}, T.SetName{i}, T.SetNameMismatchDetail{i}); %#ok<AGROW>
    end
end
report_lines{end+1} = '';
report_lines{end+1} = 'Output:';
report_lines{end+1} = sprintf('  %s', fullfile(META_DIR, 'intervention_metadata.mat'));
report_lines{end+1} = sprintf('  %s', fullfile(META_DIR, 'intervention_metadata.xlsx'));
report_lines{end+1} = sprintf('  %s', fullfile(META_DIR, 'missing_recordings_phase1.csv'));
report_lines{end+1} = sprintf('  %s', fullfile(META_DIR, 'qc_excel_inspection.csv'));
report_lines{end+1} = sprintf('  %s', rpt);
report_lines{end+1} = sprintf('  %s', qcRpt);
report_lines{end+1} = '========================================================';

for L = 1:numel(report_lines)
    fprintf('%s\n', report_lines{L});
    fprintf(fid, '%s\n', report_lines{L});
end
fclose(fid);

%% ========================================================================
%  LOCAL FUNCTIONS
%  (must appear here, after all top-level script code, per MATLAB rules
%  for scripts that mix executable statements with local functions)
% ========================================================================

function [subj, intnum, condRaw, ok] = parse_filename(fname)
    % Expected convention: sub<N>_in<K>_<COND>_pp.set
    % Fully case-insensitive; tolerant of "in"/"int" and leading zeros.
    % subj is CHAR (e.g. '16'), matching 01_audit_intervention_data.m.
    FNAME_PATTERN = '^[Ss]ub0*(\d+)_[Ii][Nn][Tt]?0*(\d+)_([A-Za-z]+)_pp$';
    [~, base, ~] = fileparts(fname);
    tok = regexp(base, FNAME_PATTERN, 'tokens', 'once');
    if isempty(tok)
        subj = ''; intnum = NaN; condRaw = ''; ok = false;
    else
        subj    = tok{1};
        intnum  = str2double(tok{2});
        condRaw = tok{3};
        ok = true;
    end
end

function cond = normalize_condition(raw)
    if isempty(raw)
        cond = '';
        return;
    end
    c = upper(raw);
    switch c
        case {'E', 'I'}
            cond = 'E';
        case 'N'
            cond = 'N';
        case 'R'
            cond = 'R';
        otherwise
            cond = '';
    end
end

function [s_subj, s_int, s_cond, s_ok] = parse_setname(setname)
    s_subj = ''; s_int = NaN; s_cond = ''; s_ok = false;
    if isempty(setname)
        return;
    end
    core = regexprep(setname, '_ICA_clean$', '', 'ignorecase');
    pat = '^[Ss]ub0*(\d+)_[Ii][Nn]_?0*(\d+)(?:_([A-Za-z]))?$';
    tok = regexp(core, pat, 'tokens', 'once');
    if isempty(tok)
        return;
    end
    s_subj = tok{1};
    s_int  = str2double(tok{2});
    if numel(tok) >= 3 && ~isempty(tok{3})
        s_cond = tok{3};
    end
    s_ok = true;
end

function [qcPath, qcFileName, howFound] = find_qc_file(folderPath, intnum)
    qcPath = ''; qcFileName = ''; howFound = '';
    if isnan(intnum)
        return;
    end
    exactName = sprintf('Preprocessing_QC_%d.xlsx', intnum);
    exactPath = fullfile(folderPath, exactName);
    if exist(exactPath, 'file') == 2
        qcPath = exactPath; qcFileName = exactName; howFound = 'exact';
        return;
    end
    d = dir(fullfile(folderPath, '*.xlsx'));
    for k = 1:numel(d)
        nm = d(k).name;
        if ~isempty(regexpi(nm, ['QC.*[_]?0*' num2str(intnum) '(\D|$)'], 'once'))
            qcPath = fullfile(folderPath, nm);
            qcFileName = nm;
            howFound = 'fuzzy';
            return;
        end
    end
end

function qc = empty_qc_result(canonFields)
    qc = struct();
    qc.SheetNames = {};
    qc.ChosenSheet = '';
    qc.ColumnNames = {};
    qc.RowCount = NaN;
    qc.RowMatched = false;
    qc.Ambiguous = false;
    qc.ReadOK = false;
    qc.UnmappedColumns = {};
    qc.Notes = '';
    qc.Mapped = struct();
    for k = 1:numel(canonFields)
        % All QC_* fields are stored as TEXT (char), '' when missing --
        % even ones that are conceptually numeric (e.g. NBadChannels).
        % This is deliberate: if some recordings have a QC value for a
        % field and others don't, a numeric-vs-char mix across rows of
        % the same table column makes struct2table fail (or behave
        % inconsistently) when the master table is assembled. Storing
        % everything as text avoids that; values are converted with
        % str2double(...) where they are used.
        qc.Mapped.(canonFields{k}) = '';
    end
end

function qc = inspect_qc_excel(qcPath, subjID, intnum, fieldMap)
    canon = fieldnames(fieldMap);
    qc = empty_qc_result(canon);

    try
        sheets = sheetnames(qcPath);
        if isstring(sheets)
            sheets = cellstr(sheets);
        end
    catch
        try
            [~, sheets] = xlsfinfo(qcPath);
        catch
            sheets = {};
        end
    end
    qc.SheetNames = sheets;

    chosen = '';
    if ~isempty(sheets)
        exactSheetIdx = find(strcmpi(sheets, 'Preprocessing_QC'), 1);
        if ~isempty(exactSheetIdx)
            chosen = sheets{exactSheetIdx};
        elseif numel(sheets) == 1
            chosen = sheets{1};
        else
            knownNames = lower(all_known_source_names(fieldMap));
            bestScore = -1; bestSheet = sheets{1};
            for si = 1:numel(sheets)
                try
                    Ttmp = readtable(qcPath, 'Sheet', sheets{si});
                    score = numel(intersect(lower(Ttmp.Properties.VariableNames), knownNames));
                catch
                    score = -1;
                end
                if score > bestScore
                    bestScore = score; bestSheet = sheets{si};
                end
            end
            chosen = bestSheet;
            if bestScore <= 0
                qc.Notes = 'Multiple sheets found, none clearly matched known QC columns; used first/best-guess sheet -- verify manually.';
            end
        end
    end
    qc.ChosenSheet = chosen;

    try
        if isempty(chosen)
            Tqc = readtable(qcPath);
        else
            Tqc = readtable(qcPath, 'Sheet', chosen);
        end
        qc.ReadOK = true;
    catch ME
        qc.Notes = strtrim([qc.Notes ' | Failed to read QC file: ' ME.message]);
        return;
    end

    qc.ColumnNames = Tqc.Properties.VariableNames;
    qc.RowCount = height(Tqc);

    if qc.RowCount == 0
        qc.Notes = strtrim([qc.Notes ' | QC sheet has zero data rows.']);
        return;
    end

    if qc.RowCount == 1
        rowIdx = 1;
    else
        vn = qc.ColumnNames;
        candidateIdx = 1:qc.RowCount;
        subjCol = find(strcmpi(vn, 'Subject'), 1);
        if ~isempty(subjCol)
            subjVals = arrayfun(@(r) extract_scalar_as_str(Tqc.(vn{subjCol})(r)), candidateIdx, 'UniformOutput', false);
            keep = strcmpi(subjVals, subjID);
            if any(keep)
                candidateIdx = candidateIdx(keep);
            end
        end
        if numel(candidateIdx) > 1
            taskNumCol = find(strcmpi(vn, 'TaskNumber'), 1);
            if ~isempty(taskNumCol)
                tVals = arrayfun(@(r) extract_scalar_as_str(Tqc.(vn{taskNumCol})(r)), candidateIdx, 'UniformOutput', false);
                keep2 = strcmp(tVals, num2str(intnum));
                if any(keep2)
                    candidateIdx = candidateIdx(keep2);
                end
            end
        end
        if numel(candidateIdx) == 1
            rowIdx = candidateIdx;
        else
            qc.Ambiguous = true;
            qc.Notes = strtrim([qc.Notes sprintf(' | QC sheet has %d rows; could not uniquely match to this recording.', qc.RowCount)]);
            return;
        end
    end

    qc.RowMatched = true;
    row1 = Tqc(rowIdx, :);

    mappedCols = {};
    for k = 1:numel(canon)
        candidates = fieldMap.(canon{k});
        hitCol = '';
        for c = 1:numel(candidates)
            hi = find(strcmpi(qc.ColumnNames, candidates{c}), 1);
            if ~isempty(hi)
                hitCol = qc.ColumnNames{hi};
                break;
            end
        end
        if ~isempty(hitCol)
            qc.Mapped.(canon{k}) = extract_scalar_as_str(row1.(hitCol));
            mappedCols{end+1} = hitCol; %#ok<AGROW>
        end
    end
    qc.UnmappedColumns = setdiff(qc.ColumnNames, mappedCols, 'stable');
end

function v = extract_scalar(val)
    if iscell(val) && isscalar(val)
        v = val{1};
    elseif isstring(val) && isscalar(val)
        v = char(val);
    else
        v = val;
    end
end

function s = extract_scalar_as_str(val)
    v = extract_scalar(val);
    if isnumeric(v) || islogical(v)
        s = num2str(v);
    elseif ischar(v)
        s = v;
    else
        s = '';
    end
end

function names = all_known_source_names(fieldMap)
    canon = fieldnames(fieldMap);
    names = {};
    for k = 1:numel(canon)
        names = [names, fieldMap.(canon{k})]; %#ok<AGROW>
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
    % EEG.event.type may be numeric or char -- handle both rather than
    % assuming, since num2str() errors on char input.
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
