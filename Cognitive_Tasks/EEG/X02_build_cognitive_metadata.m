%% ========================================================================
% X02_build_cognitive_metadata.m
%
% Purpose: build a single, clean, auditable metadata table (one row per
% non-VI .set recording, with Task/BlockNumber columns), plus a QC-Excel
% inspection report, so that X03 (segmentation + PSD) can be designed
% from verified facts rather than assumptions. Read-only; no new
% preprocessing is performed.
%
% Method: for each non-VI .set file, extracts EEGLAB metadata and
% cross-checks SetName against the filename; discovers and best-effort
% maps any QC-Excel file keyed on (BlockNumber, TaskType), never guessing
% unmapped columns (most/all recordings are expected to have no QC-Excel
% file, which is handled as a non-fatal, reported condition); applies
% objective (never exclusionary) validation flags, including
% Check_ConditionValid (Baseline + E/N/R) and Check_BlockValid
% (BlockNumber in {1,2,3}); flags duration outliers and missing
% recordings (Subject x TaskType x Block, 9 expected combinations per
% subject), scoped to the documented 1-30 subject range, separately per
% task type (Full/N-back/SART are structurally different recording
% types). *_VI.set files are discovered, counted and logged, never
% loaded. The condition-letter convention: B/b -> Baseline (Block 1,
% pre-intervention), E/e/I/i -> E (Instagram Reels), N/n -> N (Neutral
% images), R/r -> R (Rest/eyes closed); the raw letter is always
% preserved in ConditionLetterRaw. Baseline rows are built and reported
% here like any other condition; exclusion from the E/N/R inferential
% model happens only downstream, in X05/X06/X08/X09.
%
% This script does not filter/re-reference/run ICA/remove ICs/
% interpolate/resample/remove channels, segment EEG into analysis
% windows, compute PSD or band power, run any statistics, or modify any
% original .set file.
%
% Requirements: EEGLAB on the MATLAB path, initialized enough for
% pop_loadset; base MATLAB table/Excel I/O (readtable, writetable,
% sheetnames, R2019a+, with an xlsfinfo fallback); MATLAB R2016b or later
% (local functions in a script file).
%
% Inputs: DATA_DIR (preprocessed .set/.fdt files), X01's audit output.
%
% Outputs (nothing under DATA_DIR is ever written to):
%   <ANALYSIS_ROOT>/02_metadata/
%       cognitive_metadata.mat      (master table + QC inspection log)
%       cognitive_metadata.xlsx     (master table, human-readable)
%       missing_recordings_phase1.csv
%       qc_excel_inspection.csv
%   <ANALYSIS_ROOT>/09_reports/
%       phase1_metadata_report_cognitive.txt
%       phase1_qc_excel_inspection_report_cognitive.txt
%
% RUN_MODE: not applicable; this script runs once.
% ========================================================================

clear; clc;

%% ------------------------------------------------------------------
%  0. USER-EDITABLE SETTINGS
% ------------------------------------------------------------------

% >>> EDIT THIS to point at your cognitive-task data root <<<
config;
DATA_ROOT = DATA_DIR;

EXPECTED_SAMPLING_RATE = 512;
EXPECTED_NUM_CHANNELS  = 62;
EXPECTED_N_SUBJECTS    = 30;
EXPECTED_BLOCKS_PER_SUBJECT = 3;
EXPECTED_TASK_TYPES = {'full', 'nback', 'sart'};
DURATION_OUTLIER_FRACTION = 0.10;   % computed per task type -- see header item 5

[DATA_ROOT_PARENT, ~] = fileparts(DATA_ROOT);
ANALYSIS_ROOT = OUT_DIR;
META_DIR   = fullfile(ANALYSIS_ROOT, '02_metadata');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');
if ~exist(META_DIR, 'dir');   mkdir(META_DIR);   end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

fprintf('========================================================\n');
fprintf('PHASE 1 -- BUILD COGNITIVE-TASK METADATA (READ-ONLY)\n');
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
%  1. KNOWN QC-EXCEL FIELD MAP (identical candidate list to the
%     intervention pipeline; see header item 3 for why no QC files are
%     expected to be found for this dataset, and why that is non-fatal)
% ------------------------------------------------------------------

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
%  2. DISCOVER .set FILES (recursive) -- exclude *_VI.set (out of scope
%     for this analysis; logged only, identical policy to X01)
% ------------------------------------------------------------------

all_set_list = dir(fullfile(DATA_ROOT, '**', '*.set'));
all_set_list = all_set_list(~[all_set_list.isdir]);

is_vi = false(numel(all_set_list), 1);
for i = 1:numel(all_set_list)
    [~, base, ~] = fileparts(all_set_list(i).name);
    if ~isempty(regexp(base, '_VI$', 'once'))
        is_vi(i) = true;
    end
end
vi_list  = all_set_list(is_vi);
set_list = all_set_list(~is_vi);

fprintf('Found %d total .set file(s) under DATA_ROOT (%d VI.set excluded as out of scope; %d remain).\n\n', ...
    numel(all_set_list), numel(vi_list), numel(set_list));
if isempty(set_list)
    error('No non-VI .set files found under DATA_ROOT. Check the path.');
end

%% ------------------------------------------------------------------
%  3. MAIN LOOP -- PARSE, LOAD METADATA, LOCATE/INSPECT QC, VALIDATE
% ------------------------------------------------------------------

n = numel(set_list);
rows = cell(n, 1);
chanlabel_store = cell(n, 1);
qc_inspect_log = cell(0, 10);

for i = 1:n

    fpath = fullfile(set_list(i).folder, set_list(i).name);
    fprintf('[%3d/%3d] %s ... ', i, n, set_list(i).name);

    row = struct();
    row.FileName   = set_list(i).name;
    row.FolderPath = set_list(i).folder;
    row.FullPath   = fpath;

    % ---- filename parsing ----
    [subj, seqRaw, condRaw, blockNum, taskType, parseOK] = parse_cognitive_filename(set_list(i).name);
    row.ParseOK = parseOK;
    row.Subject = subj;
    row.TaskType = taskType;
    row.BlockNumber = blockNum;
    row.TaskOrderSequenceRaw = seqRaw;
    row.ConditionLetterRaw = condRaw;
    row.Condition = normalize_condition(condRaw);   % 'Baseline'/'E'/'N'/'R'/''

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

        % ---- SetName vs. filename cross-check (report only, never fix) ----
        [s_subj, s_cond, s_ok] = parse_setname_cognitive(EEG.setname);
        row.SetNameParsed = s_ok;
        mismatchParts = {};
        if s_ok
            if ~isempty(s_subj) && ~isempty(subj) && ~strcmp(s_subj, subj)
                mismatchParts{end+1} = sprintf('subject(file=%s,setname=%s)', subj, s_subj); %#ok<AGROW>
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

        fprintf('OK  (Fs=%g, chans=%d, dur=%.1f s, task=%s, block=%d, cond=%s)\n', ...
            EEG.srate, EEG.nbchan, row.DurationSec, taskType, blockNum, row.Condition);

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

    % ---- locate + inspect the corresponding QC Excel file (keyed on
    %      BlockNumber + TaskType instead of InterventionNumber -- see
    %      header item 3) ----
    [qcPath, qcFileName, qcHowFound] = find_qc_file(set_list(i).folder, blockNum, taskType);
    if isempty(qcPath)
        row.QCFileExists = false;
        row.QCFileName = '';
        row.QCFilePath = '';
        qcResult = empty_qc_result(QC_CANON_FIELDS);
        qcResult.Notes = 'No matching Preprocessing_QC file found for this block/task in this folder.';
    else
        row.QCFileExists = true;
        row.QCFileName = qcFileName;
        row.QCFilePath = qcPath;
        try
            qcResult = inspect_qc_excel(qcPath, subj, blockNum, QC_FIELD_MAP);
        catch ME2
            qcResult = empty_qc_result(QC_CANON_FIELDS);
            qcResult.Notes = ['Unhandled error inspecting QC file: ' ME2.message];
        end
        if ~strcmp(qcHowFound, 'exact')
            qcResult.Notes = strtrim([qcResult.Notes ' | Matched by fuzzy filename search, not the exact expected name.']);
        end
    end

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
    row.Check_ConditionValid = ismember(row.Condition, {'Baseline','E','N','R'});
    row.Check_BlockValid     = ismember(row.BlockNumber, [1 2 3]);
    row.MetadataOK = row.LoadOK && row.Check_Fs512 && row.Check_NumChannels62 && ...
        row.Check_NoNaNInf && row.Check_HasICA && row.Check_ParseOK && ...
        row.Check_ConditionValid && row.Check_BlockValid;

    row.Notes = strjoin(notes, ' | ');
    rows{i} = row;
end

%% ------------------------------------------------------------------
%  4. DURATION OUTLIER FLAGS -- SEPARATELY PER TASK TYPE (identical
%     rationale to X01)
% ------------------------------------------------------------------

med_dur_by_task = struct('full', NaN, 'nback', NaN, 'sart', NaN);
for t = 1:numel(EXPECTED_TASK_TYPES)
    tt = EXPECTED_TASK_TYPES{t};
    idxThisTask = find(cellfun(@(r) strcmp(r.TaskType, tt), rows));
    dThisTask = cellfun(@(r) r.DurationSec, rows(idxThisTask));
    validDThisTask = dThisTask(~isnan(dThisTask));
    if isempty(validDThisTask)
        continue;
    end
    medThisTask = median(validDThisTask);
    med_dur_by_task.(tt) = medThisTask;
    for k = 1:numel(idxThisTask)
        ii = idxThisTask(k);
        d = rows{ii}.DurationSec;
        if ~isnan(d) && medThisTask > 0
            rows{ii}.DurationDeviationFromMedianPct = 100 * (d - medThisTask) / medThisTask;
            rows{ii}.DurationOutlier = abs(d - medThisTask) / medThisTask > DURATION_OUTLIER_FRACTION;
        else
            rows{ii}.DurationDeviationFromMedianPct = NaN;
            rows{ii}.DurationOutlier = false;
        end
    end
end
for i = 1:n
    if ~isfield(rows{i}, 'DurationOutlier')
        rows{i}.DurationDeviationFromMedianPct = NaN;
        rows{i}.DurationOutlier = false;
    end
end

%% ------------------------------------------------------------------
%  5. CHANNEL-LABEL CONSISTENCY CHECK (vs. first successfully loaded
%     file, across ALL task types -- identical policy to X01)
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
%  6. MISSING RECORDINGS -- Subject x TaskType x Block (based ONLY on
%     observed files, scanning the DOCUMENTED 1..EXPECTED_N_SUBJECTS
%     range, identical policy to X01/the intervention pipeline)
% ------------------------------------------------------------------

subj_all  = cellfun(@(r) r.Subject, rows, 'UniformOutput', false);
task_all  = cellfun(@(r) r.TaskType, rows, 'UniformOutput', false);
block_all = cellfun(@(r) r.BlockNumber, rows);

validSubjIdx = ~cellfun(@isempty, subj_all);
observedSubjects = unique(subj_all(validSubjIdx));
[~, sortOrd] = sort(str2double(observedSubjects));
observedSubjects = observedSubjects(sortOrd);

missing_rows = cell(0, 3);
for sidNum = 1:EXPECTED_N_SUBJECTS
    sidStr = num2str(sidNum);
    for t = 1:numel(EXPECTED_TASK_TYPES)
        tt = EXPECTED_TASK_TYPES{t};
        for b = 1:EXPECTED_BLOCKS_PER_SUBJECT
            found = any(strcmp(subj_all, sidStr) & strcmp(task_all, tt) & block_all == b);
            if ~found
                missing_rows(end+1, :) = {sidStr, tt, b}; %#ok<AGROW>
            end
        end
    end
end
fullyAbsentSubjects = setdiff(arrayfun(@num2str, 1:EXPECTED_N_SUBJECTS, 'UniformOutput', false), observedSubjects, 'stable');

numericObserved = str2double(observedSubjects);
outOfRangeSubjects = observedSubjects(numericObserved < 1 | numericObserved > EXPECTED_N_SUBJECTS | isnan(numericObserved));

for i = 1:n
    sid = rows{i}.Subject; tt = rows{i}.TaskType;
    if isempty(sid) || isempty(tt)
        rows{i}.MissingSiblingBlocksForThisSubjectTask = false;
        continue;
    end
    countForSubjTask = sum(strcmp(subj_all, sid) & strcmp(task_all, tt));
    rows{i}.MissingSiblingBlocksForThisSubjectTask = countForSubjTask < EXPECTED_BLOCKS_PER_SUBJECT;
end

%% ------------------------------------------------------------------
%  7. ASSEMBLE MASTER TABLE AND SAVE
% ------------------------------------------------------------------

T = struct2table([rows{:}]);
T = normalize_text_columns(T);

writetable(T, fullfile(META_DIR, 'cognitive_metadata.xlsx'));

MissingT = cell2table(missing_rows, 'VariableNames', {'Subject', 'TaskType', 'MissingBlockNumber'});
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

% ref_labels is saved under this same plain name for downstream scripts
% (X03/X07) to load, consistent with the intervention pipeline's own
% variable naming.
save(fullfile(META_DIR, 'cognitive_metadata.mat'), ...
    'T', 'MissingT', 'QCInspectT', 'QC_FIELD_MAP', 'chanlabel_store', ...
    'ref_labels', 'med_dur_by_task', 'fullyAbsentSubjects', 'outOfRangeSubjects');

%% ------------------------------------------------------------------
%  8. QC-EXCEL INSPECTION REPORT
% ------------------------------------------------------------------

qcRpt = fullfile(REPORT_DIR, 'phase1_qc_excel_inspection_report_cognitive.txt');
fid2 = fopen(qcRpt, 'w');
fprintf(fid2, '========================================================\n');
fprintf(fid2, 'PHASE 1 -- QC EXCEL INSPECTION REPORT (COGNITIVE-TASK DATA)\n');
fprintf(fid2, 'Generated: %s\n', datestr(now));
fprintf(fid2, '========================================================\n\n');
fprintf(fid2, 'Known QC field candidates searched for (mapped ONLY if actually present):\n');
for k = 1:numel(QC_CANON_FIELDS)
    fprintf(fid2, '  %-22s <- %s\n', QC_CANON_FIELDS{k}, strjoin(QC_FIELD_MAP.(QC_CANON_FIELDS{k}), ' / '));
end
fprintf(fid2, '\n');

nWithQC = sum(cellfun(@(r) r.QCFileExists, rows));
fprintf(fid2, 'Recordings with a QC file found : %d / %d\n', nWithQC, n);
fprintf(fid2, 'Recordings with NO QC file      : %d / %d\n', n - nWithQC, n);
fprintf(fid2, '(If this is 0 found / %d missing: no Preprocessing_QC*.xlsx files were visible\n', n);
fprintf(fid2, ' in the folder listings the user provided for this dataset -- see script header\n');
fprintf(fid2, ' item 3. This is expected and non-fatal; every other output in this pipeline is\n');
fprintf(fid2, ' unaffected by a missing QC file.)\n\n');

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

rpt = fullfile(REPORT_DIR, 'phase1_metadata_report_cognitive.txt');
fid = fopen(rpt, 'w');

report_lines = {};
report_lines{end+1} = '========================================================';
report_lines{end+1} = 'PHASE 1 -- COGNITIVE-TASK METADATA REPORT';
report_lines{end+1} = '========================================================';
report_lines{end+1} = '';
report_lines{end+1} = sprintf('Subjects found         : %d', numel(observedSubjects));
report_lines{end+1} = sprintf('Recordings found (non-VI): %d', n);
report_lines{end+1} = sprintf('VI.set files (out of scope): %d', numel(vi_list));
report_lines{end+1} = sprintf('Expected recordings    : %d (%d subjects x %d task types x %d blocks)', ...
    EXPECTED_N_SUBJECTS * numel(EXPECTED_TASK_TYPES) * EXPECTED_BLOCKS_PER_SUBJECT, ...
    EXPECTED_N_SUBJECTS, numel(EXPECTED_TASK_TYPES), EXPECTED_BLOCKS_PER_SUBJECT);
report_lines{end+1} = '';
report_lines{end+1} = 'Condition counts (filename-derived, normalized; includes Baseline):';
for c = {'Baseline', 'E', 'N', 'R'}
    report_lines{end+1} = sprintf('  %s : %d', c{1}, sum(strcmp(T.Condition, c{1}))); %#ok<AGROW>
end
report_lines{end+1} = '';
report_lines{end+1} = 'Task type / block counts:';
for t = 1:numel(EXPECTED_TASK_TYPES)
    tt = EXPECTED_TASK_TYPES{t};
    for b = 1:EXPECTED_BLOCKS_PER_SUBJECT
        report_lines{end+1} = sprintf('  %s, Block %d : %d', tt, b, sum(strcmp(T.TaskType, tt) & T.BlockNumber == b)); %#ok<AGROW>
    end
end
report_lines{end+1} = '';
report_lines{end+1} = sprintf('Missing Subject x TaskType x Block combinations (%d total):', height(MissingT));
if isempty(MissingT)
    report_lines{end+1} = '  None.';
else
    for i = 1:height(MissingT)
        report_lines{end+1} = sprintf('  Sub%s, %s, Block %d', MissingT.Subject{i}, MissingT.TaskType{i}, MissingT.MissingBlockNumber(i)); %#ok<AGROW>
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
report_lines{end+1} = sprintf('  QC files that were ambiguous : %d', nAmbiguous);
report_lines{end+1} = '';
report_lines{end+1} = 'Duration (computed separately PER TASK TYPE -- see header item 5):';
for t = 1:numel(EXPECTED_TASK_TYPES)
    tt = EXPECTED_TASK_TYPES{t};
    dThis = T.DurationSec(strcmp(T.TaskType, tt) & strcmp(T.LoadStatus, 'OK'));
    if isempty(dThis)
        report_lines{end+1} = sprintf('  %s: no successfully loaded files.', tt); %#ok<AGROW>
        continue;
    end
    report_lines{end+1} = sprintf('  %s: median=%.2fs min=%.2fs max=%.2fs outliers(>%.0f%%)=%d', ...
        tt, med_dur_by_task.(tt), min(dThis), max(dThis), DURATION_OUTLIER_FRACTION*100, ...
        sum(T.DurationOutlier & strcmp(T.TaskType, tt))); %#ok<AGROW>
end
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
report_lines{end+1} = sprintf('  %s', fullfile(META_DIR, 'cognitive_metadata.mat'));
report_lines{end+1} = sprintf('  %s', fullfile(META_DIR, 'cognitive_metadata.xlsx'));
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

fprintf('\nSTOP: review the metadata table and reports before running X03.\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%  (must appear here, after all top-level script code, per MATLAB rules)
% ========================================================================

function [subj, seqRaw, condRaw, blockNum, taskType, ok] = parse_cognitive_filename(fname)
    % IDENTICAL regex/logic to X01_audit_cognitive_data.m's
    % parse_cognitive_filename() -- duplicated here (rather than shared
    % via a path-added function file) to keep this script fully
    % self-contained and runnable on its own, matching the intervention
    % pipeline's own convention of each phase script being standalone.
    FNAME_PATTERN = '^[Ss]ub0*(\d+)_([A-Za-z]{3})_([1-3])_([A-Za-z]+)_pp$';
    [~, base, ~] = fileparts(fname);
    tok = regexp(base, FNAME_PATTERN, 'tokens', 'once');
    if isempty(tok)
        subj = ''; seqRaw = ''; condRaw = ''; blockNum = NaN; taskType = ''; ok = false;
        return;
    end
    subj = tok{1};
    code3 = tok{2};
    blockNum = str2double(tok{3});
    taskTypeRaw = lower(tok{4});
    if numel(code3) ~= 3
        subj = ''; seqRaw = ''; condRaw = ''; blockNum = NaN; taskType = ''; ok = false;
        return;
    end
    seqRaw = code3(1:2);
    condRaw = code3(3);
    validTaskTypes = {'full','nback','sart'};
    if ~ismember(taskTypeRaw, validTaskTypes)
        ok = false;
        taskType = taskTypeRaw;
        return;
    end
    taskType = taskTypeRaw;
    ok = true;
end

function cond = normalize_condition(raw)
    % Mapping (extends the intervention pipeline's own convention):
    % B/b -> Baseline, E/e/I/i -> E, N/n -> N, R/r -> R.
    if isempty(raw)
        cond = '';
        return;
    end
    c = upper(raw);
    switch c
        case 'B'
            cond = 'Baseline';
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

function [s_subj, s_cond, s_ok] = parse_setname_cognitive(setname)
    % Best-effort parse of EEG.setname for the cognitive-task convention.
    % Never corrects anything -- only extracts what it can so the calling
    % code can compare it against the filename-derived values. If the
    % actual SetName convention used by this dataset's preprocessing
    % differs from the filename convention, this will simply report
    % s_ok=false and no cross-check is attempted for that file (never
    % guessed).
    s_subj = ''; s_cond = ''; s_ok = false;
    if isempty(setname)
        return;
    end
    core = regexprep(setname, '_ICA_clean$', '', 'ignorecase');
    pat = '^[Ss]ub0*(\d+)_([A-Za-z]{3})_[1-3]_[A-Za-z]+$';
    tok = regexp(core, pat, 'tokens', 'once');
    if isempty(tok)
        return;
    end
    s_subj = tok{1};
    code3 = tok{2};
    if numel(code3) == 3
        s_cond = code3(3);
    end
    s_ok = true;
end

function [qcPath, qcFileName, howFound] = find_qc_file(folderPath, blockNum, taskType)
    % Locates the QC Excel file for a given block/task in the SAME folder
    % as the .set file. Keyed on BlockNumber (analogous to the
    % intervention pipeline's InterventionNumber) with TaskType as an
    % additional fuzzy-match hint. See header item 3: no QC files were
    % observed for this dataset in the materials provided, so this is
    % expected to return empty for most/all files -- handled gracefully
    % by the caller, exactly like the intervention pipeline.
    qcPath = ''; qcFileName = ''; howFound = '';
    if isnan(blockNum)
        return;
    end
    exactName = sprintf('Preprocessing_QC_%d.xlsx', blockNum);
    exactPath = fullfile(folderPath, exactName);
    if exist(exactPath, 'file') == 2
        qcPath = exactPath; qcFileName = exactName; howFound = 'exact';
        return;
    end
    d = dir(fullfile(folderPath, '*.xlsx'));
    for k = 1:numel(d)
        nm = d(k).name;
        if ~isempty(regexpi(nm, ['QC.*[_]?0*' num2str(blockNum) '(\D|$)'], 'once')) || ...
           (~isempty(taskType) && ~isempty(regexpi(nm, ['QC.*' taskType], 'once')))
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
        qc.Mapped.(canonFields{k}) = '';
    end
end

function qc = inspect_qc_excel(qcPath, subjID, blockNum, fieldMap)
    % IDENTICAL logic to the intervention pipeline's inspect_qc_excel(),
    % with "TaskNumber" matched against BlockNumber instead of
    % InterventionNumber (the QC schema itself, if such files exist for
    % this dataset, is not independently verified -- this candidate
    % mapping is preserved verbatim rather than guessed anew).
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
                qc.Notes = 'Multiple sheets found, none clearly matched known QC columns; used first/best-guess sheet.';
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
                keep2 = strcmp(tVals, num2str(blockNum));
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
