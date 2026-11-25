%% ========================================================================
% X01_audit_cognitive_data.m
%
% Purpose: read-only audit of the cognitive-task EEG data (Full-session /
% N-back / SART recordings, each measured under Baseline + Instagram
% Reels [E] / Neutral Images [N] / Rest [R]) before any decision is made
% about segmentation, PSD estimation, or statistics.
%
% Method: recursively discovers all *_pp.set files and *_VI.set files
% under DATA_ROOT (VI files are out of scope for this analysis and are
% logged only, never loaded or parsed further); loads each *_pp.set
% dataset's metadata (and full data, to check NaN/Inf) via EEGLAB without
% modifying or re-saving anything; parses Subject / TaskBlockSequence /
% ConditionLetter / BlockNumber / TaskType from the filename (see Section
% 3 for the naming convention: sub<N>_<3-letter code>_<block 1-3>_
% <full|nback|sart>_pp.set); cross-references the 30-subject roadmap
% table; builds a master audit table and a human-readable report; flags
% (but does not silently resolve) any contradictions, missing recordings,
% duplicates, or unexpected properties. Duration-outlier flags are
% computed separately per task type (full/nback/sart), since durations
% differ structurally between a full session and a single N-back or SART
% block; the check itself (flag only, never reject, relative to the
% group's own median) is unchanged from the intervention-data audit.
% Baseline (Block 1) is audited and reported like any other condition; it
% is excluded from the E/N/R inferential model only downstream, in
% X05/X06/X08/X09.
%
% This script does not segment the data into windows, compute any PSD or
% spectral quantity, run any statistics, re-filter/re-reference/re-run
% ICA/interpolate anything, or overwrite/modify any original .set file.
%
% Requirements: EEGLAB on the MATLAB path and initialized enough that
% pop_loadset is available (`eeglab('nogui')` or `eeglab;`); MATLAB
% R2016b or later (uses local functions within a script file).
%
% Inputs: DATA_DIR (preprocessed .set/.fdt files).
%
% Outputs (nothing under DATA_DIR is ever written to):
%   <ANALYSIS_ROOT>/01_audit/
%       audit_table_cognitive.csv
%       audit_table_cognitive.mat   (also contains channel-label lists,
%                                     event-type detail, missing/duplicate
%                                     sub-tables, and the VI-file list)
%       duplicate_recordings.csv
%       vi_files_out_of_scope.csv
%   <ANALYSIS_ROOT>/09_reports/
%       phase0_audit_report_cognitive.txt
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

% Reference expectations, used ONLY to flag deviations for your review.
% These are NOT enforced/assumed anywhere in the processing below.
EXPECTED_SAMPLING_RATE = 512;   % Hz, per the preprocessing pipeline (SR = 512), unchanged
EXPECTED_NUM_CHANNELS  = 62;    % scalp channels after A1/A2 removal, unchanged
EXPECTED_N_SUBJECTS    = 30;
EXPECTED_BLOCKS_PER_SUBJECT = 3;      % Block 1 (Baseline), Block 2, Block 3
EXPECTED_TASK_TYPES = {'full', 'nback', 'sart'};   % each block should have all three
DURATION_OUTLIER_FRACTION = 0.10;  % flag files whose duration deviates from
                                    % THAT TASK TYPE's own median by >10%
                                    % (see header item 5 -- computed per task type)

% Whether the hard-coded roadmap table in Section 1 is trusted to
% correspond, row-for-row, to the Subject numbers embedded in the actual
% filenames. IDENTICAL caveat and default as the intervention pipeline's
% 01_audit_intervention_data.m -- this project's subjects and roadmap are
% shared with the intervention recordings, but the actual cognitive-task
% recording order is not independently verified to follow this table's
% row order either. Cross-checking is therefore DISABLED by default.
ALLOC_CROSS_VALIDATION_TRUSTED = false;

% Output location: sibling of DATA_ROOT, mirroring the intervention
% pipeline's layout convention exactly, but under a NEW root name so the
% original intervention_analysis output is never touched or overwritten.
[DATA_ROOT_PARENT, ~] = fileparts(DATA_ROOT);
ANALYSIS_ROOT = OUT_DIR;
AUDIT_DIR  = fullfile(ANALYSIS_ROOT, '01_audit');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');

if ~exist(AUDIT_DIR, 'dir');  mkdir(AUDIT_DIR);  end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

fprintf('========================================================\n');
fprintf('PHASE 0 -- COGNITIVE-TASK DATA AUDIT (READ-ONLY)\n');
fprintf('========================================================\n');
fprintf('DATA_ROOT     : %s\n', DATA_ROOT);
fprintf('ANALYSIS_ROOT : %s  (NEW output root -- separate from intervention_analysis)\n', ANALYSIS_ROOT);
fprintf('--------------------------------------------------------\n\n');

if exist('pop_loadset', 'file') ~= 2
    error(['EEGLAB (pop_loadset) not found on the MATLAB path. ' ...
           'Start EEGLAB first, e.g.:  eeglab(''nogui'');  or  eeglab;']);
end

if ~exist(DATA_ROOT, 'dir')
    error('DATA_ROOT does not exist: %s', DATA_ROOT);
end

if ~ALLOC_CROSS_VALIDATION_TRUSTED
    fprintf(['NOTE: ALLOC_CROSS_VALIDATION_TRUSTED = false (identical policy to the\n' ...
             'intervention pipeline). The roadmap table in Section 1 is NOT used to\n' ...
             'validate/override anything parsed from the cognitive-task filenames.\n\n']);
end

%% ------------------------------------------------------------------
%  1. PLANNED COUNTERBALANCING ROADMAP (NOT ESTABLISHED AS GROUND TRUTH)
%     IDENTICAL table to 01_audit_intervention_data.m / X01_audit_
%     intervention_data.m -- same 30 subjects, same study, transcribed
%     directly from "Experimental_Design_30_Subjects_3.pdf" (columns:
%     Subject, Gender, T1, Int1, T2, Int2, T3). Never re-typed or edited.
%
%     NOTE ON RELEVANCE TO THE COGNITIVE-TASK FILENAME CONVENTION: the T1/
%     T2/T3 columns (each 'SN' or 'NS') are exactly the two-letter task-
%     order-sequence part of the 3-letter cognitive-task filename code
%     used here (see Section 3), and Int1/Int2 are exactly
%     the condition letter expected for Block 2 / Block 3 respectively
%     (Block 1 is always Baseline, which this table has no column for
%     because the table predates the cognitive-task recordings). This is
%     useful corroborating context for the naming convention, but -- per
%     the SAME caveat as the intervention pipeline -- it does NOT make the
%     table trustworthy for row-to-actual-subject correspondence, so
%     ALLOC_CROSS_VALIDATION_TRUSTED stays false and this table is not
%     used to override anything parsed from the actual filenames.
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
ALLOC = normalize_text_columns(ALLOC);

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
%  2. DISCOVER .set FILES (recursive) -- *_pp.set (to be processed) and
%     *_VI.set (out of scope, logged only)
% ------------------------------------------------------------------

all_set_list = dir(fullfile(DATA_ROOT, '**', '*.set'));
all_set_list = all_set_list(~[all_set_list.isdir]);
fprintf('Found %d total .set file(s) under DATA_ROOT.\n', numel(all_set_list));

is_vi = false(numel(all_set_list), 1);
for i = 1:numel(all_set_list)
    [~, base, ~] = fileparts(all_set_list(i).name);
    if ~isempty(regexp(base, '_VI$', 'once'))
        is_vi(i) = true;
    end
end
vi_list  = all_set_list(is_vi);
set_list = all_set_list(~is_vi);

fprintf('  %d file(s) matched "*_VI.set" -- OUT OF SCOPE for this analysis: logged only, never loaded or processed.\n', numel(vi_list));
fprintf('  %d file(s) remain for parsing/loading (full_pp / nback_pp / sart_pp).\n\n', numel(set_list));

if isempty(set_list)
    error('No non-VI .set files found under DATA_ROOT. Check the path.');
end

%% ------------------------------------------------------------------
%  3. FILENAME PARSING
% ------------------------------------------------------------------
%
% Naming convention:
%
%   sub<N>_<SEQCOND>_<BLOCK>_<TASKTYPE>_pp.set
%
%   e.g. sub1_NSB_1_full_pp.set, sub1_NSB_1_nback_pp.set,
%        sub1_NSB_1_sart_pp.set, sub1_NSB_1_VI.set (out of scope)
%
%   SEQCOND (3 letters): first 2 letters = task-order sequence for that
%     block (NS = N-back-then-SART, SN = SART-then-N-back -- same
%     vocabulary as the roadmap table's T1/T2/T3 columns); 3rd letter =
%     condition code for that block (B = Baseline; I/E/N/R = the
%     intervention condition that block follows, identical letter
%     vocabulary and I->E mapping as the intervention pipeline).
%   BLOCK (1-3): 1 = Baseline (before any intervention), 2 = after
%     Intervention 1, 3 = after Intervention 2 (per the study timeline in
%     Experimental Design.md).
%   TASKTYPE: full / nback / sart -- three PARALLEL, already-preprocessed
%     dataset types for the SAME block (same underlying recording
%     session, different task segments/processing).
%
% The functions parse_cognitive_filename() and map_condition_letter() and
% get_expected_condition() that implement this (including the exact
% regex) are defined as LOCAL FUNCTIONS at the end of this file (MATLAB
% requires all local functions in a script to follow all top-level code).

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
    [subj, seqRaw, condRaw, blockNum, taskType, parseOK] = parse_cognitive_filename(set_list(i).name);
    row.ParseOK = parseOK;
    row.Subject = subj;
    row.TaskOrderSequenceRaw = seqRaw;
    row.ConditionLetterRaw = condRaw;
    row.BlockNumber = blockNum;
    row.TaskType = taskType;

    if parseOK
        [condCandidate, needsConfirm] = map_condition_letter(condRaw);
    else
        condCandidate = ''; needsConfirm = true;
    end
    row.ConditionCandidate = condCandidate;
    row.ConditionNeedsConfirmation = needsConfirm;

    % ---- block-vs-condition structural consistency check (independent
    %      of ALLOC; a property of the naming convention itself: Block 1
    %      should always carry the Baseline code, Blocks 2/3 should
    %      always carry an intervention letter, never Baseline) ----
    blockCondOK = true;
    if parseOK
        if blockNum == 1 && ~strcmp(condCandidate, 'Baseline')
            blockCondOK = false;
        elseif ismember(blockNum, [2 3]) && strcmp(condCandidate, 'Baseline')
            blockCondOK = false;
        end
    end
    row.BlockConditionConsistent = blockCondOK;

    % ---- cross-validate against allocation table (only if trusted;
    %      identical policy to the intervention pipeline) ----
    if ALLOC_CROSS_VALIDATION_TRUSTED && parseOK && ~isempty(subj) && ismember(blockNum, [2 3])
        intnum = blockNum - 1;   % Block 2 -> Int1, Block 3 -> Int2
        expectedCond = get_expected_condition(ALLOC, subj, intnum);
    else
        expectedCond = '';
    end
    row.ConditionExpectedFromAllocation = expectedCond;

    if ~ALLOC_CROSS_VALIDATION_TRUSTED
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
        event_detail_store{i} = row.EventTypesSummary;

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
            notes{end+1} = 'EMPTY_DATA';
        end

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
        if parseOK && ~blockCondOK
            notes{end+1} = sprintf('BLOCK_CONDITION_INCONSISTENT(block=%d,cond=%s)', blockNum, condCandidate);
        end

        fprintf('OK  (Fs=%g, chans=%d, dur=%.1f s, task=%s, block=%d, cond=%s)\n', ...
            EEG.srate, EEG.nbchan, row.DurationSec, taskType, blockNum, condCandidate);

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
%  5. DURATION OUTLIER CHECK -- SEPARATELY PER TASK TYPE (see header
%     item 5: pooling across full/nback/sart would be meaningless since
%     these are structurally different recording lengths by design).
%     Flag only; NEVER reject a file for duration, per explicit
%     instruction, regardless of task type.
% ------------------------------------------------------------------

med_dur_by_task = struct('full', NaN, 'nback', NaN, 'sart', NaN);
for t = 1:numel(EXPECTED_TASK_TYPES)
    tt = EXPECTED_TASK_TYPES{t};
    idxThisTask = find(cellfun(@(r) strcmp(r.TaskType, tt), rows));
    durationsThisTask = cellfun(@(r) r.DurationSec, rows(idxThisTask));
    validDurThisTask = durationsThisTask(~isnan(durationsThisTask));
    if isempty(validDurThisTask)
        continue;
    end
    medThisTask = median(validDurThisTask);
    med_dur_by_task.(tt) = medThisTask;
    for k = 1:numel(idxThisTask)
        ii = idxThisTask(k);
        d = rows{ii}.DurationSec;
        if ~isnan(d) && medThisTask > 0 && abs(d - medThisTask) / medThisTask > DURATION_OUTLIER_FRACTION
            rows{ii}.Notes = strtrim([rows{ii}.Notes ' | ' sprintf('DURATION_OUTLIER_WITHIN_TASKTYPE(%.1fs vs %s median %.1fs)', d, tt, medThisTask)]);
            if strcmp(rows{ii}.Status, 'OK')
                rows{ii}.Status = 'REVIEW_NEEDED';
            end
        end
    end
end

%% ------------------------------------------------------------------
%  6. DUPLICATE DETECTION -- keyed on Subject x TaskType x BlockNumber
%     (the correct uniqueness key for this file layout; the intervention
%     pipeline used Subject x InterventionNumber because it had no Task
%     or Block dimension).
% ------------------------------------------------------------------

subj_all  = cellfun(@(r) r.Subject, rows, 'UniformOutput', false);
task_all  = cellfun(@(r) r.TaskType, rows, 'UniformOutput', false);
block_all = cellfun(@(r) r.BlockNumber, rows);

dup_rows = cell(0, 5);  % Subject, TaskType, BlockNumber, FileName, FullPath
valid_idx = find(~cellfun(@isempty, subj_all) & ~cellfun(@isempty, task_all) & ~isnan(block_all));
combo_keys = arrayfun(@(i) sprintf('%s_%s_%d', subj_all{i}, task_all{i}, block_all(i)), valid_idx, 'UniformOutput', false);
[u_keys, ~, ic] = unique(combo_keys);
for k = 1:numel(u_keys)
    members = valid_idx(ic == k);
    if numel(members) > 1
        for m = 1:numel(members)
            r = rows{members(m)};
            dup_rows(end+1, :) = {r.Subject, r.TaskType, r.BlockNumber, r.FileName, r.FullPath}; %#ok<AGROW>
            rows{members(m)}.Notes = strtrim([rows{members(m)}.Notes ' | DUPLICATE_RECORDING']);
            rows{members(m)}.Status = 'REVIEW_NEEDED';
        end
    end
end

%% ------------------------------------------------------------------
%  7. COMPLETENESS CHECK -- for every Subject (1..EXPECTED_N_SUBJECTS)
%     and every TaskType, are all 3 blocks present? This does NOT rely
%     on ALLOC (identical caution as the intervention pipeline's Section
%     7, which only computed "missing" against a trusted allocation --
%     here the expected key set is simply "3 blocks per Subject x
%     TaskType", which needs no allocation table at all).
% ------------------------------------------------------------------

missing_rows = cell(0, 3);  % Subject, TaskType, MissingBlockNumber
for s = 1:EXPECTED_N_SUBJECTS
    subjID = num2str(s);
    for t = 1:numel(EXPECTED_TASK_TYPES)
        tt = EXPECTED_TASK_TYPES{t};
        for b = 1:EXPECTED_BLOCKS_PER_SUBJECT
            found = any(strcmp(subj_all, subjID) & strcmp(task_all, tt) & block_all == b);
            if ~found
                missing_rows(end+1, :) = {subjID, tt, b}; %#ok<AGROW>
            end
        end
    end
end
observedSubjNums = unique(str2double(subj_all(~cellfun(@isempty, subj_all))));
outOfRangeSubj = observedSubjNums(observedSubjNums < 1 | observedSubjNums > EXPECTED_N_SUBJECTS | isnan(observedSubjNums));

%% ------------------------------------------------------------------
%  8. CONDITION MISMATCH SUB-TABLE (only meaningful while
%     ALLOC_CROSS_VALIDATION_TRUSTED = true; identical to intervention
%     pipeline's Section 8)
% ------------------------------------------------------------------

mismatch_rows = cell(0, 6);
for i = 1:n
    if strcmp(rows{i}.ConditionMatch, 'MISMATCH')
        r = rows{i};
        mismatch_rows(end+1, :) = {r.Subject, r.BlockNumber, r.FileName, ...
            r.ConditionLetterRaw, r.ConditionCandidate, r.ConditionExpectedFromAllocation}; %#ok<AGROW>
    end
end

%% ------------------------------------------------------------------
%  9. CHANNEL LABEL CONSISTENCY CHECK (identical policy to the
%     intervention pipeline: one reference set, from the first
%     successfully loaded file, compared across ALL files regardless of
%     task type -- all three task types are recorded on the same
%     65-channel g.tec system, so a single reference is still correct.)
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
            rows{i}.Notes = strtrim([rows{i}.Notes ' | CHANNEL_LABEL_SET_DIFFERS_FROM_REFERENCE']);
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

preferredOrder = {'Subject','TaskType','BlockNumber','TaskOrderSequenceRaw','ConditionLetterRaw', ...
    'ConditionCandidate','ConditionNeedsConfirmation','BlockConditionConsistent', ...
    'ConditionExpectedFromAllocation','ConditionMatch', ...
    'FileName','FolderPath','LoadStatus','LoadErrorMessage', ...
    'DurationSec','DurationMin','SamplingRateHz','NumChannels', ...
    'NumSamples','NumChannelLabels','ChanlocsHaveXYZ','NumEvents', ...
    'EventTypesSummary','HasICAWeights','NumICs','NaNCount','InfCount', ...
    'SetName','ParseOK','FullPath','Notes','Status'};
existingCols = preferredOrder(ismember(preferredOrder, T.Properties.VariableNames));
T = T(:, existingCols);
T = normalize_text_columns(T);

writetable(T, fullfile(AUDIT_DIR, 'audit_table_cognitive.csv'));

DupT = cell2table(dup_rows, 'VariableNames', {'Subject','TaskType','BlockNumber','FileName','FullPath'});
DupT = normalize_text_columns(DupT);
writetable(DupT, fullfile(AUDIT_DIR, 'duplicate_recordings.csv'));

MissingT = cell2table(missing_rows, 'VariableNames', {'Subject','TaskType','MissingBlockNumber'});
MissingT = normalize_text_columns(MissingT);
writetable(MissingT, fullfile(AUDIT_DIR, 'missing_recordings.csv'));

MismatchT = cell2table(mismatch_rows, 'VariableNames', ...
    {'Subject','BlockNumber','FileName','ConditionLetterRaw','ConditionCandidate','ConditionExpectedFromAllocation'});
MismatchT = normalize_text_columns(MismatchT);
writetable(MismatchT, fullfile(AUDIT_DIR, 'condition_mismatches.csv'));

viRows = cell(numel(vi_list), 2);
for i = 1:numel(vi_list)
    viRows(i, :) = {vi_list(i).name, fullfile(vi_list(i).folder, vi_list(i).name)};
end
ViT = cell2table(viRows, 'VariableNames', {'FileName','FullPath'});
ViT = normalize_text_columns(ViT);
writetable(ViT, fullfile(AUDIT_DIR, 'vi_files_out_of_scope.csv'));

save(fullfile(AUDIT_DIR, 'audit_table_cognitive.mat'), ...
    'T', 'DupT', 'MissingT', 'MismatchT', 'ViT', 'ALLOC', 'chanlabel_store', ...
    'ref_labels', 'label_mismatch_files', 'med_dur_by_task', 'outOfRangeSubj');

%% ------------------------------------------------------------------
%  11. HUMAN-READABLE REPORT
% ------------------------------------------------------------------

rpt = fullfile(REPORT_DIR, 'phase0_audit_report_cognitive.txt');
fid = fopen(rpt, 'w');

fprintf(fid, '========================================================\n');
fprintf(fid, 'PHASE 0 -- COGNITIVE-TASK DATA AUDIT REPORT\n');
fprintf(fid, 'Generated: %s\n', datestr(now));
fprintf(fid, '========================================================\n\n');

fprintf(fid, 'DATA_ROOT     : %s\n', DATA_ROOT);
fprintf(fid, 'ANALYSIS_ROOT : %s (new output root; intervention_analysis untouched)\n\n', ANALYSIS_ROOT);

fprintf(fid, '--- ALLOCATION-TABLE CROSS-VALIDATION STATUS ---\n');
if ALLOC_CROSS_VALIDATION_TRUSTED
    fprintf(fid, '  ENABLED.\n');
else
    fprintf(fid, '  DISABLED (identical policy to the intervention pipeline -- see script header).\n');
end
fprintf(fid, '\n');

fprintf(fid, '--- OVERVIEW ---\n');
fprintf(fid, 'Total .set files found (incl. VI)  : %d\n', numel(all_set_list));
fprintf(fid, 'VI.set files (out of scope)        : %d\n', numel(vi_list));
fprintf(fid, 'Non-VI files parsed/loaded          : %d\n', n);
fprintf(fid, 'Expected subjects                  : %d\n', EXPECTED_N_SUBJECTS);
fprintf(fid, 'Expected blocks/subject/task        : %d\n', EXPECTED_BLOCKS_PER_SUBJECT);
fprintf(fid, 'Expected task types                 : %s\n', strjoin(EXPECTED_TASK_TYPES, ', '));
fprintf(fid, 'Expected total pp recordings        : %d\n', EXPECTED_N_SUBJECTS * EXPECTED_BLOCKS_PER_SUBJECT * numel(EXPECTED_TASK_TYPES));
fprintf(fid, 'Files successfully loaded           : %d\n', sum(strcmp(T.LoadStatus, 'OK')));
fprintf(fid, 'Files with LOAD_ERROR               : %d\n', sum(strcmp(T.LoadStatus, 'LOAD_ERROR')));
fprintf(fid, 'Files with unparseable filenames    : %d\n', sum(~T.ParseOK));
fprintf(fid, 'Files with block/condition inconsistency: %d\n', sum(~T.BlockConditionConsistent));
fprintf(fid, 'Duplicate recordings (file rows)    : %d\n', height(DupT));
fprintf(fid, 'Missing Subject x Task x Block combos: %d\n', height(MissingT));
fprintf(fid, 'Subject numbers outside 1..%d observed: %s\n', EXPECTED_N_SUBJECTS, ...
    iif(isempty(outOfRangeSubj), 'none', mat2str(outOfRangeSubj')));
fprintf(fid, 'Files flagged CONTAINS_NAN          : %d\n', sum(contains(T.Notes, 'CONTAINS_NAN')));
fprintf(fid, 'Files flagged CONTAINS_INF          : %d\n', sum(contains(T.Notes, 'CONTAINS_INF')));
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

fprintf(fid, '--- DURATION (computed separately PER TASK TYPE -- see header item 5) ---\n');
for t = 1:numel(EXPECTED_TASK_TYPES)
    tt = EXPECTED_TASK_TYPES{t};
    idxThisTask = strcmp(T.TaskType, tt) & strcmp(T.LoadStatus, 'OK');
    dThis = T.DurationSec(idxThisTask);
    if isempty(dThis)
        fprintf(fid, '  %s: no successfully loaded files.\n', tt);
        continue;
    end
    fprintf(fid, '  %s: median=%.2fs min=%.2fs max=%.2fs n=%d outliers(>±%.0f%% of that median)=%d\n', ...
        tt, med_dur_by_task.(tt), min(dThis), max(dThis), numel(dThis), DURATION_OUTLIER_FRACTION*100, ...
        sum(contains(T.Notes, 'DURATION_OUTLIER_WITHIN_TASKTYPE') & strcmp(T.TaskType, tt)));
end
fprintf(fid, '\n');

fprintf(fid, '--- CONDITION COUNTS (PROVISIONAL / filename-derived; includes Baseline) ---\n');
condsFound = T.ConditionCandidate(strcmp(T.LoadStatus,'OK'));
for c = {'Baseline','E','N','R'}
    fprintf(fid, '  %s : %d recording(s)\n', c{1}, sum(strcmp(condsFound, c{1})));
end
fprintf(fid, '\n');

fprintf(fid, '--- TASK TYPE / BLOCK COUNTS ---\n');
for t = 1:numel(EXPECTED_TASK_TYPES)
    tt = EXPECTED_TASK_TYPES{t};
    for b = 1:EXPECTED_BLOCKS_PER_SUBJECT
        fprintf(fid, '  %s, Block %d : %d recording(s)\n', tt, b, sum(strcmp(T.TaskType, tt) & T.BlockNumber == b));
    end
end
fprintf(fid, '\n');

fprintf(fid, '--- BLOCK / CONDITION STRUCTURAL CONSISTENCY (Block 1 should be Baseline; Blocks 2/3 should not be) ---\n');
inconsistIdx = find(~T.BlockConditionConsistent & T.ParseOK);
if isempty(inconsistIdx)
    fprintf(fid, '  None -- every parsed file''s block number and condition letter are structurally consistent.\n');
else
    for i = inconsistIdx'
        fprintf(fid, '  %s : block=%d cond=%s (unexpected combination)\n', T.FileName{i}, T.BlockNumber(i), T.ConditionCandidate{i});
    end
end
fprintf(fid, '\n');

fprintf(fid, '--- MISSING Subject x TaskType x Block COMBINATIONS ---\n');
if isempty(MissingT)
    fprintf(fid, '  None.\n');
else
    for i = 1:height(MissingT)
        fprintf(fid, '  Sub%s, %s, Block %d\n', MissingT.Subject{i}, MissingT.TaskType{i}, MissingT.MissingBlockNumber(i));
    end
end
fprintf(fid, '\n');

fprintf(fid, '--- DUPLICATE RECORDINGS ---\n');
if isempty(DupT)
    fprintf(fid, '  None.\n');
else
    for i = 1:height(DupT)
        fprintf(fid, '  Sub%s, %s, Block %d : %s\n', DupT.Subject{i}, DupT.TaskType{i}, DupT.BlockNumber(i), DupT.FileName{i});
    end
end
fprintf(fid, '\n');

fprintf(fid, '--- CONDITION MISMATCHES (filename-derived vs. allocation-table) ---\n');
if ~ALLOC_CROSS_VALIDATION_TRUSTED
    fprintf(fid, '  N/A -- cross-validation disabled (see status note above).\n');
elseif isempty(MismatchT)
    fprintf(fid, '  None.\n');
else
    for i = 1:height(MismatchT)
        fprintf(fid, '  Sub%s, Block%d, file=%s : raw="%s" -> candidate=%s, but allocation table expects %s\n', ...
            MismatchT.Subject{i}, MismatchT.BlockNumber(i), MismatchT.FileName{i}, ...
            MismatchT.ConditionLetterRaw{i}, MismatchT.ConditionCandidate{i}, MismatchT.ConditionExpectedFromAllocation{i});
    end
end
fprintf(fid, '\n');

fprintf(fid, '--- FILES NEEDING CONDITION-LETTER CONFIRMATION ---\n');
confirmIdx = find(contains(T.Notes, 'CONDITION_LETTER_NEEDS_CONFIRMATION'));
if isempty(confirmIdx)
    fprintf(fid, '  None.\n');
else
    for i = confirmIdx'
        fprintf(fid, '  %s : raw="%s" -> assumed %s\n', T.FileName{i}, T.ConditionLetterRaw{i}, T.ConditionCandidate{i});
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
    fprintf(fid, '  N/A -- no file loaded successfully.\n');
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

fprintf(fid, '--- VI.set FILES (out of scope, logged only, never loaded) ---\n');
fprintf(fid, '  %d VI.set file(s) found.\n', numel(vi_list));
for i = 1:height(ViT)
    fprintf(fid, '  %s\n', ViT.FullPath{i});
end

fclose(fid);

fprintf('\n========================================================\n');
fprintf('AUDIT COMPLETE.\n');
fprintf('Master table : %s\n', fullfile(AUDIT_DIR, 'audit_table_cognitive.csv'));
fprintf('Report       : %s\n', rpt);
fprintf('========================================================\n');
fprintf('\nSTOP: review the audit table and report before running X02.\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%  (must appear here, after all top-level script code, per MATLAB rules)
% ========================================================================

function [subj, seqRaw, condRaw, blockNum, taskType, ok] = parse_cognitive_filename(fname)
    % Expected convention:
    %   sub<N>_<SEQCOND(3 letters)>_<BLOCK(1-3)>_<full|nback|sart>_pp.set
    % Fully case-insensitive for the "sub" prefix and task-type token;
    % the 3-letter SEQCOND code's own case is preserved verbatim in
    % seqRaw/condRaw for downstream reporting (case is not semantically
    % meaningful for E/N/R/B/I per map_condition_letter below, but is
    % preserved rather than silently normalized here).
    FNAME_PATTERN = '^[Ss]ub0*(\d+)_([A-Za-z]{3})_([1-3])_([A-Za-z]+)_pp$';
    [~, base, ~] = fileparts(fname);
    tok = regexp(base, FNAME_PATTERN, 'tokens', 'once');
    if isempty(tok)
        subj = ''; seqRaw = ''; condRaw = ''; blockNum = NaN; taskType = ''; ok = false;
        return;
    end
    subj = tok{1};                     % CHAR, e.g. '16'
    code3 = tok{2};                    % e.g. 'NSB'
    blockNum = str2double(tok{3});     % 1, 2, or 3
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
        taskType = taskTypeRaw;   % preserved for reporting even though unrecognized
        return;
    end
    taskType = taskTypeRaw;
    ok = true;
end

function [condCandidate, needsConfirm] = map_condition_letter(condRaw)
    % Maps the 3rd letter of the SEQCOND code to a candidate condition
    % label. IDENTICAL E/I/N/R mapping to the intervention pipeline (I->E),
    % plus a B->Baseline mapping for the cognitive-task convention's
    % pre-intervention block.
    c = upper(condRaw);
    switch c
        case 'B'
            condCandidate = 'Baseline'; needsConfirm = false;
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
    if cond
        out = a;
    else
        out = b;
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
