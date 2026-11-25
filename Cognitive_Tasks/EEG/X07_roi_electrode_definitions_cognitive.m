%% ========================================================================
% X07_roi_electrode_definitions_cognitive.m
%
% Purpose: defines the 62-electrode montage and the 5-region ROI
% membership table ("Option B") shared by every downstream script
% (X08/X09/X10) for every Task. This script is not run per Task -- the
% montage and ROI definition are properties of the recording
% hardware/montage, not of Task, so exactly one shared definition is
% produced here. ROI/electrode definitions themselves are not re-derived
% here; the same definitions are used throughout ("same ROI definitions,
% no new ROIs, same 62 electrodes").
%
% Method: loads the reference channel order from X02's
% cognitive_metadata.mat (ref_labels), cross-checks it against X03's
% cognitive_psd.mat (ChannelLabelsReference), and builds the 5-ROI
% membership table (P9/P10 are excluded from ROI averaging only) by
% matching real channel XYZ coordinates (from the channel-location
% template) to each channel's label -- never fabricated or substituted
% with a generic montage.
%
% Requirements: base MATLAB only (no EEGLAB, no toolboxes). X02 output
% must exist (cognitive_metadata.mat, ref_labels); X03 output must exist
% (cognitive_psd.mat, ChannelLabelsReference); channel_locations_template.mat
% must exist at CHANLOCS_TEMPLATE_FILE below (or its fallback candidate).
%
% Inputs: DATA_DIR (for the channel-location template), X02 metadata, X03 PSD output.
%
% Outputs (shared across all three Tasks):
%   <ANALYSIS_ROOT>/07_roi_electrode_definitions/
%       phase5_channel_roi_definitions.mat   (ChannelDefT + ROI_DEFINITION + PARAMS5D)
%       phase5_channel_roi_definitions.xlsx
%       phase5_channel_roi_definitions.csv
%   <ANALYSIS_ROOT>/09_reports/
%       phase5_roi_electrode_definitions_report_cognitive.txt
%
% RUN_MODE: not applicable; this script runs once.
% ========================================================================

clear; clc;

%% ------------------------------------------------------------------
%  0. USER-EDITABLE SETTINGS
% ------------------------------------------------------------------

config;
DATA_ROOT_FOR_SIBLING = DATA_DIR;
[DATA_ROOT_PARENT, ~] = fileparts(DATA_ROOT_FOR_SIBLING);
ANALYSIS_ROOT = OUT_DIR;

META_DIR = fullfile(ANALYSIS_ROOT, '02_metadata');
PSD_DIR  = fullfile(ANALYSIS_ROOT, '03_psd');
DEF_DIR  = fullfile(ANALYSIS_ROOT, '07_roi_electrode_definitions');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');

META_MAT_FILE = fullfile(META_DIR, 'cognitive_metadata.mat');
PSD_MAT_FILE  = fullfile(PSD_DIR, 'cognitive_psd.mat');

% Primary candidate: the channel-location template shipped alongside the
% preprocessed data (DATA_DIR). Fallback candidate: the intervention
% pipeline's fileparts(fileparts(...))-derived formula, in case the same
% physical file is reachable that way instead. Whichever resolves is
% logged explicitly; if neither resolves, this script stops with a hard
% error rather than fabricating/substituting a generic montage.
CHANLOCS_CANDIDATE_1 = fullfile(DATA_DIR, 'channel_locations_template.mat');
CHANLOCS_CANDIDATE_2 = fullfile(fileparts(fileparts(ANALYSIS_ROOT)), 'Code', 'results', 'metadata', 'channel_locations_template.mat');

if exist(CHANLOCS_CANDIDATE_1, 'file')
    CHANLOCS_TEMPLATE_FILE = CHANLOCS_CANDIDATE_1;
    CHANLOCS_SOURCE_NOTE = 'candidate 1 (DATA_DIR path)';
elseif exist(CHANLOCS_CANDIDATE_2, 'file')
    CHANLOCS_TEMPLATE_FILE = CHANLOCS_CANDIDATE_2;
    CHANLOCS_SOURCE_NOTE = 'candidate 2 (intervention-pipeline-derived fallback path)';
else
    CHANLOCS_TEMPLATE_FILE = '';
    CHANLOCS_SOURCE_NOTE = 'NEITHER CANDIDATE RESOLVED';
end

if ~exist(DEF_DIR, 'dir');    mkdir(DEF_DIR);    end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

fprintf('========================================================\n');
fprintf('PHASE 5D -- ROI / ELECTRODE DEFINITION AND VALIDATION (COGNITIVE-TASK DATA)\n');
fprintf('========================================================\n');
fprintf('ANALYSIS_ROOT           : %s\n', ANALYSIS_ROOT);
fprintf('META_MAT_FILE            : %s\n', META_MAT_FILE);
fprintf('PSD_MAT_FILE             : %s\n', PSD_MAT_FILE);
fprintf('CHANLOCS_TEMPLATE_FILE   : %s  [%s]\n', CHANLOCS_TEMPLATE_FILE, CHANLOCS_SOURCE_NOTE);
fprintf('OUTPUT                   : %s\n', DEF_DIR);
fprintf('--------------------------------------------------------\n\n');

if ~exist(META_MAT_FILE, 'file')
    error(['X02 metadata file not found: %s\nRun ' ...
           'X02_build_cognitive_metadata.m first.'], META_MAT_FILE);
end
if ~exist(PSD_MAT_FILE, 'file')
    error(['X03 PSD file not found: %s\nRun ' ...
           'X03_compute_cognitive_psd.m first.'], PSD_MAT_FILE);
end
if isempty(CHANLOCS_TEMPLATE_FILE) || ~exist(CHANLOCS_TEMPLATE_FILE, 'file')
    error(['Channel-location template not found at either candidate path:\n  %s\n  %s\n' ...
           'Edit CHANLOCS_CANDIDATE_1/2 above to point at your copy of ' ...
           'channel_locations_template.mat. Real XYZ coordinates are ' ...
           'required -- this script will not fabricate or substitute ' ...
           'a generic montage.'], CHANLOCS_CANDIDATE_1, CHANLOCS_CANDIDATE_2);
end

%% ------------------------------------------------------------------
%  1. LOAD THE AUTHORITATIVE 62-CHANNEL REFERENCE ORDER
%     (cross-checked between its two independent save points)
% ------------------------------------------------------------------

metaFile = load(META_MAT_FILE, 'ref_labels');
ref_labels = metaFile.ref_labels(:);

psdFile = load(PSD_MAT_FILE, 'ChannelLabelsReference');
ChannelLabelsReference = psdFile.ChannelLabelsReference(:);

if numel(ref_labels) ~= 62
    error('Phase5D:UnexpectedChannelCount', ...
        'ref_labels has %d entries; this project is defined for a 62-channel montage. Stopping rather than proceeding on an unexpected count.', ...
        numel(ref_labels));
end
if numel(ChannelLabelsReference) ~= numel(ref_labels) || ...
        ~isequal(ChannelLabelsReference, ref_labels)
    error('Phase5D:ReferenceOrderMismatch', ...
        ['X02''s ref_labels and X03''s ChannelLabelsReference are not ' ...
         'identical. This is exactly the silent-divergence case this project ' ...
         'never allows to pass unnoticed. Stopping -- investigate before ' ...
         'proceeding (compare cognitive_metadata.mat and cognitive_psd.mat).']);
end

ChannelLabelsReference = ref_labels;
nRefChan = numel(ChannelLabelsReference);
fprintf('Loaded and cross-checked %d-channel reference order (X02 ref_labels == X03 ChannelLabelsReference).\n\n', nRefChan);

%% ------------------------------------------------------------------
%  2. LOAD REAL CHANNEL XYZ COORDINATES AND MATCH TO THE REFERENCE ORDER
%     (matched by LABEL, case-insensitively -- never assumed positional)
% ------------------------------------------------------------------

locFile = load(CHANLOCS_TEMPLATE_FILE, 'Chanlocs', 'Labels', 'ref_file');
Chanlocs = locFile.Chanlocs;
locLabelsUC = upper(locFile.Labels(:));

X = nan(nRefChan, 1); Y = nan(nRefChan, 1); Z = nan(nRefChan, 1);
Theta = nan(nRefChan, 1); Radius = nan(nRefChan, 1);
CoordLabelMismatch = cell(0, 1);

for i = 1:nRefChan
    target = upper(ChannelLabelsReference{i});
    hit = find(strcmp(locLabelsUC, target));
    if numel(hit) ~= 1
        CoordLabelMismatch{end+1, 1} = ChannelLabelsReference{i}; %#ok<AGROW>
        continue;
    end
    X(i) = Chanlocs(hit).X;
    Y(i) = Chanlocs(hit).Y;
    Z(i) = Chanlocs(hit).Z;
    Theta(i) = Chanlocs(hit).theta;
    Radius(i) = Chanlocs(hit).radius;
end

usedTemplateLabels = false(numel(locLabelsUC), 1);
for i = 1:nRefChan
    hit = find(strcmp(locLabelsUC, upper(ChannelLabelsReference{i})));
    if numel(hit) == 1
        usedTemplateLabels(hit) = true;
    end
end
UnusedTemplateLabels = locFile.Labels(~usedTemplateLabels);

fprintf('Channel-location template: %s\n', locFile.ref_file);
fprintf('Matched %d / %d reference channels to real XYZ coordinates.\n', ...
    nRefChan - numel(CoordLabelMismatch), nRefChan);
if ~isempty(CoordLabelMismatch)
    fprintf('  *** %d reference channel(s) have NO coordinate match: %s\n', ...
        numel(CoordLabelMismatch), strjoin(CoordLabelMismatch, ', '));
end
if ~isempty(UnusedTemplateLabels)
    fprintf('  NOTE: %d template label(s) did not match any reference channel: %s\n', ...
        numel(UnusedTemplateLabels), strjoin(UnusedTemplateLabels, ', '));
end
fprintf('\n');

%% ------------------------------------------------------------------
%  3. ROI MEMBERSHIP -- FIXED CONSTANTS, VERBATIM FROM 07_roi_electrode_
%     definitions.m (see header: "ROI DEFINITION -- LOCKED"). NEVER
%     derive, adjust, or overwrite these lists from any EEG result.
% ------------------------------------------------------------------

ROI_DEFINITION = struct();
ROI_DEFINITION.Frontal  = {'AFZ','FP1','FP2','AF3','AF4','F7','F3','FZ','F4','F8', ...
                            'FC5','FC1','FC2','FC6','AF7','AF8','F5','F1','F2','F6', ...
                            'FC3','FCz','FC4'};
ROI_DEFINITION.Central  = {'C5','C3','C1','Cz','C2','C4','C6'};
ROI_DEFINITION.Temporal = {'T7','T8','FT7','FT8','TP7','TP8'};
ROI_DEFINITION.Parietal = {'CP5','CP1','CP2','CP6','P7','P3','PZ','P4','P8', ...
                            'CP3','CPz','CP4','P5','P1','P2','P6'};
ROI_DEFINITION.Occipital = {'PO3','PO4','O1','O2','PO7','POz','PO8','Oz'};

ROI_NAMES = {'Frontal','Central','Temporal','Parietal','Occipital'};

UNASSIGNED_FOR_ROI = {'P9','P10'};

EXPECTED_ROI_COUNTS = containers.Map( ...
    {'Frontal','Central','Temporal','Parietal','Occipital'}, ...
    {23, 7, 6, 16, 8});

%% ------------------------------------------------------------------
%  4. VALIDATION -- every one of the 62 reference channels must be
%     accounted for EXACTLY ONCE (one ROI, or explicitly unassigned)
% ------------------------------------------------------------------

[RoiLabelPerChannel, DuplicateAssignments, MissingChannels, UnknownLabelsInDefinition] = ...
    local_validate_roi_membership(ChannelLabelsReference, ROI_DEFINITION, ROI_NAMES, UNASSIGNED_FOR_ROI);

fails = 0;

fprintf('--------------------------------------------------------\n');
fprintf('VALIDATION\n');
fprintf('--------------------------------------------------------\n');

fprintf('Channels per ROI (n=%d total assigned to a ROI):\n', sum(~strcmp(RoiLabelPerChannel, 'UNASSIGNED')));
for r = 1:numel(ROI_NAMES)
    n = sum(strcmp(RoiLabelPerChannel, ROI_NAMES{r}));
    expected = EXPECTED_ROI_COUNTS(ROI_NAMES{r});
    match = '';
    if n ~= expected
        match = sprintf('   *** MISMATCH: expected %d ***', expected);
        fails = fails + 1;
    end
    fprintf('  %-10s : %3d%s\n', ROI_NAMES{r}, n, match);
end
nUnassigned = sum(strcmp(RoiLabelPerChannel, 'UNASSIGNED'));
fprintf('  %-10s : %3d  (%s)\n', 'Unassigned', nUnassigned, strjoin(UNASSIGNED_FOR_ROI, ', '));
if nUnassigned ~= numel(UNASSIGNED_FOR_ROI)
    fprintf('  *** MISMATCH: expected %d unassigned channel(s) ***\n', numel(UNASSIGNED_FOR_ROI));
    fails = fails + 1;
end

fprintf('\nDuplicate assignments: ');
if isempty(DuplicateAssignments)
    fprintf('NONE\n');
else
    fprintf('%d -- %s\n', numel(DuplicateAssignments), strjoin(DuplicateAssignments, ', '));
    fails = fails + 1;
end

fprintf('\nMissing channels: ');
if isempty(MissingChannels)
    fprintf('NONE\n');
else
    fprintf('%d -- %s\n', numel(MissingChannels), strjoin(MissingChannels, ', '));
    fails = fails + 1;
end

fprintf('\nUnknown labels in the ROI/unassigned definition: ');
if isempty(UnknownLabelsInDefinition)
    fprintf('NONE\n');
else
    fprintf('%d -- %s\n', numel(UnknownLabelsInDefinition), strjoin(UnknownLabelsInDefinition, ', '));
    fails = fails + 1;
end

fprintf('\nCoordinate-label mismatches: ');
if isempty(CoordLabelMismatch)
    fprintf('NONE\n');
else
    fprintf('%d -- %s\n', numel(CoordLabelMismatch), strjoin(CoordLabelMismatch, ', '));
    fails = fails + 1;
end

fprintf('\n=== TOTAL VALIDATION FAILURES: %d ===\n\n', fails);
if fails > 0
    error('Phase5D:ValidationFailed', ...
        'ROI/electrode definition validation failed (%d issue(s) above). Fix before saving output.', fails);
end

%% ------------------------------------------------------------------
%  5. BUILD AND SAVE THE CHANNEL DEFINITION TABLE
% ------------------------------------------------------------------

ChannelDefT = table((1:nRefChan)', ChannelLabelsReference, RoiLabelPerChannel, ...
    X, Y, Z, Theta, Radius, ...
    'VariableNames', {'Index','Label','ROI','X','Y','Z','Theta','Radius'});

PARAMS5D = struct();
PARAMS5D.ScriptName = 'X07_roi_electrode_definitions_cognitive.m';
PARAMS5D.GeneratedOn = datestr(now);
PARAMS5D.SourceMetaFile = META_MAT_FILE;
PARAMS5D.SourcePsdFile = PSD_MAT_FILE;
PARAMS5D.SourceChanlocsFile = CHANLOCS_TEMPLATE_FILE;
PARAMS5D.SourceChanlocsFileNote = CHANLOCS_SOURCE_NOTE;
PARAMS5D.ChanlocsRefRecording = locFile.ref_file;
PARAMS5D.ROIDecision = 'Option B (practical 5-ROI resolution) -- identical to the intervention pipeline, never re-derived for this dataset.';
PARAMS5D.UnassignedForROI = {UNASSIGNED_FOR_ROI};
PARAMS5D.SharedAcrossTasks = true;

save(fullfile(DEF_DIR, 'phase5_channel_roi_definitions.mat'), ...
    'ChannelDefT', 'ROI_DEFINITION', 'ROI_NAMES', 'UNASSIGNED_FOR_ROI', 'PARAMS5D', '-v7.3');
writetable(ChannelDefT, fullfile(DEF_DIR, 'phase5_channel_roi_definitions.xlsx'));
writetable(ChannelDefT, fullfile(DEF_DIR, 'phase5_channel_roi_definitions.csv'));

fprintf('Saved: %s\n', fullfile(DEF_DIR, 'phase5_channel_roi_definitions.mat'));
fprintf('Saved: %s\n', fullfile(DEF_DIR, 'phase5_channel_roi_definitions.xlsx'));
fprintf('Saved: %s\n\n', fullfile(DEF_DIR, 'phase5_channel_roi_definitions.csv'));

%% ------------------------------------------------------------------
%  6. WRITE REPORT
% ------------------------------------------------------------------

reportLines = {};
reportLines{end+1} = 'PHASE 5D -- ROI / ELECTRODE DEFINITION AND VALIDATION REPORT (COGNITIVE-TASK DATA)';
reportLines{end+1} = repmat('=', 1, 60);
reportLines{end+1} = sprintf('Generated: %s', datestr(now));
reportLines{end+1} = '';
reportLines{end+1} = 'This definition is SHARED across all three Tasks (Full/Nback/SART) -- the';
reportLines{end+1} = '62-electrode montage and 5-ROI membership do not depend on Task.';
reportLines{end+1} = '';
reportLines{end+1} = sprintf('Reference channel order source : %s (ref_labels)', META_MAT_FILE);
reportLines{end+1} = sprintf('                                  cross-checked against ChannelLabelsReference in %s', PSD_MAT_FILE);
reportLines{end+1} = sprintf('Coordinate source               : %s [%s]', CHANLOCS_TEMPLATE_FILE, CHANLOCS_SOURCE_NOTE);
reportLines{end+1} = sprintf('  (real EEG.chanlocs read from one representative recording: %s)', locFile.ref_file);
reportLines{end+1} = '';
reportLines{end+1} = 'CHANNELS PER ROI:';
for r = 1:numel(ROI_NAMES)
    n = sum(strcmp(RoiLabelPerChannel, ROI_NAMES{r}));
    reportLines{end+1} = sprintf('  %-10s : %d', ROI_NAMES{r}, n); %#ok<AGROW>
end
reportLines{end+1} = sprintf('  %-10s : %d (%s)', 'Unassigned', nUnassigned, strjoin(UNASSIGNED_FOR_ROI, ', '));
reportLines{end+1} = '';
reportLines{end+1} = sprintf('Duplicate assignments      : %d', numel(DuplicateAssignments));
reportLines{end+1} = sprintf('Missing channels            : %d', numel(MissingChannels));
reportLines{end+1} = sprintf('Unknown labels in definition : %d', numel(UnknownLabelsInDefinition));
reportLines{end+1} = sprintf('Coordinate-label mismatches : %d', numel(CoordLabelMismatch));
reportLines{end+1} = '';
reportLines{end+1} = sprintf('TOTAL VALIDATION FAILURES: %d', fails);
reportLines{end+1} = '';
reportLines{end+1} = 'This script performs NO inferential statistics. No p-values, FDR,';
reportLines{end+1} = 'Holm correction, effect sizes, or E/N/R comparisons are computed here.';
reportLines{end+1} = 'ROI membership above is fixed by the same explicit user sign-off as the';
reportLines{end+1} = 'intervention pipeline and is never derived from, or adjusted in response';
reportLines{end+1} = 'to, any EEG result.';

fid = fopen(fullfile(REPORT_DIR, 'phase5_roi_electrode_definitions_report_cognitive.txt'), 'w');
fprintf(fid, '%s\n', reportLines{:});
fclose(fid);

fprintf('Saved: %s\n', fullfile(REPORT_DIR, 'phase5_roi_electrode_definitions_report_cognitive.txt'));
fprintf('\n=== X07 complete. STOP for review before writing/running X08/X09/X10. ===\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%% ========================================================================

function [RoiLabelPerChannel, DuplicateAssignments, MissingChannels, UnknownLabelsInDefinition] = ...
    local_validate_roi_membership(refLabels, roiDef, roiNames, unassignedList)
    n = numel(refLabels);
    refUC = upper(refLabels);
    RoiLabelPerChannel = repmat({'UNASSIGNED'}, n, 1);
    hitCount = zeros(n, 1);
    UnknownLabelsInDefinition = {};

    for r = 1:numel(roiNames)
        theseLabels = roiDef.(roiNames{r});
        for k = 1:numel(theseLabels)
            idx = find(strcmp(refUC, upper(theseLabels{k})));
            if isempty(idx)
                UnknownLabelsInDefinition{end+1} = sprintf('%s (ROI=%s)', theseLabels{k}, roiNames{r}); %#ok<AGROW>
                continue;
            end
            if hitCount(idx) == 0
                RoiLabelPerChannel{idx} = roiNames{r};
            end
            hitCount(idx) = hitCount(idx) + 1;
        end
    end

    for k = 1:numel(unassignedList)
        idx = find(strcmp(refUC, upper(unassignedList{k})));
        if isempty(idx)
            UnknownLabelsInDefinition{end+1} = sprintf('%s (unassigned list)', unassignedList{k}); %#ok<AGROW>
            continue;
        end
        if hitCount(idx) == 0
            RoiLabelPerChannel{idx} = 'UNASSIGNED';
        end
        hitCount(idx) = hitCount(idx) + 1;
    end

    DuplicateAssignments = refLabels(hitCount > 1);
    MissingChannels = refLabels(hitCount == 0);
end
