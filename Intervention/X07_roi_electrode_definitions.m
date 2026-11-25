%% ========================================================================
% X07_roi_electrode_definitions.m
%
% Channel table (labels and XYZ coordinates) and the five regions of interest used for
% the regional analyses: frontal (23), central (7), temporal (6), parietal (16), occipital (8);
% P9 and P10 are not assigned to a region.
% Inputs : 02_metadata, 03_psd, channel_locations_template.mat (channel coordinates).
% Outputs: 07_roi_electrode_definitions/phase5_channel_roi_definitions.mat/.csv/.xlsx, 09_reports/.
% ========================================================================

clear; clc;

%% ------------------------------------------------------------------

% >>> Must match the ANALYSIS_ROOT used in Phases 1-4 <<<
config;
ANALYSIS_ROOT = OUT_DIR;

META_DIR = fullfile(ANALYSIS_ROOT, '02_metadata');
PSD_DIR  = fullfile(ANALYSIS_ROOT, '03_psd');
DEF_DIR  = fullfile(ANALYSIS_ROOT, '07_roi_electrode_definitions');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');

META_MAT_FILE = fullfile(META_DIR, 'intervention_metadata.mat');
PSD_MAT_FILE  = fullfile(PSD_DIR, 'intervention_psd.mat');

% Channel coordinates: channel_locations_template.mat holds EEG.chanlocs
% (variables Chanlocs, Labels, ref_file) of one representative preprocessed
% recording and is expected in DATA_DIR.
CHANLOCS_TEMPLATE_FILE = fullfile(DATA_DIR, 'channel_locations_template.mat');

if ~exist(DEF_DIR, 'dir');    mkdir(DEF_DIR);    end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

fprintf('========================================================\n');
fprintf('PHASE 5D -- ROI / ELECTRODE DEFINITION AND VALIDATION\n');
fprintf('========================================================\n');
fprintf('ANALYSIS_ROOT           : %s\n', ANALYSIS_ROOT);
fprintf('META_MAT_FILE            : %s\n', META_MAT_FILE);
fprintf('PSD_MAT_FILE             : %s\n', PSD_MAT_FILE);
fprintf('CHANLOCS_TEMPLATE_FILE   : %s\n', CHANLOCS_TEMPLATE_FILE);
fprintf('OUTPUT                   : %s\n', DEF_DIR);
fprintf('--------------------------------------------------------\n\n');

if ~exist(META_MAT_FILE, 'file')
    error(['Phase 1 metadata file not found: %s\nRun ' ...
           '02_build_intervention_metadata.m first.'], META_MAT_FILE);
end
if ~exist(PSD_MAT_FILE, 'file')
    error(['Phase 2 PSD file not found: %s\nRun ' ...
           '03_compute_intervention_psd.m first.'], PSD_MAT_FILE);
end
if ~exist(CHANLOCS_TEMPLATE_FILE, 'file')
    error(['Channel-location template not found: %s\nEdit ' ...
           'CHANLOCS_TEMPLATE_FILE above to point at your copy of ' ...
           'channel_locations_template.mat. Real XYZ coordinates are ' ...
           'required -- this script will not fabricate or substitute ' ...
           'a generic montage.'], CHANLOCS_TEMPLATE_FILE);
end

%% ------------------------------------------------------------------
%  1. LOAD THE AUTHORITATIVE 62-CHANNEL REFERENCE ORDER
%     (cross-checked between its two independent save points rather
%     than trusted from either one alone)
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
        ['Phase 1''s ref_labels and Phase 2''s ChannelLabelsReference are not ' ...
         'identical. This is exactly the silent-divergence case this project ' ...
         'never allows to pass unnoticed. Stopping -- investigate before ' ...
         'proceeding (compare intervention_metadata.mat and intervention_psd.mat).']);
end

ChannelLabelsReference = ref_labels;   % identical; ref_labels used as the single working copy from here on
nRefChan = numel(ChannelLabelsReference);
fprintf('Loaded and cross-checked %d-channel reference order (Phase 1 ref_labels == Phase 2 ChannelLabelsReference).\n\n', nRefChan);

%% ------------------------------------------------------------------
%  2. LOAD REAL CHANNEL XYZ COORDINATES AND MATCH TO THE REFERENCE ORDER
%     (matched by LABEL, case-insensitively -- never assumed positional)
% ------------------------------------------------------------------

locFile = load(CHANLOCS_TEMPLATE_FILE, 'Chanlocs', 'Labels', 'ref_file');
Chanlocs = locFile.Chanlocs;
locLabelsUC = upper(locFile.Labels(:));

X = nan(nRefChan, 1); Y = nan(nRefChan, 1); Z = nan(nRefChan, 1);
Theta = nan(nRefChan, 1); Radius = nan(nRefChan, 1);
CoordLabelMismatch = cell(0, 1);   % reference channels with no coordinate match

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

% Coordinate-template entries that never matched ANY reference channel
% (would indicate the template has extra/misnamed channels -- report,
% don't silently ignore).
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

% Excluded from the 5-ROI averages only; still analyzed individually at
% the electrode level (09_electrode_level_statistics.m must NOT drop
% these two).
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

fprintf('\nDuplicate assignments (channel listed in more than one ROI, or in both\n');
fprintf('a ROI and the unassigned list): ');
if isempty(DuplicateAssignments)
    fprintf('NONE\n');
else
    fprintf('%d -- %s\n', numel(DuplicateAssignments), strjoin(DuplicateAssignments, ', '));
    fails = fails + 1;
end

fprintf('\nMissing channels (in the 62-channel reference order but not assigned\n');
fprintf('to any ROI and not in the unassigned list): ');
if isempty(MissingChannels)
    fprintf('NONE\n');
else
    fprintf('%d -- %s\n', numel(MissingChannels), strjoin(MissingChannels, ', '));
    fails = fails + 1;
end

fprintf('\nUnknown labels in the ROI/unassigned definition (listed in ROI_DEFINITION or\n');
fprintf('UNASSIGNED_FOR_ROI but not found among the 62 reference channels -- would\n');
fprintf('indicate a typo in this script): ');
if isempty(UnknownLabelsInDefinition)
    fprintf('NONE\n');
else
    fprintf('%d -- %s\n', numel(UnknownLabelsInDefinition), strjoin(UnknownLabelsInDefinition, ', '));
    fails = fails + 1;
end

fprintf('\nCoordinate-label mismatches (reference channels with no XYZ match in\n');
fprintf('the channel-location template): ');
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
PARAMS5D.ScriptName = '07_roi_electrode_definitions.m';
PARAMS5D.GeneratedOn = datestr(now);
PARAMS5D.SourceMetaFile = META_MAT_FILE;
PARAMS5D.SourcePsdFile = PSD_MAT_FILE;
PARAMS5D.SourceChanlocsFile = CHANLOCS_TEMPLATE_FILE;
PARAMS5D.ChanlocsRefRecording = locFile.ref_file;
PARAMS5D.ROIDecision = 'Option B (practical 5-ROI resolution), user-confirmed after reviewing actual coordinate export';
PARAMS5D.UnassignedForROI = {UNASSIGNED_FOR_ROI};

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
reportLines{end+1} = 'PHASE 5D -- ROI / ELECTRODE DEFINITION AND VALIDATION REPORT';
reportLines{end+1} = repmat('=', 1, 60);
reportLines{end+1} = sprintf('Generated: %s', datestr(now));
reportLines{end+1} = '';
reportLines{end+1} = sprintf('Reference channel order source : %s (ref_labels)', META_MAT_FILE);
reportLines{end+1} = sprintf('                                  cross-checked against ChannelLabelsReference in %s', PSD_MAT_FILE);
reportLines{end+1} = sprintf('Coordinate source               : %s', CHANLOCS_TEMPLATE_FILE);
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
reportLines{end+1} = 'ROI membership above is fixed by explicit user sign-off and is never';
reportLines{end+1} = 'derived from, or adjusted in response to, any EEG result.';

fid = fopen(fullfile(REPORT_DIR, 'phase5_roi_electrode_definitions_report.txt'), 'w');
fprintf(fid, '%s\n', reportLines{:});
fclose(fid);

fprintf('Saved: %s\n', fullfile(REPORT_DIR, 'phase5_roi_electrode_definitions_report.txt'));
fprintf('\n=== 07_roi_electrode_definitions.m complete. STOP for review before writing 08/09/10. ===\n');

%% ========================================================================
%  LOCAL FUNCTIONS
%% ========================================================================

function [RoiLabelPerChannel, DuplicateAssignments, MissingChannels, UnknownLabelsInDefinition] = ...
    local_validate_roi_membership(refLabels, roiDef, roiNames, unassignedList)

    n = numel(refLabels);
    refUC = upper(refLabels);
    RoiLabelPerChannel = repmat({'UNASSIGNED'}, n, 1);
    hitCount = zeros(n, 1);   % how many places (ROIs + unassigned) claim this channel
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
