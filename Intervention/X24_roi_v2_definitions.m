%% ========================================================================
% X24_roi_v2_definitions.m
%
% Defines the four literature-based regions used for the planned contrasts and for the
% regional band-power results during the interventions: frontal midline (Fz, FCz, F1, F2,
% FC1, FC2), central (C3, C1, Cz, C2, C4), parietal (P3, P1, Pz, P2, P4) and occipital
% (PO3, POz, PO4, O1, Oz, O2).
% Inputs : 07_roi_electrode_definitions.
% Outputs: 24_roi_v2_definitions/roi_v2_definitions.csv, 09_reports/.
% ========================================================================

clear; clc;

%% 0. SETTINGS
config;
ANALYSIS_ROOT = OUT_DIR;
DEF_DIR    = fullfile(ANALYSIS_ROOT, '07_roi_electrode_definitions');
OUT_DIR    = fullfile(ANALYSIS_ROOT, '24_roi_v2_definitions');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');
roiDefFile = fullfile(DEF_DIR, 'phase5_channel_roi_definitions.mat');

if ~exist(roiDefFile, 'file')
    error('Phase6A:MissingInput', 'Phase 5D file not found: %s\nRun 07_roi_electrode_definitions.m first.', roiDefFile);
end
if ~exist(OUT_DIR, 'dir');    mkdir(OUT_DIR);    end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

%% 1. ROI DEFINITIONS (fixed constants)
ROI2_DEFINITION = struct();
ROI2_DEFINITION.FrontalMidline = {'FZ','FCz','F1','F2','FC1','FC2'};
ROI2_DEFINITION.Central        = {'C3','C1','Cz','C2','C4'};
ROI2_DEFINITION.Parietal       = {'P3','P1','PZ','P2','P4'};
ROI2_DEFINITION.Occipital      = {'PO3','POz','PO4','O1','Oz','O2'};
ROI2_NAMES = {'FrontalMidline','Central','Parietal','Occipital'};

ROI2_RATIONALE = struct( ...
    'FrontalMidline', 'Frontal-midline theta (attention / cognitive control): Ishihara & Yoshii 1972; Mitchell et al. 2008; Cavanagh & Frank 2014', ...
    'Central',        'Sensorimotor mu/beta (finger movement while scrolling): Pfurtscheller & Lopes da Silva 1999', ...
    'Parietal',       'Parietal alpha (attention allocation): Klimesch 1999, 2012', ...
    'Occipital',      'Posterior alpha (visual processing / engagement): Berger 1929; Klimesch 1999');

EDGE_EXCLUDED = {'FP1','FP2','AF7','AF8','F7','F8','FT7','FT8','T7','T8','TP7','TP8','P9','P10'};

% Confirmatory hypotheses: planned E - N contrasts (both metrics each).
CONFIRMATORY = cell2table({ ...
    'H1', 'FrontalMidline', 'Theta', 'E - N'; ...
    'H2', 'Occipital',      'Alpha', 'E - N'}, ...
    'VariableNames', {'Hypothesis', 'ROI', 'Band', 'Contrast'});
CONFIRMATORY_METRICS = {'LogAbsolutePower', 'RelativePower'};

PRESPECIFICATION_STATUS = ['Literature-based ROIs and hypotheses, written AFTER the Phase 5B ' ...
    'electrode-level results were seen (not preregistered). Report them as literature-based, not preregistered.'];

%% 2. LOAD REFERENCE CHANNELS + COORDINATES
S = load(roiDefFile, 'ChannelDefT');
ChannelDefT = S.ChannelDefT;
refUC = upper(ChannelDefT.Label);
fails = 0;

fprintf('PHASE 6A -- ROI v2 DEFINITIONS\n');
fprintf('------------------------------------------------\n');

%% 3. VALIDATION
roiOfChannel = repmat({''}, height(ChannelDefT), 1);
for r = 1:numel(ROI2_NAMES)
    labs = ROI2_DEFINITION.(ROI2_NAMES{r});
    for k = 1:numel(labs)
        idx = find(strcmp(refUC, upper(labs{k})));
        if numel(idx) ~= 1
            fprintf('  *** FAIL: %s (ROI %s) not found exactly once among the 62 reference channels.\n', labs{k}, ROI2_NAMES{r});
            fails = fails + 1;
            continue;
        end
        if ~isempty(roiOfChannel{idx})
            fprintf('  *** FAIL: %s is in both %s and %s (ROIs must not overlap).\n', labs{k}, roiOfChannel{idx}, ROI2_NAMES{r});
            fails = fails + 1;
        end
        roiOfChannel{idx} = ROI2_NAMES{r};
    end
end
for k = 1:numel(EDGE_EXCLUDED)
    idx = find(strcmp(refUC, upper(EDGE_EXCLUDED{k})));
    if numel(idx) ~= 1
        fprintf('  *** FAIL: edge electrode %s not found exactly once.\n', EDGE_EXCLUDED{k});
        fails = fails + 1;
    elseif ~isempty(roiOfChannel{idx})
        fprintf('  *** FAIL: edge electrode %s is inside ROI %s.\n', EDGE_EXCLUDED{k}, roiOfChannel{idx});
        fails = fails + 1;
    end
end
% Confirmatory list refers to real ROIs / bands
for i = 1:height(CONFIRMATORY)
    if ~ismember(CONFIRMATORY.ROI{i}, ROI2_NAMES)
        fprintf('  *** FAIL: confirmatory %s refers to unknown ROI %s.\n', CONFIRMATORY.Hypothesis{i}, CONFIRMATORY.ROI{i});
        fails = fails + 1;
    end
    if ~ismember(CONFIRMATORY.Contrast{i}, {'E - N','E - R','N - R'})
        fprintf('  *** FAIL: confirmatory %s has unknown contrast %s.\n', CONFIRMATORY.Hypothesis{i}, CONFIRMATORY.Contrast{i});
        fails = fails + 1;
    end
end

% Coordinate sanity (EEGLAB convention: +X anterior, +Y left)
meanX = nan(numel(ROI2_NAMES), 1); meanY = nan(numel(ROI2_NAMES), 1);
for r = 1:numel(ROI2_NAMES)
    sel = strcmp(roiOfChannel, ROI2_NAMES{r});
    meanX(r) = mean(ChannelDefT.X(sel)); meanY(r) = mean(ChannelDefT.Y(sel));
    fprintf('  %-15s n=%d  mean X (anterior+) = %+.3f   mean Y (left+) = %+.3f\n', ROI2_NAMES{r}, sum(sel), meanX(r), meanY(r));
    if abs(meanY(r)) > 0.05
        fprintf('  *** FAIL: %s is not left/right symmetric (|mean Y| = %.3f > 0.05).\n', ROI2_NAMES{r}, abs(meanY(r)));
        fails = fails + 1;
    end
end
if ~issorted(-meanX)   % FrontalMidline > Central > Parietal > Occipital along X
    fprintf('  *** FAIL: ROIs are not ordered front-to-back by mean X as expected.\n');
    fails = fails + 1;
end

fprintf('\nValidation failures: %d\n', fails);
if fails > 0
    error('Phase6A:ValidationFailed', 'ROI v2 validation failed (%d issue(s)). Nothing saved.', fails);
end

%% 4. SAVE
isEdge = ismember(refUC, upper(EDGE_EXCLUDED));
Roi2DefT = table((1:height(ChannelDefT))', ChannelDefT.Label, roiOfChannel, isEdge, ...
    ChannelDefT.X, ChannelDefT.Y, ChannelDefT.Z, ...
    'VariableNames', {'Index','Label','ROI_v2','EdgeExcluded','X','Y','Z'});

PARAMS11 = struct();
PARAMS11.ScriptName = 'X24_roi_v2_definitions.m';
PARAMS11.GeneratedOn = datestr(now);
PARAMS11.SourcePhase5DFile = roiDefFile;
PARAMS11.Rationale = ROI2_RATIONALE;
PARAMS11.EdgeExcluded = {EDGE_EXCLUDED};
PARAMS11.ConfirmatoryMetrics = {CONFIRMATORY_METRICS};
PARAMS11.PrespecificationStatus = PRESPECIFICATION_STATUS;

save(fullfile(OUT_DIR, 'roi_v2_definitions.mat'), 'Roi2DefT', 'ROI2_DEFINITION', 'ROI2_NAMES', ...
    'CONFIRMATORY', 'CONFIRMATORY_METRICS', 'PARAMS11');
writetable(Roi2DefT, fullfile(OUT_DIR, 'roi_v2_definitions.csv'));

%% 5. REPORT
L = {};
L{end+1} = 'PHASE 6A -- LITERATURE-BASED ROI v2 DEFINITIONS';
L{end+1} = repmat('=', 1, 60);
L{end+1} = sprintf('Generated: %s', datestr(now));
L{end+1} = '';
for r = 1:numel(ROI2_NAMES)
    L{end+1} = sprintf('%-15s (%d): %s', ROI2_NAMES{r}, numel(ROI2_DEFINITION.(ROI2_NAMES{r})), strjoin(ROI2_DEFINITION.(ROI2_NAMES{r}), ' ')); %#ok<SAGROW>
    L{end+1} = sprintf('    rationale: %s', ROI2_RATIONALE.(ROI2_NAMES{r})); %#ok<SAGROW>
    L{end+1} = sprintf('    mean X=%+.3f  mean Y=%+.3f', meanX(r), meanY(r)); %#ok<SAGROW>
end
L{end+1} = '';
L{end+1} = sprintf('Edge electrodes excluded from all ROIs (artifact-prone): %s', strjoin(EDGE_EXCLUDED, ' '));
L{end+1} = sprintf('Electrodes in no ROI: %d of %d', sum(cellfun(@isempty, roiOfChannel)), height(ChannelDefT));
L{end+1} = '';
L{end+1} = 'CONFIRMATORY HYPOTHESES (planned contrasts, tested in 12, Holm across all rows x metrics):';
for i = 1:height(CONFIRMATORY)
    L{end+1} = sprintf('  %s: %s  %s  %s   metrics: %s', CONFIRMATORY.Hypothesis{i}, CONFIRMATORY.ROI{i}, ...
        CONFIRMATORY.Band{i}, CONFIRMATORY.Contrast{i}, strjoin(CONFIRMATORY_METRICS, ', ')); %#ok<SAGROW>
end
L{end+1} = '';
L{end+1} = ['PRESPECIFICATION STATUS: ' PRESPECIFICATION_STATUS];
L{end+1} = '';
L{end+1} = 'Validation failures: 0';

fid = fopen(fullfile(REPORT_DIR, 'phase6a_roi_v2_definitions_report.txt'), 'w');
fprintf(fid, '%s\n', L{:});
fclose(fid);
fprintf('%s\n', L{:});
fprintf('\n=== X24_roi_v2_definitions.m complete. Run X25_roi_v2_statistics.m next. ===\n');
