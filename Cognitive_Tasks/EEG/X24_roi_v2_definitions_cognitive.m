%% ========================================================================
% X24_roi_v2_definitions_cognitive.m
%
% PORTED from the intervention pipeline's X24_roi_v2_definitions.m so the
% cognitive-task results have the same confirmatory-analysis layer as the
% intervention results (same literature-based ROI set, same two
% confirmatory hypotheses' ROI/Band, so the two datasets are directly
% comparable in the shared paper). What changed in the port:
%   - Input/output root: <ANALYSIS_ROOT> = cognitive_task_analysis (this
%     pipeline's own X07 output), not the intervention folders.
%   - ROI_DEFINITION, EDGE_EXCLUDED and the coordinate-sanity checks are
%     UNCHANGED from the intervention X24 -- the 62-electrode montage and
%     every electrode label used here (FZ, FCz, F1, F2, FC1, FC2, C3, C1,
%     Cz, C2, C4, P3, P1, PZ, P2, P4, PO3, POz, PO4, O1, Oz, O2, plus the
%     14 edge electrodes) are verified present, by this same script's own
%     Section 3 validation, in THIS pipeline's X07 output
%     (X07_roi_electrode_definitions_cognitive.m's CANONICAL_ORDER /
%     ChannelDefT) before anything is saved -- the script stops rather
%     than assuming the two montages match.
%   - CONFIRMATORY HYPOTHESES: this pipeline's own E/N/R contrast
%     convention is N-E / R-E / R-N (X05/X06/X11-X20's convention -- see
%     those scripts' headers), the OPPOSITE sign/order of the
%     intervention pipeline's E-N/E-R/N-R convention. The two confirmatory
%     hypotheses below were chosen, disclosed exactly like the
%     intervention X24, AFTER seeing this pipeline's own Phase 4/5 results
%     (the FullSample-vs-SensitivityExclNoisy comparison already
%     delivered): Occipital/Parietal Alpha-and-Beta-family
%     LogAbsolutePower effects via the R-E/R-N contrasts were the single
%     most consistent pattern across every RUN_MODE-aware script reviewed,
%     and Frontal-midline theta is the literature's own primary marker of
%     the sustained-attention/cognitive-control construct that Nback/SART
%     are designed to tax -- making it at least as well-motivated here as
%     for the intervention dataset. See H1/H2 below. As with the
%     intervention X24: DO
%     NOT edit the ROI list or the hypothesis list after seeing this
%     Task's own X25_roi_v2_statistics_cognitive.m output.
%   - X24_roi_v2_definitions.m's own header origin note ("renamed from
%     X11, number reused by the median-robust script") does not apply
%     here -- this pipeline's X11 is its own, unrelated, already-delivered
%     median-robust-family script (X11_median_robust_family_cognitive.m);
%     X24/X25 are new numbers, chosen because 16 and 21 are already used
%     by this pipeline's own X16/X21, and 22/23 by its own topography
%     scripts, so 24/25 are the next free slot, kept IDENTICAL to the
%     intervention pipeline's own X24/X25 numbering for cross-dataset
%     readability.
%   - Unlike X07_roi_electrode_definitions_cognitive.m (also shared across
%     all three Tasks), this script performs NO Task loop and writes NO
%     per-Task output -- ROI v2 membership and the two hypotheses do not
%     depend on Task, exactly as the intervention X24 does not depend on
%     InterventionNumber. X25_roi_v2_statistics_cognitive.m (which DOES
%     loop over Full/Nback/SART) reads this single, shared output file.
% ========================================================================
%
% PHASE 6A -- LITERATURE-BASED ROI DEFINITIONS (v2) + CONFIRMATORY TESTS,
% COGNITIVE-TASK DATA.
%
% PURPOSE
%   Define a SMALL set of compact, symmetric, literature-based ROIs (4-6
%   electrodes each) and lock the short list of confirmatory hypotheses
%   that X25_roi_v2_statistics_cognitive.m will test as planned contrasts,
%   independently for each of Full/Nback/SART.
%
%   This replaces neither X07 nor X08: the original 5-ROI partition stays
%   as it is. This is an additional, hypothesis-driven ROI set.
%
% WHY THESE ROIs (identical rationale and literature to the intervention
% pipeline's X24 -- the montage is the same 62-channel system)
%   - Every ROI comes from the literature, NOT from this dataset's results.
%   - Each ROI is small (4-6 electrodes), so a focal effect is not diluted
%     by averaging with distant electrodes (the problem with the
%     23-electrode v1 Frontal ROI).
%   - Edge electrodes (FP1/FP2, AF7/AF8, F7/F8, FT7/FT8, T7/T8, TP7/TP8,
%     P9/P10) are not used in any ROI: they carry the most eye-movement and
%     muscle artifact, which is an independent, a priori reason.
%   - ROIs are left/right symmetric and do not overlap.
%
%     FrontalMidline (6): FZ FCz F1 F2 FC1 FC2
%         Frontal-midline theta (FMtheta): sustained attention / cognitive
%         control. Ishihara & Yoshii 1972; Mitchell et al. 2008 (Prog
%         Neurobiol); Cavanagh & Frank 2014 (Trends Cogn Sci).
%     Central (5): C3 C1 Cz C2 C4
%         Sensorimotor (mu/beta).  Pfurtscheller & Lopes da Silva 1999.
%     Parietal (5): P3 P1 PZ P2 P4
%         Parietal alpha -- attention allocation. Klimesch 1999, 2012.
%     Occipital (6): PO3 POz PO4 O1 Oz O2
%         Posterior alpha -- visual processing / visual engagement.
%         Klimesch 1999; Berger 1929.
%
% CONFIRMATORY HYPOTHESES (tested in X25 as PLANNED contrasts, not gated
% by an omnibus test; Holm-corrected across this whole list x both metrics
% x each Task separately -- Tasks are NEVER pooled, see X25):
%     H1  Occipital       Alpha  R - E
%     H2  FrontalMidline  Theta  N - E
%
% HONESTY NOTE -- READ BEFORE REPORTING (identical disclosure to the
% intervention pipeline's X24)
%   These ROIs and hypotheses are taken from the literature, but they were
%   written down AFTER this pipeline's own Phase 4/5 electrode/ROI results
%   were seen (including the FullSample-vs-SensitivityExclNoisy
%   comparison). They are therefore "literature-based, not preregistered".
%   Say so in the paper. Do NOT edit the ROI lists or the hypothesis list
%   after seeing X25's output for any Task -- doing that turns a
%   confirmatory test into a data-driven one and the p-values stop meaning
%   what they claim. PRESPECIFICATION_STATUS below is saved into the
%   output and printed in X25's report so this note travels with the
%   results.
%
% REQUIREMENTS
%   - Base MATLAB only.
%   - X07 output: 07_roi_electrode_definitions/phase5_channel_roi_definitions.mat
%     (source of the locked 62-channel order and real X/Y/Z coordinates,
%     SHARED across all three Tasks -- see X07_roi_electrode_definitions_
%     cognitive.m).
%
% OUTPUT
%   <ANALYSIS_ROOT>/24_roi_v2_definitions/
%       roi_v2_definitions.mat   (Roi2DefT, ROI2_DEFINITION, ROI2_NAMES,
%                                 CONFIRMATORY, CONFIRMATORY_METRICS, PARAMS24)
%       roi_v2_definitions.csv
%   <ANALYSIS_ROOT>/09_reports/phase6a_roi_v2_definitions_report.txt
% ========================================================================

clear; clc;

%% 0. SETTINGS
config;
DATA_ROOT_FOR_SIBLING = DATA_DIR;
[DATA_ROOT_PARENT, ~] = fileparts(DATA_ROOT_FOR_SIBLING);
ANALYSIS_ROOT = OUT_DIR;

DEF_DIR    = fullfile(ANALYSIS_ROOT, '07_roi_electrode_definitions');
OUT_DIR    = fullfile(ANALYSIS_ROOT, '24_roi_v2_definitions');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');
roiDefFile = fullfile(DEF_DIR, 'phase5_channel_roi_definitions.mat');

if ~exist(roiDefFile, 'file')
    error('Phase6A:MissingInput', 'X07 output file not found: %s\nRun X07_roi_electrode_definitions_cognitive.m first.', roiDefFile);
end
if ~exist(OUT_DIR, 'dir');    mkdir(OUT_DIR);    end
if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

%% 1. ROI DEFINITIONS (fixed constants -- see header; identical to the
%     intervention pipeline's X24, same 62-channel montage)
ROI2_DEFINITION = struct();
ROI2_DEFINITION.FrontalMidline = {'FZ','FCz','F1','F2','FC1','FC2'};
ROI2_DEFINITION.Central        = {'C3','C1','Cz','C2','C4'};
ROI2_DEFINITION.Parietal       = {'P3','P1','PZ','P2','P4'};
ROI2_DEFINITION.Occipital      = {'PO3','POz','PO4','O1','Oz','O2'};
ROI2_NAMES = {'FrontalMidline','Central','Parietal','Occipital'};

ROI2_RATIONALE = struct( ...
    'FrontalMidline', 'Frontal-midline theta (attention / cognitive control): Ishihara & Yoshii 1972; Mitchell et al. 2008; Cavanagh & Frank 2014', ...
    'Central',        'Sensorimotor mu/beta: Pfurtscheller & Lopes da Silva 1999', ...
    'Parietal',       'Parietal alpha (attention allocation): Klimesch 1999, 2012', ...
    'Occipital',      'Posterior alpha (visual processing / engagement): Berger 1929; Klimesch 1999');

EDGE_EXCLUDED = {'FP1','FP2','AF7','AF8','F7','F8','FT7','FT8','T7','T8','TP7','TP8','P9','P10'};

% Confirmatory hypotheses: planned contrasts, this pipeline's own N-E/R-E/R-N
% convention (see X05/X06/X11-X20 headers), both metrics each, tested
% independently per Task in X25 -- never pooled across Full/Nback/SART.
CONFIRMATORY = cell2table({ ...
    'H1', 'Occipital',      'Alpha', 'R - E'; ...
    'H2', 'FrontalMidline', 'Theta', 'N - E'}, ...
    'VariableNames', {'Hypothesis', 'ROI', 'Band', 'Contrast'});
CONFIRMATORY_METRICS = {'LogAbsolutePower', 'RelativePower'};

PRESPECIFICATION_STATUS = ['Literature-based ROIs and hypotheses, written AFTER this pipeline''s own ' ...
    'Phase 4/5 electrode/ROI results (including the FullSample-vs-SensitivityExclNoisy comparison) ' ...
    'were seen (not preregistered). Report them as literature-based, not preregistered. Chosen to ' ...
    'directly parallel the intervention pipeline''s X24 (Occipital/Alpha <-> its H2; FrontalMidline/' ...
    'Theta <-> its H1) for cross-dataset comparability in the shared paper.'];

%% 2. LOAD REFERENCE CHANNELS + COORDINATES (shared X07 output)
S = load(roiDefFile, 'ChannelDefT');
ChannelDefT = S.ChannelDefT;
refUC = upper(ChannelDefT.Label);
fails = 0;

fprintf('PHASE 6A -- ROI v2 DEFINITIONS (COGNITIVE-TASK DATA)\n');
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
% Confirmatory list refers to real ROIs / bands / contrasts
for i = 1:height(CONFIRMATORY)
    if ~ismember(CONFIRMATORY.ROI{i}, ROI2_NAMES)
        fprintf('  *** FAIL: confirmatory %s refers to unknown ROI %s.\n', CONFIRMATORY.Hypothesis{i}, CONFIRMATORY.ROI{i});
        fails = fails + 1;
    end
    if ~ismember(CONFIRMATORY.Contrast{i}, {'N - E','R - E','R - N'})
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

PARAMS24 = struct();
PARAMS24.ScriptName = 'X24_roi_v2_definitions_cognitive.m';
PARAMS24.GeneratedOn = datestr(now); %#ok<TNOW1,DATST>
PARAMS24.SourceX07File = roiDefFile;
PARAMS24.Rationale = ROI2_RATIONALE;
PARAMS24.EdgeExcluded = {EDGE_EXCLUDED};
PARAMS24.ConfirmatoryMetrics = {CONFIRMATORY_METRICS};
PARAMS24.ContrastConvention = 'N-E / R-E / R-N (this pipeline''s own convention; X05/X06/X11-X20 -- OPPOSITE sign/order of the intervention pipeline''s E-N/E-R/N-R).';
PARAMS24.PrespecificationStatus = PRESPECIFICATION_STATUS;

save(fullfile(OUT_DIR, 'roi_v2_definitions.mat'), 'Roi2DefT', 'ROI2_DEFINITION', 'ROI2_NAMES', ...
    'CONFIRMATORY', 'CONFIRMATORY_METRICS', 'PARAMS24');
writetable(Roi2DefT, fullfile(OUT_DIR, 'roi_v2_definitions.csv'));

%% 5. REPORT
L = {};
L{end+1} = 'PHASE 6A -- LITERATURE-BASED ROI v2 DEFINITIONS (COGNITIVE-TASK DATA)';
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
L{end+1} = 'CONFIRMATORY HYPOTHESES (planned contrasts, tested independently per Task in X25, Holm across all rows x metrics x Task):';
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
fprintf('\n=== X24_roi_v2_definitions_cognitive.m complete. Run X25_roi_v2_statistics_cognitive.m next. ===\n');
