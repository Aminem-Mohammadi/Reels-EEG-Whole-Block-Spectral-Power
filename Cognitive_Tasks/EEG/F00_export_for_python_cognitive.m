%% ========================================================================
% F00_export_for_python_cognitive.m   (cognitive_task pipeline -- figure hand-off)
%
% TASK-SIDE PORT OF F00_export_for_python.m (Intervention pipeline). Writes
% the few inputs the Python figure script cannot read directly, because
% they live in MATLAB tables / struct arrays or in -v7.3 .mat files.
% Everything else (X05/X06, X08, X09, X09B, X12, X24/X25 CSV/XLSX outputs)
% is read by Python straight from their own, already-Task-looped folders.
%
% Run once, after X03, X04 and X07. No statistics are computed here.
%
% WHAT IS DIFFERENT FROM THE INTERVENTION VERSION, AND WHY (full disclosure)
%   - 'InterventionNumber' does not exist in this pipeline. Every recording
%     here carries BOTH 'TaskType' (full/nback/sart) and 'BlockNumber'
%     instead; both columns are exported, in place of InterventionNumber,
%     wherever the Intervention version exported InterventionNumber.
%   - This script does NOT loop over Task and does NOT produce separate
%     per-Task export folders/files. This is not a simplification choice
%     made for this port -- it reflects how X03/X04 themselves already
%     store this pipeline's data: BandPowerT, GlobalBandPowerT (from X04)
%     and psd_data (from X03) are each a SINGLE table/struct array
%     containing all three Tasks' rows together, with TaskType as an
%     ordinary column -- never three separate Task-specific tables.
%     Exporting three separate copies would not reflect any split that
%     exists upstream, so a single export (with TaskType included, exactly
%     as X04/X03 already have it) is the one-to-one port of the
%     Intervention version's own "export exactly what X03/X04 already
%     are" principle. The Task-adapted F01-equivalent Python script is
%     responsible for filtering by TaskType itself when building each
%     Task's own figures, the same way X05/X06/X08/X09/X25 already filter
%     by TaskType internally before fitting any model.
%   - Source file names/paths differ to match this pipeline's own X03/X04
%     output (cognitive_psd.mat instead of intervention_psd.mat;
%     ANALYSIS_ROOT derived the same way X05/X08/X09/X24/X25_cognitive
%     already derive it, from DATA_ROOT_FOR_SIBLING).
%   - Everything else -- which files are written, their column sets (minus
%     the InterventionNumber-to-TaskType+BlockNumber substitution), the
%     PSD export shape and -v7 scipy-readable format, and the "no
%     statistics computed here" scope -- is reproduced verbatim.
%
% OUTPUT  <ANALYSIS_ROOT>/30_python_inputs/
%   channel_locations.csv   Label, X, Y, Z, Theta, Radius   (X07 ChannelDefT;
%                           EEGLAB convention: +X = nose, +Y = left ear, +Z = up)
%   bandpower_electrode.csv Subject, TaskType, BlockNumber, Condition, Channel,
%                           Band, AbsolutePower, LogAbsolutePower, RelativePower
%                           (X04; all three Tasks together, TaskType column
%                           distinguishes them -- see disclosure above)
%   bandpower_global.csv    same without Channel                     (X04)
%   psd_recordings.mat      (-v7, readable by scipy) PSD [nRec x nChan x nFreq]
%                           (mean over valid windows, uV^2/Hz), FreqHz,
%                           ChannelLabels, RecSubject, RecTaskType,
%                           RecBlockNumber, RecCondition   (X03; all three
%                           Tasks together, RecTaskType column distinguishes
%                           them)
% ========================================================================

clear; clc;

config;
DATA_ROOT_FOR_SIBLING = DATA_DIR;
[DATA_ROOT_PARENT, ~] = fileparts(DATA_ROOT_FOR_SIBLING);
ANALYSIS_ROOT = OUT_DIR;

OUT_DIR = fullfile(ANALYSIS_ROOT, '30_python_inputs');
if ~exist(OUT_DIR, 'dir'); mkdir(OUT_DIR); end

defFile = fullfile(ANALYSIS_ROOT, '07_roi_electrode_definitions', 'phase5_channel_roi_definitions.mat');
bpFile  = fullfile(ANALYSIS_ROOT, '04_subject_level', 'phase3_subject_level_bandpower.mat');
psdFile = fullfile(ANALYSIS_ROOT, '03_psd', 'cognitive_psd.mat');
for f = {defFile, bpFile, psdFile}
    if ~exist(f{1}, 'file'); error('F00:MissingInput', 'Not found: %s', f{1}); end
end

fprintf('========================================================\n');
fprintf('F00 (cognitive_task) -- EXPORT INPUTS FOR PYTHON PUBLICATION FIGURES\n');
fprintf('========================================================\n');
fprintf('ANALYSIS_ROOT : %s\n', ANALYSIS_ROOT);
fprintf('OUT_DIR        : %s\n\n', OUT_DIR);

%% 1. Channel locations (shared, Task-invariant -- same montage as Intervention)
S = load(defFile, 'ChannelDefT'); C = S.ChannelDefT;
writetable(C(:, {'Label', 'X', 'Y', 'Z', 'Theta', 'Radius'}), fullfile(OUT_DIR, 'channel_locations.csv'));
fprintf('channel_locations.csv : %d channels\n', height(C));

%% 2. Subject-level band power (X04) -- ALL three Tasks together, TaskType column
S = load(bpFile, 'BandPowerT', 'GlobalBandPowerT');
colsE = {'Subject', 'TaskType', 'BlockNumber', 'Condition', 'Channel', 'Band', 'AbsolutePower', 'LogAbsolutePower', 'RelativePower'};
colsG = {'Subject', 'TaskType', 'BlockNumber', 'Condition', 'Band', 'AbsolutePower', 'LogAbsolutePower', 'RelativePower'};
missE = colsE(~ismember(colsE, S.BandPowerT.Properties.VariableNames));
missG = colsG(~ismember(colsG, S.GlobalBandPowerT.Properties.VariableNames));
if ~isempty(missE); error('F00:MissingColumns', 'BandPowerT is missing: %s', strjoin(missE, ', ')); end
if ~isempty(missG); error('F00:MissingColumns', 'GlobalBandPowerT is missing: %s', strjoin(missG, ', ')); end
writetable(S.BandPowerT(:, colsE), fullfile(OUT_DIR, 'bandpower_electrode.csv'));
writetable(S.GlobalBandPowerT(:, colsG), fullfile(OUT_DIR, 'bandpower_global.csv'));
fprintf('bandpower_electrode.csv: %d rows (all 3 Tasks) | bandpower_global.csv: %d rows (all 3 Tasks)\n', ...
    height(S.BandPowerT), height(S.GlobalBandPowerT));
taskCountsE = varfun(@numel, S.BandPowerT, 'InputVariables', 'Subject', 'GroupingVariables', 'TaskType');
for i = 1:height(taskCountsE)
    fprintf('  TaskType=%-6s : %d electrode-level rows\n', taskCountsE.TaskType{i}, taskCountsE.GroupCount(i));
end

%% 3. Recording-level mean PSD (X03) -- ALL three Tasks together, RecTaskType column
S = load(psdFile, 'psd_data', 'FrequencyVectorHz', 'ChannelLabelsReference');
P = S.psd_data;
nRec = numel(P);
FreqHz = double(S.FrequencyVectorHz(:))';
ChannelLabels = cellstr(S.ChannelLabelsReference(:));
nChan = numel(ChannelLabels); nFreq = numel(FreqHz);
PSD = nan(nRec, nChan, nFreq);
RecSubject = nan(nRec, 1); RecBlockNumber = nan(nRec, 1);
RecTaskType = cell(nRec, 1); RecCondition = cell(nRec, 1);
for r = 1:nRec
    if ~isequal(cellstr(P(r).ChannelLabels(:)), ChannelLabels)
        error('F00:ChannelOrder', 'Recording %d has a different channel order.', r);
    end
    M = P(r).MeanPSD;                         % nChan x nFreq
    if ~isequal(size(M), [nChan nFreq]); error('F00:PSDSize', 'Recording %d: MeanPSD size mismatch.', r); end
    PSD(r, :, :) = reshape(M, [1 nChan nFreq]);
    RecSubject(r) = str2double(regexprep(char(string(P(r).Subject)), '\D', ''));
    RecBlockNumber(r) = double(P(r).BlockNumber);
    RecTaskType{r} = char(P(r).TaskType);
    RecCondition{r} = char(P(r).Condition);
end
save(fullfile(OUT_DIR, 'psd_recordings.mat'), 'PSD', 'FreqHz', 'ChannelLabels', ...
    'RecSubject', 'RecTaskType', 'RecBlockNumber', 'RecCondition', '-v7');
fprintf('psd_recordings.mat     : %d recordings x %d channels x %d frequencies (all 3 Tasks together)\n', nRec, nChan, nFreq);
for tk = {'full','nback','sart'}
    fprintf('  TaskType=%-6s : %d recordings\n', tk{1}, sum(strcmpi(RecTaskType, tk{1})));
end

fprintf('\n=== F00 (cognitive_task) complete. Now run:  python3 F01_publication_figures_cognitive.py ===\n');
