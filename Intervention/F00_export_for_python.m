%% ========================================================================
% F00_export_for_python.m
%
% Exports the inputs of the Python figure script (F01_publication_figures.py).
% Inputs : 03_psd, 04_subject_level, 07_roi_electrode_definitions.
% Outputs: 30_python_inputs/ -- channel_locations.csv (Label, X, Y, Z, Theta, Radius;
%          EEGLAB convention +X = nose, +Y = left ear), bandpower_electrode.csv,
%          bandpower_global.csv, psd_recordings.mat (recording-mean PSD, -v7).
% No statistics are computed.
% ========================================================================

clear; clc;
config;
ANALYSIS_ROOT = OUT_DIR;
OUT_DIR = fullfile(ANALYSIS_ROOT, '30_python_inputs');
if ~exist(OUT_DIR, 'dir'); mkdir(OUT_DIR); end

defFile = fullfile(ANALYSIS_ROOT, '07_roi_electrode_definitions', 'phase5_channel_roi_definitions.mat');
bpFile  = fullfile(ANALYSIS_ROOT, '04_subject_level', 'phase3_subject_level_bandpower.mat');
psdFile = fullfile(ANALYSIS_ROOT, '03_psd', 'intervention_psd.mat');
for f = {defFile, bpFile, psdFile}
    if ~exist(f{1}, 'file'); error('F00:MissingInput', 'Not found: %s', f{1}); end
end

%% 1. Channel locations
S = load(defFile, 'ChannelDefT'); C = S.ChannelDefT;
writetable(C(:, {'Label', 'X', 'Y', 'Z', 'Theta', 'Radius'}), fullfile(OUT_DIR, 'channel_locations.csv'));
fprintf('channel_locations.csv : %d channels\n', height(C));

%% 2. Subject-level band power (X04)
S = load(bpFile, 'BandPowerT', 'GlobalBandPowerT');
colsE = {'Subject', 'InterventionNumber', 'Condition', 'Channel', 'Band', 'AbsolutePower', 'LogAbsolutePower', 'RelativePower'};
colsG = {'Subject', 'InterventionNumber', 'Condition', 'Band', 'AbsolutePower', 'LogAbsolutePower', 'RelativePower'};
writetable(S.BandPowerT(:, colsE), fullfile(OUT_DIR, 'bandpower_electrode.csv'));
writetable(S.GlobalBandPowerT(:, colsG), fullfile(OUT_DIR, 'bandpower_global.csv'));
fprintf('bandpower_electrode.csv: %d rows | bandpower_global.csv: %d rows\n', height(S.BandPowerT), height(S.GlobalBandPowerT));

%% 3. Recording-level mean PSD (X03)
S = load(psdFile, 'psd_data', 'FrequencyVectorHz', 'ChannelLabelsReference');
P = S.psd_data;
nRec = numel(P);
FreqHz = double(S.FrequencyVectorHz(:))';
ChannelLabels = cellstr(S.ChannelLabelsReference(:));
nChan = numel(ChannelLabels); nFreq = numel(FreqHz);
PSD = nan(nRec, nChan, nFreq);
RecSubject = nan(nRec, 1); RecIntervention = nan(nRec, 1); RecCondition = cell(nRec, 1);
for r = 1:nRec
    if ~isequal(cellstr(P(r).ChannelLabels(:)), ChannelLabels)
        error('F00:ChannelOrder', 'Recording %d has a different channel order.', r);
    end
    M = P(r).MeanPSD;                         % nChan x nFreq
    if ~isequal(size(M), [nChan nFreq]); error('F00:PSDSize', 'Recording %d: MeanPSD size mismatch.', r); end
    PSD(r, :, :) = reshape(M, [1 nChan nFreq]);
    RecSubject(r) = str2double(regexprep(char(string(P(r).Subject)), '\D', ''));
    RecIntervention(r) = double(P(r).InterventionNumber);
    RecCondition{r} = char(P(r).Condition);
end
save(fullfile(OUT_DIR, 'psd_recordings.mat'), 'PSD', 'FreqHz', 'ChannelLabels', ...
    'RecSubject', 'RecIntervention', 'RecCondition', '-v7');
fprintf('psd_recordings.mat     : %d recordings x %d channels x %d frequencies\n', nRec, nChan, nFreq);
fprintf('\n=== F00 complete. Now run:  python3 F01_publication_figures.py ===\n');
