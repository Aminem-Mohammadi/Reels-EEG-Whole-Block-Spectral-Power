%% ========================================================================
% X15_aperiodic_1f_family_intervention.m
%
% Aperiodic (1/f) component of the power spectrum: exponent and offset estimated with
% FOOOF (through X15_fit_fooof.py, called with python3) and with a robust log-log linear
% fit; global, regional and electrode level.
% Inputs : 03_psd, 07_roi_electrode_definitions.
% Outputs: 15_aperiodic_1f_family/<RUN_MODE>/, 09_reports/.
% Run twice, with RUN_MODE = 'full_sample' (N = 30) and RUN_MODE = 'sensitivity_excl_noisy'
% (cleaned sample, N = 28; excluded subjects listed in EXCLUSION_CONFIG_intervention.m).
% Outputs are written to <OUT_DIR>/<folder>/<RUN_MODE>/ and reports to <OUT_DIR>/09_reports/.
% ========================================================================

clear; clc;

%% ------------------------------------------------------------------
%  SECTION 0: PATHS, SETTINGS, EXISTENCE CHECKS
% ------------------------------------------------------------------

config;
DATA_ROOT_FOR_SIBLING = DATA_DIR;
[DATA_ROOT_PARENT, ~] = fileparts(DATA_ROOT_FOR_SIBLING);
ANALYSIS_ROOT = OUT_DIR;

THIS_SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(THIS_SCRIPT_DIR);
run(fullfile(THIS_SCRIPT_DIR, 'EXCLUSION_CONFIG_intervention.m'));
RUN_MODE = 'full_sample';   % 'full_sample' (N = 30) or 'sensitivity_excl_noisy' (N = 28)
if ~ismember(RUN_MODE, RUN_MODES_intervention)
    error('RunMode:Bad', 'RUN_MODE must be one of: %s', strjoin(RUN_MODES_intervention, ', '));
end
RUN_MODE_LABEL = RUN_MODES_DISPLAY_intervention.(RUN_MODE);
RUN_MODE_EXCL_SUBJ = [];   % filled in after the data are loaded
fprintf('RUN_MODE = %s\n', RUN_MODE_LABEL);
if strcmp(RUN_MODE, 'sensitivity_excl_noisy')
    fprintf('Excluding subject(s) {%s} -- see EXCLUSION_CONFIG_intervention.m for the reasons.\n', ...
        strjoin(cellstr(num2str(EXCLUDE_NOISY_SUBJECTS(:))), ','));
end

PSD_DIR    = fullfile(ANALYSIS_ROOT, '03_psd');
DEF_DIR    = fullfile(ANALYSIS_ROOT, '07_roi_electrode_definitions');
REPORT_DIR = fullfile(ANALYSIS_ROOT, '09_reports');

psdMatFile = fullfile(PSD_DIR, 'intervention_psd.mat');
roiDefFile = fullfile(DEF_DIR, 'phase5_channel_roi_definitions.mat');

THIS_SCRIPT_DIR = fileparts(mfilename('fullpath'));
PYTHON_SCRIPT_PATH = fullfile(THIS_SCRIPT_DIR, 'X15_fit_fooof.py');

% ----------------------------------------------------------------------
% Resolve a working `python3` executable. A GUI-launched MATLAB (opened
% from Finder/Dock/Spotlight, not a Terminal) does NOT inherit the
% shell's PATH, so `system('python3 ...')` can fail with "command not
% found" even though `python3` works fine in Terminal. This is
% especially common on Apple Silicon Macs with Homebrew (python3 lives
% at /opt/homebrew/bin/python3, which a GUI-launched MATLAB's PATH often
% does not include). Same dual/multi-candidate resolve-and-log pattern
% already used elsewhere in this pipeline (e.g. CHANLOCS_TEMPLATE_FILE
% in X07/X10): try the bare `python3` on MATLAB's own PATH first, then
% fall back to the standard Homebrew locations, logging whichever one
% actually works. Never silently guess -- hard error with an actionable
% message if none work.
% ----------------------------------------------------------------------
PYTHON_CANDIDATES = {
    'python3', ...                        % candidate 1: rely on MATLAB's inherited PATH
    '/opt/homebrew/bin/python3', ...      % candidate 2: Homebrew on Apple Silicon
    '/usr/local/bin/python3', ...         % candidate 3: Homebrew on Intel Mac / manual install
    '/usr/bin/python3' ...                % candidate 4: macOS system python3, if present
};
PYTHON_EXECUTABLE = '';
PYTHON_RESOLUTION_LOG = {};
for pc = 1:numel(PYTHON_CANDIDATES)
    cand = PYTHON_CANDIDATES{pc};
    [chkStatus, chkOut] = system(sprintf('"%s" --version', cand));
    PYTHON_RESOLUTION_LOG{end+1} = sprintf('  Candidate %d (%s): exit=%d, output=%s', ...
        pc, cand, chkStatus, strtrim(chkOut)); %#ok<AGROW>
    if chkStatus == 0
        PYTHON_EXECUTABLE = cand;
        break;
    end
end
if isempty(PYTHON_EXECUTABLE)
    error('Phase_Aperiodic:PythonNotFound', ...
        ['Could not find a working python3 executable via any candidate path. Tried:\n%s\n' ...
         'Fix: install Python 3 (e.g. `brew install python3` on macOS), then either (a) launch ' ...
         'MATLAB from a Terminal that has python3 on its PATH, or (b) add the exact candidate ' ...
         'path that works in Terminal (`which python3`) to PYTHON_CANDIDATES at the top of this ' ...
         'script.'], strjoin(PYTHON_RESOLUTION_LOG, '\n'));
end

if ~exist(REPORT_DIR, 'dir'); mkdir(REPORT_DIR); end

fprintf('========================================================\n');
fprintf('PHASE "APERIODIC 1/F FAMILY" -- INTERVENTION DATA (X15, exploratory extension)\n');
fprintf('Single intervention analysis (conditions E / N / R).\n');
fprintf('========================================================\n');
fprintf('ANALYSIS_ROOT      : %s\n', ANALYSIS_ROOT);
fprintf('psdMatFile          : %s\n', psdMatFile);
fprintf('roiDefFile           : %s\n', roiDefFile);
fprintf('PYTHON_SCRIPT_PATH   : %s\n', PYTHON_SCRIPT_PATH);
fprintf('PYTHON_EXECUTABLE    : %s  (resolved from %d candidate(s); see PYTHON_RESOLUTION_LOG)\n', PYTHON_EXECUTABLE, numel(PYTHON_CANDIDATES));
fprintf('%s\n', strjoin(PYTHON_RESOLUTION_LOG, '\n'));
fprintf('--------------------------------------------------------\n\n');

if ~exist(psdMatFile, 'file')
    error('Phase_Aperiodic:MissingInput', ['X03 PSD file not found: %s\n' ...
        'Run X03_compute_intervention_psd.m first.'], psdMatFile);
end
if ~exist(roiDefFile, 'file')
    error('Phase_Aperiodic:MissingInput', ['X07 ROI definition file not found: %s\n' ...
        'Run X07_roi_electrode_definitions.m first. This script will NOT ' ...
        'reconstruct or guess ROI membership on its own.'], roiDefFile);
end
if ~exist(PYTHON_SCRIPT_PATH, 'file')
    error('Phase_Aperiodic:MissingCompanionScript', ...
        'X15_fit_fooof.py not found next to this .m file (expected: %s).', PYTHON_SCRIPT_PATH);
end
if isempty(ver('stats'))
    error('Phase_Aperiodic:MissingToolbox', ...
        'Statistics and Machine Learning Toolbox not detected. This script requires fitlme, anova(LinearMixedModel), coefTest, and qqplot.');
end

%% ------------------------------------------------------------------

fprintf('Loading X03 PSD data (this may take a moment -- intervention_psd.mat is large) ...\n');
S3 = intervention_load(psdMatFile, 'psd_data', 'FrequencyVectorHz', 'ChannelLabelsReference', 'PARAMS');
psd_data = S3.psd_data;
FrequencyVectorHz = S3.FrequencyVectorHz(:);
ChannelLabelsReference = S3.ChannelLabelsReference(:);
PARAMS3 = S3.PARAMS;
clear S3;

nChan = numel(ChannelLabelsReference);
if nChan ~= 62
    error('Phase_Aperiodic:UnexpectedChannelCount', ...
        'ChannelLabelsReference has %d entries; this project is defined for a 62-channel montage.', nChan);
end

S5D = load(roiDefFile, 'ChannelDefT', 'ROI_DEFINITION', 'ROI_NAMES', 'UNASSIGNED_FOR_ROI');
ChannelDefT = S5D.ChannelDefT;
ROI_DEFINITION = S5D.ROI_DEFINITION;
ROI_NAMES = S5D.ROI_NAMES;
UNASSIGNED_FOR_ROI = S5D.UNASSIGNED_FOR_ROI;
clear S5D;
nROI = numel(ROI_NAMES);

if ~isequal(sort(upper(ChannelDefT.Label)), sort(upper(ChannelLabelsReference)))
    error('Phase_Aperiodic:ChannelMismatch', ...
        'X07''s ChannelDefT.Label does not match X03''s ChannelLabelsReference -- investigate before proceeding.');
end

% --- 1-45 Hz retained portion only (the 55-100 Hz portion is separated
%     by the 45-55 Hz gap and is not used). ---
freqMask45 = FrequencyVectorHz >= 1 & FrequencyVectorHz <= 45;
freqVec45 = FrequencyVectorHz(freqMask45);
nFreq45 = numel(freqVec45);
if nFreq45 < 10
    error('Phase_Aperiodic:TooFewFrequencyBins', ...
        'Only %d frequency bin(s) fall in 1-45 Hz -- too few for a meaningful aperiodic fit. Check FrequencyVectorHz.', nFreq45);
end
fprintf('1-45 Hz retained portion: %d frequency bins (%.4f-%.4f Hz).\n', nFreq45, freqVec45(1), freqVec45(end));

% --- Auto-detect the MedianPSDBased variant (same logic as X14). ---
hasMedianVariant = isfield(psd_data, 'PSD_WindowLevel');
if hasMedianVariant
    fprintf('psd_data carries the per-window full-PSD field (PSD_WindowLevel) -- MedianPSDBased variant WILL be computed alongside MeanPSDBased.\n');
else
    fprintf('psd_data does NOT carry a per-window full-PSD field (PARAMS.SaveFullWindowLevelPSD was false in X03, the default) -- MedianPSDBased variant is SKIPPED cleanly (not an error). Only MeanPSDBased is computed.\n');
end

% --- robustfit() availability check (Method 2). ---
USE_ROBUSTFIT = true;
robustFitFallbackNote = '';
try
    robustfit([1;2;3;4], [1.1;2.0;2.9;4.2]); %#ok<NASGU>
catch
    USE_ROBUSTFIT = false;
    robustFitFallbackNote = ['robustfit() was unavailable in this MATLAB session (Statistics and ' ...
        'Machine Learning Toolbox not fully licensed/installed for this function) -- Method 2 ' ...
        '(the sensitivity/fallback log-log linear fit) fell back to ordinary least-squares ' ...
        'polyfit() for EVERY row, in every Task. This is a REDUCED-ROBUSTNESS fallback ' ...
        '(no down-weighting of outlier frequency bins), disclosed here and in every Task''s report.'];
    warning('Phase_Aperiodic:NoRobustfit', '%s', robustFitFallbackNote);
end
if USE_ROBUSTFIT
    fprintf('robustfit() is available -- Method 2 uses robust (outlier-downweighted) linear regression.\n');
else
    fprintf('*** robustfit() UNAVAILABLE -- Method 2 falls back to plain least-squares polyfit() (see warning above). ***\n');
end

VARIANTS = {'MeanPSDBased'};
if hasMedianVariant
    VARIANTS{end+1} = 'MedianPSDBased';
end
METHODS = {'FOOOF', 'RobustFit'};
METRIC_TYPES = {'AperiodicExponent', 'AperiodicOffset'};
CONDITION_ORDER = {'E', 'N', 'R'};
REF_LEVEL = CONDITION_ORDER{1};
FDR_ALPHA = 0.05;

TASK_TYPES = {'intervention'};
TASK_DISPLAY = struct('intervention', 'Intervention');

fprintf('\n');

%% ========================================================================
for taskI = 1:numel(TASK_TYPES)
    TASK_KEY = TASK_TYPES{taskI};
    TASK_LABEL = TASK_DISPLAY.(TASK_KEY);

    STATS_DIR = fullfile(ANALYSIS_ROOT, '15_aperiodic_1f_family', RUN_MODE);
    if ~exist(STATS_DIR, 'dir'); mkdir(STATS_DIR); end

    fprintf('\n########################################################\n');
    fprintf('### TASK = %s | RUN_MODE = %s\n', TASK_LABEL, RUN_MODE_LABEL);
    fprintf('########################################################\n\n');

    integrityIssues = {};

    %% --------------------------------------------------------------
    taskMaskAll = strcmp({psd_data.TaskType}, TASK_KEY);
    thisTaskEntries = psd_data(taskMaskAll);
    nBaselineExcluded = sum(strcmp({thisTaskEntries.Condition}, 'Baseline'));
    entries = thisTaskEntries(ismember({thisTaskEntries.Condition}, CONDITION_ORDER));
    [entries, RUN_MODE_N_EXCL_ROWS, RUN_MODE_EXCL_SUBJ] = apply_run_mode_filter(entries, RUN_MODE, EXCLUDE_NOISY_SUBJECTS);
    if RUN_MODE_N_EXCL_ROWS > 0; fprintf('RUN_MODE=%s: removed %d recording(s) of subject(s) {%s}.\n', RUN_MODE_LABEL, RUN_MODE_N_EXCL_ROWS, num2str(RUN_MODE_EXCL_SUBJ)); end
    fprintf('Task=%s: %d total recordings loaded, %d Baseline recordings excluded, %d E/N/R recordings retained.\n', ...
        TASK_LABEL, numel(thisTaskEntries), nBaselineExcluded, numel(entries));

    if isempty(entries)
        warning('Task=%s: no E/N/R recordings found after filtering -- skipping this Task entirely.', TASK_LABEL);
        continue;
    end

    nEntries = numel(entries);
    for e = 1:nEntries
        if numel(entries(e).FrequenciesHz) ~= numel(FrequencyVectorHz) || ...
                ~isequal(entries(e).FrequenciesHz(:), FrequencyVectorHz)
            error('Phase_Aperiodic:FrequencyVectorMismatch', ...
                'Recording Sub%s Task=%s Block%g Cond%s has a FrequenciesHz vector that does not match the shared FrequencyVectorHz.', ...
                entries(e).Subject, TASK_LABEL, entries(e).BlockNumber, entries(e).Condition);
        end
        if numel(entries(e).ChannelLabels) ~= nChan || ~isequal(upper(entries(e).ChannelLabels(:)), upper(ChannelLabelsReference))
            error('Phase_Aperiodic:ChannelLabelMismatch', ...
                'Recording Sub%s Task=%s Block%g Cond%s has ChannelLabels that do not match ChannelLabelsReference.', ...
                entries(e).Subject, TASK_LABEL, entries(e).BlockNumber, entries(e).Condition);
        end
    end

    %% --------------------------------------------------------------
    RecMetaT = table({entries.Subject}', [entries.BlockNumber]', {entries.Condition}', ...
        'VariableNames', {'Subject', 'BlockNumber', 'Condition'});
    comboT = unique(RecMetaT, 'rows');
    nRecordings = height(comboT);
    fprintf('Section 3: %d unique Subject x BlockNumber x Condition recordings found.\n', nRecordings);

    subjNumeric = local_subject_numeric(comboT.Subject);
    uniqSubj = unique(subjNumeric);
    nSubjects = numel(uniqSubj);

    subjPairs = cell(nSubjects, 1);
    pairTypeCount = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for i = 1:nSubjects
        rowsI = subjNumeric == uniqSubj(i);
        condsI = local_to_cellstr(comboT.Condition(rowsI));
        if numel(condsI) ~= 2
            integrityIssues{end+1} = sprintf('Subject %g has %d recordings (expected 2) for Task=%s.', uniqSubj(i), numel(condsI), TASK_LABEL); %#ok<AGROW>
        elseif strcmp(condsI{1}, condsI{2})
            integrityIssues{end+1} = sprintf('Subject %g has the SAME condition (%s) in both recordings for Task=%s -- not a valid pair.', uniqSubj(i), condsI{1}, TASK_LABEL); %#ok<AGROW>
        else
            pairKey = strjoin(sort(condsI), '_');
            if isKey(pairTypeCount, pairKey)
                pairTypeCount(pairKey) = pairTypeCount(pairKey) + 1;
            else
                pairTypeCount(pairKey) = 1;
            end
        end
        subjPairs{i} = condsI;
    end
    fprintf('Condition-pair coverage across %d subjects (Task=%s):\n', nSubjects, TASK_LABEL);
    pairKeys = keys(pairTypeCount);
    for k = 1:numel(pairKeys)
        fprintf('  %s : %d subjects\n', strrep(pairKeys{k}, '_', ' & '), pairTypeCount(pairKeys{k}));
    end

    allSubjectsComplete3 = false;
    condsPerSubjCounts = cellfun(@numel, subjPairs);
    if all(condsPerSubjCounts >= 3)
        allSubjectsComplete3 = true;
        integrityIssues{end+1} = sprintf(['Task=%s: every subject appears to have 3+ recordings/conditions -- the incomplete-design ' ...
            'justification for using an LMM instead of RM-ANOVA no longer applies. STOP and re-examine.'], TASK_LABEL); %#ok<AGROW>
    end
    if allSubjectsComplete3
        error('Phase_Aperiodic:DesignAssumptionViolated', ...
            'Task=%s: design verification found a complete 3-condition-per-subject structure. Stopping before fitting -- see integrityIssues.', TASK_LABEL);
    end

    %% --------------------------------------------------------------
    %  SECTION 4: BUILD CHANNEL-LEVEL POWER ROWS (1-45 Hz), COMPUTE
    %     METHOD 2 (RobustFit) DIRECTLY, AND STAGE METHOD 1's CSV INPUT
    % --------------------------------------------------------------
    fprintf('\nSection 4: building %d recordings x %d channels = %d channel-level rows (Task=%s) ...\n', ...
        nEntries, nChan, nEntries * nChan, TASK_LABEL);

    nRows = nEntries * nChan;
    SubjectCol   = cell(nRows, 1);
    TaskTypeCol  = cell(nRows, 1);
    BlockNumberCol = zeros(nRows, 1);
    ConditionCol = cell(nRows, 1);
    ChannelCol   = cell(nRows, 1);
    RowIDCol     = (1:nRows)';

    PowerMean45 = nan(nRows, nFreq45);
    RobustExp_Mean = nan(nRows, 1);
    RobustOff_Mean = nan(nRows, 1);
    RobustStatus_Mean = cell(nRows, 1);

    if hasMedianVariant
        PowerMedian45 = nan(nRows, nFreq45);
        RobustExp_Median = nan(nRows, 1);
        RobustOff_Median = nan(nRows, 1);
        RobustStatus_Median = cell(nRows, 1);
    end

    nZeroValidWindowRecordings = 0;
    rowCounter = 0;
    for e = 1:nEntries
        entry = entries(e);
        powerMean45_e = entry.MeanPSD(:, freqMask45);   % nChan x nFreq45

        if hasMedianVariant
            validW = entry.WindowValid;
            if any(validW)
                medPSD = squeeze(median(entry.PSD_WindowLevel(validW, :, :), 1));  % nChan x nFreqSel
                if size(medPSD, 1) ~= nChan
                    medPSD = medPSD';   % defend against squeeze collapsing a singleton channel dim
                end
                powerMedian45_e = medPSD(:, freqMask45);
            else
                nZeroValidWindowRecordings = nZeroValidWindowRecordings + 1;
                powerMedian45_e = nan(nChan, nFreq45);
            end
        end

        for c = 1:nChan
            rowCounter = rowCounter + 1;
            SubjectCol{rowCounter} = entry.Subject;
            TaskTypeCol{rowCounter} = entry.TaskType;
            BlockNumberCol(rowCounter) = entry.BlockNumber;
            ConditionCol{rowCounter} = entry.Condition;
            ChannelCol{rowCounter} = ChannelLabelsReference{c};

            PowerMean45(rowCounter, :) = powerMean45_e(c, :);
            [slopeM, offM, statM] = local_fit_robust_1f(freqVec45, powerMean45_e(c, :), USE_ROBUSTFIT);
            RobustExp_Mean(rowCounter) = -slopeM;
            RobustOff_Mean(rowCounter) = offM;
            RobustStatus_Mean{rowCounter} = statM;

            if hasMedianVariant
                PowerMedian45(rowCounter, :) = powerMedian45_e(c, :);
                [slopeD, offD, statD] = local_fit_robust_1f(freqVec45, powerMedian45_e(c, :), USE_ROBUSTFIT);
                RobustExp_Median(rowCounter) = -slopeD;
                RobustOff_Median(rowCounter) = offD;
                RobustStatus_Median{rowCounter} = statD;
            end
        end
    end
    if nZeroValidWindowRecordings > 0
        integrityIssues{end+1} = sprintf('%d recording(s) had 0 valid windows -- MedianPSDBased rows for those recordings are NaN (excluded list-wise downstream), never fabricated.', nZeroValidWindowRecordings); %#ok<AGROW>
    end
    fprintf('Section 4 complete (Task=%s): %d channel-level rows built; Method 2 (RobustFit) computed for MeanPSDBased%s.\n', ...
        TASK_LABEL, nRows, local_yn_str(hasMedianVariant, ' and MedianPSDBased', ''));

    %% --------------------------------------------------------------
    %  SECTION 5: METHOD 1 (FOOOF) -- WRITE CSV, CALL PYTHON, READ BACK
    % --------------------------------------------------------------
    fprintf('\nSection 5: Method 1 (FOOOF) round-trip via Python (Task=%s) ...\n', TASK_LABEL);

    FOOOF_RESULTS = struct();   % FOOOF_RESULTS.(VariantName) = table with RowID, exponent, offset, r2, status

    for vI = 1:numel(VARIANTS)
        variantName = VARIANTS{vI};
        if strcmp(variantName, 'MeanPSDBased')
            powerBlock = PowerMean45;
        else
            powerBlock = PowerMedian45;
        end

        inputCsvPath = fullfile(STATS_DIR, sprintf('fooof_input_%s_%s.csv', TASK_LABEL, variantName));
        outputCsvPath = fullfile(STATS_DIR, sprintf('fooof_output_%s_%s.csv', TASK_LABEL, variantName));

        local_write_wide_psd_csv(inputCsvPath, SubjectCol, TaskTypeCol, BlockNumberCol, ...
            ConditionCol, ChannelCol, RowIDCol, freqVec45, powerBlock);
        fprintf('  [%s / %s] Wrote FOOOF input CSV: %s (%d rows x %d freq bins)\n', ...
            TASK_LABEL, variantName, inputCsvPath, nRows, nFreq45);

        cmd = sprintf('"%s" "%s" --input "%s" --output "%s" --freq-range %g %g', ...
            PYTHON_EXECUTABLE, PYTHON_SCRIPT_PATH, inputCsvPath, outputCsvPath, freqVec45(1), freqVec45(end));
        fprintf('  [%s / %s] Calling: %s\n', TASK_LABEL, variantName, cmd);
        [sysStatus, sysOut] = system(cmd);
        fprintf('%s\n', sysOut);
        if sysStatus ~= 0
            error('Phase_Aperiodic:FOOOFFailed', ...
                ['FOOOF fitting failed -- ensure Python 3 is installed and `pip install fooof` ' ...
                 'has been run. Command: %s. Exit status: %d. Command output:\n%s'], ...
                cmd, sysStatus, sysOut);
        end
        if ~exist(outputCsvPath, 'file')
            error('Phase_Aperiodic:FOOOFOutputMissing', ...
                ['FOOOF process reported success (exit status 0) but did not write the expected ' ...
                 'output file: %s. Command: %s'], outputCsvPath, cmd);
        end

        outT = readtable(outputCsvPath);
        requiredOutCols = {'RowID', 'AperiodicExponent_FOOOF', 'AperiodicOffset_FOOOF', 'FOOOF_RSquared', 'FitStatus'};
        missingOutCols = requiredOutCols(~ismember(requiredOutCols, outT.Properties.VariableNames));
        if ~isempty(missingOutCols)
            error('Phase_Aperiodic:FOOOFOutputColumnsMissing', ...
                'FOOOF output CSV %s is missing expected column(s): %s', outputCsvPath, strjoin(missingOutCols, ', '));
        end
        if height(outT) ~= nRows
            error('Phase_Aperiodic:FOOOFOutputRowCountMismatch', ...
                'FOOOF output CSV %s has %d rows; expected %d (one per input row). It must never silently drop rows.', ...
                outputCsvPath, height(outT), nRows);
        end
        outT = sortrows(outT, 'RowID');
        if ~isequal(outT.RowID, RowIDCol)
            error('Phase_Aperiodic:FOOOFOutputRowIDMismatch', ...
                'FOOOF output CSV %s RowID values do not match the input RowID sequence after sorting.', outputCsvPath);
        end
        FOOOF_RESULTS.(variantName) = outT;

        nFOOOFOK = sum(strcmp(local_to_cellstr(outT.FitStatus), 'OK'));
        nFOOOFFailed = height(outT) - nFOOOFOK;
        fprintf('  [%s / %s] FOOOF results read back: %d/%d OK, %d/%d FAILED.\n', ...
            TASK_LABEL, variantName, nFOOOFOK, height(outT), nFOOOFFailed, height(outT));
        if nFOOOFFailed > 0
            integrityIssues{end+1} = sprintf('%s / %s: %d/%d FOOOF fits FAILED (see fooof_output_%s_%s.csv FitStatus column and this Task''s stdout log above for per-row reasons; those rows are NaN in FeatureT, never dropped).', ...
                TASK_LABEL, variantName, nFOOOFFailed, height(outT), TASK_LABEL, variantName); %#ok<AGROW>
        end
    end

    %% --------------------------------------------------------------
    %  SECTION 6: ASSEMBLE THE CHANNEL-LEVEL LONG FEATURE TABLE
    %     (Method x Variant stacked; Level='Electrode')
    % --------------------------------------------------------------
    subTables = {};
    for vI = 1:numel(VARIANTS)
        variantName = VARIANTS{vI};
        outT = FOOOF_RESULTS.(variantName);

        % --- FOOOF rows ---
        tF = table(repmat({'Electrode'}, nRows, 1), ChannelCol, SubjectCol, TaskTypeCol, ...
            BlockNumberCol, ConditionCol, repmat({'FOOOF'}, nRows, 1), repmat({variantName}, nRows, 1), ...
            outT.AperiodicExponent_FOOOF, outT.AperiodicOffset_FOOOF, outT.FOOOF_RSquared, local_to_cellstr(outT.FitStatus), ...
            repmat(1, nRows, 1), ...
            'VariableNames', {'Level','Unit','Subject','TaskType','BlockNumber','Condition','Method','Variant', ...
                               'AperiodicExponent','AperiodicOffset','RSquared','FitStatus','NumChannelsAveraged'});
        subTables{end+1} = tF; %#ok<AGROW>

        % --- RobustFit rows ---
        if strcmp(variantName, 'MeanPSDBased')
            robExp = RobustExp_Mean; robOff = RobustOff_Mean; robStat = RobustStatus_Mean;
        else
            robExp = RobustExp_Median; robOff = RobustOff_Median; robStat = RobustStatus_Median;
        end
        tR = table(repmat({'Electrode'}, nRows, 1), ChannelCol, SubjectCol, TaskTypeCol, ...
            BlockNumberCol, ConditionCol, repmat({'RobustFit'}, nRows, 1), repmat({variantName}, nRows, 1), ...
            robExp, robOff, nan(nRows, 1), robStat, repmat(1, nRows, 1), ...
            'VariableNames', {'Level','Unit','Subject','TaskType','BlockNumber','Condition','Method','Variant', ...
                               'AperiodicExponent','AperiodicOffset','RSquared','FitStatus','NumChannelsAveraged'});
        subTables{end+1} = tR; %#ok<AGROW>
    end
    FeatureT_Electrode = vertcat(subTables{:});
    FeatureT_Electrode.ChannelUC = upper(FeatureT_Electrode.Unit);

    fprintf('\nSection 6 complete (Task=%s): FeatureT_Electrode has %d rows (%d channel-rows x %d methods x %d variants).\n', ...
        TASK_LABEL, height(FeatureT_Electrode), nRows, numel(METHODS), numel(VARIANTS));

    %% --------------------------------------------------------------
    %  SECTION 7: AGGREGATE TO ROI AND GLOBAL LEVELS (arithmetic mean of
    %     already-fit per-channel Exponent/Offset)
    % --------------------------------------------------------------
    recCombos2 = unique(FeatureT_Electrode(:, {'Subject','BlockNumber','Condition','Method','Variant'}), 'rows');
    roiRows = cell(0, 8);      % Subject,BlockNumber,Condition,ROI,Method,Variant,Exp,Off,nAvg(added separately)
    globalRows = cell(0, 7);   % Subject,BlockNumber,Condition,Method,Variant,Exp,Off
    roiNAvg = [];
    globalNAvg = [];

    for i = 1:height(recCombos2)
        mask = strcmp(FeatureT_Electrode.Subject, recCombos2.Subject{i}) & ...
               FeatureT_Electrode.BlockNumber == recCombos2.BlockNumber(i) & ...
               strcmp(FeatureT_Electrode.Condition, recCombos2.Condition{i}) & ...
               strcmp(FeatureT_Electrode.Method, recCombos2.Method{i}) & ...
               strcmp(FeatureT_Electrode.Variant, recCombos2.Variant{i});
        subChan = FeatureT_Electrode(mask, :);

        % --- Global: all 62 reference channels (P9/P10 INCLUDED, matching
        %     X05's whole-scalp convention). ---
        globExp = mean(subChan.AperiodicExponent, 'omitnan');
        globOff = mean(subChan.AperiodicOffset, 'omitnan');
        nAvgGlob = sum(~isnan(subChan.AperiodicExponent));
        globalRows(end+1, :) = {recCombos2.Subject{i}, recCombos2.BlockNumber(i), recCombos2.Condition{i}, ...
            recCombos2.Method{i}, recCombos2.Variant{i}, globExp, globOff}; %#ok<AGROW>
        globalNAvg(end+1, 1) = nAvgGlob; %#ok<AGROW>

        % --- ROI: locked ROI_DEFINITION member channels only (P9/P10
        %     EXCLUDED, matching X07/X08's ROI-averaging convention). ---
        for r = 1:nROI
            roiChansUC = upper(ROI_DEFINITION.(ROI_NAMES{r}));
            roiMask = ismember(subChan.ChannelUC, roiChansUC);
            roiExp = mean(subChan.AperiodicExponent(roiMask), 'omitnan');
            roiOff = mean(subChan.AperiodicOffset(roiMask), 'omitnan');
            nAvgRoi = sum(roiMask & ~isnan(subChan.AperiodicExponent));
            roiRows(end+1, :) = {recCombos2.Subject{i}, recCombos2.BlockNumber(i), recCombos2.Condition{i}, ...
                ROI_NAMES{r}, recCombos2.Method{i}, recCombos2.Variant{i}, roiExp, roiOff}; %#ok<AGROW>
            roiNAvg(end+1, 1) = nAvgRoi; %#ok<AGROW>
        end
    end

    GlobalT = cell2table(globalRows, 'VariableNames', {'Subject','BlockNumber','Condition','Method','Variant','AperiodicExponent','AperiodicOffset'});
    GlobalT.BlockNumber = cell2mat(cellfun(@(x) double(x), globalRows(:,2), 'UniformOutput', false));
    GlobalT.AperiodicExponent = cell2mat(globalRows(:,6));
    GlobalT.AperiodicOffset = cell2mat(globalRows(:,7));

    RoiT = cell2table(roiRows, 'VariableNames', {'Subject','BlockNumber','Condition','ROI','Method','Variant','AperiodicExponent','AperiodicOffset'});
    RoiT.BlockNumber = cell2mat(cellfun(@(x) double(x), roiRows(:,2), 'UniformOutput', false));
    RoiT.AperiodicExponent = cell2mat(roiRows(:,7));
    RoiT.AperiodicOffset = cell2mat(roiRows(:,8));

    nGlobalRows = height(GlobalT);
    nRoiRows = height(RoiT);

    FeatureT_Global = table(repmat({'Global'}, nGlobalRows, 1), repmat({'Global'}, nGlobalRows, 1), ...
        GlobalT.Subject, repmat({TASK_LABEL}, nGlobalRows, 1), GlobalT.BlockNumber, GlobalT.Condition, ...
        GlobalT.Method, GlobalT.Variant, GlobalT.AperiodicExponent, GlobalT.AperiodicOffset, ...
        nan(nGlobalRows, 1), repmat({''}, nGlobalRows, 1), globalNAvg, ...
        'VariableNames', {'Level','Unit','Subject','TaskType','BlockNumber','Condition','Method','Variant', ...
                           'AperiodicExponent','AperiodicOffset','RSquared','FitStatus','NumChannelsAveraged'});

    FeatureT_ROI = table(repmat({'ROI'}, nRoiRows, 1), RoiT.ROI, RoiT.Subject, repmat({TASK_LABEL}, nRoiRows, 1), ...
        RoiT.BlockNumber, RoiT.Condition, RoiT.Method, RoiT.Variant, RoiT.AperiodicExponent, RoiT.AperiodicOffset, ...
        nan(nRoiRows, 1), repmat({''}, nRoiRows, 1), roiNAvg, ...
        'VariableNames', {'Level','Unit','Subject','TaskType','BlockNumber','Condition','Method','Variant', ...
                           'AperiodicExponent','AperiodicOffset','RSquared','FitStatus','NumChannelsAveraged'});

    FeatureT_Electrode.TaskType = repmat({TASK_LABEL}, height(FeatureT_Electrode), 1);
    FeatureT = [FeatureT_Electrode(:, {'Level','Unit','Subject','TaskType','BlockNumber','Condition','Method','Variant', ...
                                        'AperiodicExponent','AperiodicOffset','RSquared','FitStatus','NumChannelsAveraged'}); ...
                FeatureT_ROI; FeatureT_Global];

    fprintf('\nSection 7 complete (Task=%s): FeatureT has %d rows total (%d Electrode + %d ROI + %d Global).\n', ...
        TASK_LABEL, height(FeatureT), height(FeatureT_Electrode), nRoiRows, nGlobalRows);

    %% --------------------------------------------------------------
    %  SECTION 8: LMM FITTING -- Level x Unit x Method x Variant x
    %     MetricType (Metric ~ Condition + (1|Subject), REML, Satterthwaite)
    % --------------------------------------------------------------
    LEVELS = {'Electrode', 'ROI', 'Global'};
    unitsByLevel = struct('Electrode', {ChannelLabelsReference}, 'ROI', {ROI_NAMES(:)}, 'Global', {{'Global'}});

    nExpectedFits = (nChan + nROI + 1) * numel(METHODS) * numel(VARIANTS) * numel(METRIC_TYPES);
    fprintf('\nSection 8: fitting up to %d LMMs (Task=%s: %d Electrode + %d ROI + %d Global units, x %d Methods x %d Variants x %d MetricTypes) ...\n', ...
        nExpectedFits, TASK_LABEL, nChan, nROI, 1, numel(METHODS), numel(VARIANTS), numel(METRIC_TYPES));

    omnibusRows = cell(0, 9);
    lmeMap = containers.Map();
    convergenceWarnings = {};

    for lv = 1:numel(LEVELS)
        level = LEVELS{lv};
        units = unitsByLevel.(level);
        for u = 1:numel(units)
            unitName = units{u};
            for mI = 1:numel(METHODS)
                method = METHODS{mI};
                for vI = 1:numel(VARIANTS)
                    variantName = VARIANTS{vI};
                    for mtI = 1:numel(METRIC_TYPES)
                        metricType = METRIC_TYPES{mtI};

                        subT = FeatureT(strcmp(FeatureT.Level, level) & strcmp(FeatureT.Unit, unitName) & ...
                            strcmp(FeatureT.Method, method) & strcmp(FeatureT.Variant, variantName), :);
                        subT = subT(~isnan(subT.(metricType)), :);
                        if height(subT) < 4
                            integrityIssues{end+1} = sprintf('%s / Unit=%s / %s / %s / %s (Task=%s): only %d non-NaN rows -- LMM not fit.', ...
                                level, unitName, method, variantName, metricType, TASK_LABEL, height(subT)); %#ok<AGROW>
                            continue;
                        end
                        subT.Subject = categorical(local_subject_numeric(subT.Subject));
                        subT.Condition = categorical(local_to_cellstr(subT.Condition), CONDITION_ORDER);
                        formula = sprintf('%s ~ Condition + (1|Subject)', metricType);
                        try
                            lastwarn('');
                            lme = fitlme(subT, formula, 'FitMethod', 'REML');
                            [warnMsg, ~] = lastwarn();
                            if ~isempty(warnMsg)
                                convergenceWarnings{end+1} = sprintf('%s / %s / %s / %s / %s (Task=%s): %s', ...
                                    level, unitName, method, variantName, metricType, TASK_LABEL, warnMsg); %#ok<AGROW>
                            end
                        catch ME
                            integrityIssues{end+1} = sprintf('fitlme FAILED for %s / Unit=%s / %s / %s / %s (Task=%s): %s', ...
                                level, unitName, method, variantName, metricType, TASK_LABEL, ME.message); %#ok<AGROW>
                            continue;
                        end
                        aov = anova(lme, 'DFMethod', 'satterthwaite');
                        rowIdx = find(strcmp(aov.Term, 'Condition'), 1);
                        if isempty(rowIdx)
                            integrityIssues{end+1} = sprintf('anova(lme) has no Condition row for %s / Unit=%s / %s / %s / %s (Task=%s).', ...
                                level, unitName, method, variantName, metricType, TASK_LABEL); %#ok<AGROW>
                            continue;
                        end
                        omnibusRows(end+1, :) = {level, unitName, method, variantName, metricType, ...
                            aov.FStat(rowIdx), aov.DF1(rowIdx), aov.DF2(rowIdx), aov.pValue(rowIdx)}; %#ok<AGROW>
                        key = strjoin({level, unitName, method, variantName, metricType}, '|');
                        lmeMap(key) = lme;
                    end
                end
            end
        end
    end

    OmnibusT = local_build_table(omnibusRows, [1 2 3 4 5], [6 7 8 9], ...
        {'Level','Unit','Method','Variant','MetricType','FStat','DF1','DF2','pValue'});
    fprintf('Section 8 complete (Task=%s): %d / %d omnibus tests fit successfully.\n', TASK_LABEL, height(OmnibusT), nExpectedFits);

    %% --------------------------------------------------------------
    %  SECTION 9: BH-FDR -- separately within each (Level x Method x
    %     Variant x MetricType) family, across that Level's Units
    % --------------------------------------------------------------
    OmnibusT.pFDR = nan(height(OmnibusT), 1);
    OmnibusT.SignificantFDR = false(height(OmnibusT), 1);
    OmnibusT.FamilySize = nan(height(OmnibusT), 1);

    familyKeyAll = strcat(OmnibusT.Level, '|', OmnibusT.Method, '|', OmnibusT.Variant, '|', OmnibusT.MetricType);
    familyKeys = unique(familyKeyAll);
    for k = 1:numel(familyKeys)
        idx = strcmp(familyKeyAll, familyKeys{k});
        pvals = OmnibusT.pValue(idx);
        OmnibusT.FamilySize(idx) = numel(pvals);
        if numel(pvals) == 1
            OmnibusT.pFDR(idx) = pvals;
            OmnibusT.SignificantFDR(idx) = pvals <= FDR_ALPHA;
        else
            [adjP, reject] = local_fdr_bh(pvals, FDR_ALPHA);
            OmnibusT.pFDR(idx) = adjP;
            OmnibusT.SignificantFDR(idx) = reject;
        end
    end
    fprintf('Section 9 complete (Task=%s): BH-FDR applied within %d independent (Level x Method x Variant x MetricType) families.\n', ...
        TASK_LABEL, numel(familyKeys));

    %% --------------------------------------------------------------
    %  SECTION 10: PAIRWISE CONTRASTS (E-N, E-R, N-R -- X05/X06's
    %     convention), gated on FDR, Holm-corrected (family size 3)
    % --------------------------------------------------------------
    pairwiseRows = cell(0, 14);
    CONTRASTS = {'E','N','E - N'; 'E','R','E - R'; 'N','R','N - R'};   % {level2, level1, label}; estimate = level2 - level1

    sigIdx = find(OmnibusT.SignificantFDR)';
    for row = sigIdx
        level = OmnibusT.Level{row}; unitName = OmnibusT.Unit{row};
        method = OmnibusT.Method{row}; variantName = OmnibusT.Variant{row}; metricType = OmnibusT.MetricType{row};
        key = strjoin({level, unitName, method, variantName, metricType}, '|');
        if ~isKey(lmeMap, key); continue; end
        lme = lmeMap(key);
        coeffNames = lme.CoefficientNames;
        beta = lme.Coefficients.Estimate;
        covB = lme.CoefficientCovariance;

        contrastP = nan(3, 1);
        contrastData = cell(3, 8);
        for c = 1:3
            lvl2 = CONTRASTS{c, 1}; lvl1 = CONTRASTS{c, 2}; label = CONTRASTS{c, 3};
            C = local_condition_contrast(coeffNames, 'Condition', REF_LEVEL, lvl1, lvl2);
            est = C * beta;
            se = sqrt(C * covB * C');
            [p, ~, ~, df2c] = coefTest(lme, C, 0, 'DFMethod', 'satterthwaite');
            tstat = est / se;
            tcrit = tinv(0.975, df2c);
            ciLow = est - tcrit * se;
            ciHigh = est + tcrit * se;
            contrastP(c) = p;
            contrastData(c, :) = {label, est, se, df2c, tstat, p, ciLow, ciHigh};
        end
        [adjP, ~] = local_holm_bonferroni(contrastP, FDR_ALPHA);
        for c = 1:3
            pairwiseRows(end+1, :) = [{level, unitName, method, variantName, metricType}, contrastData(c, 1:6), {adjP(c)}, contrastData(c, 7:8)]; %#ok<AGROW>
        end
    end
    PairwiseT = local_build_table(pairwiseRows, [1 2 3 4 5 6], [7 8 9 10 11 12 13 14], ...
        {'Level','Unit','Method','Variant','MetricType','Contrast','Estimate','SE','DF','tStat','pValue','pHolm','CI_Low','CI_High'});

    fprintf('Section 10 complete (Task=%s): %d pairwise contrast rows (bands/units significant under their own FDR family).\n', ...
        TASK_LABEL, height(PairwiseT));

    %% --------------------------------------------------------------
    %  SECTION 11: FOOOF vs. RobustFit AGREEMENT (sanity/consistency
    %     check), Electrode level, per available Variant
    % --------------------------------------------------------------
    agreementLines = {};
    for vI = 1:numel(VARIANTS)
        variantName = VARIANTS{vI};
        fT = FeatureT_Electrode(strcmp(FeatureT_Electrode.Method, 'FOOOF') & strcmp(FeatureT_Electrode.Variant, variantName), ...
            {'Subject','BlockNumber','Condition','Unit','AperiodicExponent','AperiodicOffset','FitStatus'});
        rT = FeatureT_Electrode(strcmp(FeatureT_Electrode.Method, 'RobustFit') & strcmp(FeatureT_Electrode.Variant, variantName), ...
            {'Subject','BlockNumber','Condition','Unit','AperiodicExponent','AperiodicOffset'});
        fT = renamevars(fT, {'AperiodicExponent','AperiodicOffset'}, {'Exp_FOOOF','Off_FOOOF'});
        rT = renamevars(rT, {'AperiodicExponent','AperiodicOffset'}, {'Exp_Robust','Off_Robust'});
        joinedT = innerjoin(fT, rT, 'Keys', {'Subject','BlockNumber','Condition','Unit'});
        okMask = strcmp(local_to_cellstr(joinedT.FitStatus), 'OK') & ~isnan(joinedT.Exp_FOOOF) & ~isnan(joinedT.Exp_Robust);
        nComparable = sum(okMask);
        if nComparable >= 3
            corrExp = corr(joinedT.Exp_FOOOF(okMask), joinedT.Exp_Robust(okMask));
            madExp = mean(abs(joinedT.Exp_FOOOF(okMask) - joinedT.Exp_Robust(okMask)));
            corrOff = corr(joinedT.Off_FOOOF(okMask), joinedT.Off_Robust(okMask));
            madOff = mean(abs(joinedT.Off_FOOOF(okMask) - joinedT.Off_Robust(okMask)));
            agreementLines{end+1} = sprintf('%s: Exponent r=%.3f, mean|diff|=%.4f (n=%d comparable channel-rows); Offset r=%.3f, mean|diff|=%.4f.', ...
                variantName, corrExp, madExp, nComparable, corrOff, madOff); %#ok<AGROW>
        else
            agreementLines{end+1} = sprintf('%s: too few comparable rows (n=%d, need >=3) to compute a FOOOF vs. RobustFit agreement statistic.', variantName, nComparable); %#ok<AGROW>
        end
    end
    fprintf('\nSection 11 (Task=%s) -- FOOOF vs. RobustFit agreement:\n', TASK_LABEL);
    for i = 1:numel(agreementLines)
        fprintf('  %s\n', agreementLines{i});
    end

    %% --------------------------------------------------------------
    nDiagSaved = 0;
    for mI = 1:numel(METHODS)
        method = METHODS{mI};
        for vI = 1:numel(VARIANTS)
            variantName = VARIANTS{vI};
            for mtI = 1:numel(METRIC_TYPES)
                metricType = METRIC_TYPES{mtI};
                key = strjoin({'Global', 'Global', method, variantName, metricType}, '|');
                if ~isKey(lmeMap, key); continue; end
                lme = lmeMap(key);
                try
                    res = residuals(lme, 'ResidualType', 'Pearson');
                    fig = figure('Visible', 'off', 'Position', [100 100 900 400]);
                    subplot(1,2,1);
                    histogram(res, 15);
                    title(sprintf('Residuals: Global %s (%s,%s), Task=%s', metricType, method, variantName, TASK_LABEL), 'Interpreter', 'none');
                    xlabel('Pearson residual'); ylabel('Count');
                    subplot(1,2,2);
                    qqplot(res);
                    title('Normal Q-Q');
                    outPng = fullfile(STATS_DIR, sprintf('diagnostic_residuals_Global_%s_%s_%s.png', method, variantName, metricType));
                    saveas(fig, outPng);
                    close(fig);
                    nDiagSaved = nDiagSaved + 1;
                catch ME
                    integrityIssues{end+1} = sprintf('Residual diagnostic failed for Global / %s / %s / %s (Task=%s): %s', method, variantName, metricType, TASK_LABEL, ME.message); %#ok<AGROW>
                end
            end
        end
    end
    fprintf('Section 12 complete (Task=%s): %d Global-level residual diagnostic plots saved.\n', TASK_LABEL, nDiagSaved);

    %% --------------------------------------------------------------
    %  SECTION 13: PARAMETERS AND OUTPUT
    % --------------------------------------------------------------
    PARAMS15 = struct();
    PARAMS15.ScriptName = 'X15_aperiodic_1f_family_intervention.m';
    PARAMS15.Task = TASK_LABEL;
    PARAMS15.GeneratedOn = datestr(now); %#ok<TNOW1,DATST>
    PARAMS15.AnalysisLevel = sprintf('Aperiodic exponent/offset (X03 intervention_psd.mat MeanPSD, and MedianPSD if available), Task=%s only, Condition in {E,N,R} only (Baseline excluded)', TASK_LABEL);
    PARAMS15.NewAdditiveScriptNote = 'Purely additive exploratory extension (X11-X15 family); never reads or modifies any X05/X06/X08/X09/X10 statistical output; own FDR families, never pooled with the two locked co-primary metrics.';
    PARAMS15.FitRangeHz = [freqVec45(1), freqVec45(end)];
    PARAMS15.Model = 'LinearMixedEffects_RandomIntercept';
    PARAMS15.FixedEffectsPrimary = 'Condition (3-level: E, N, R; reference = E)';
    PARAMS15.RandomEffects = '(1|Subject)';
    PARAMS15.FitMethod = 'REML';
    PARAMS15.DFMethod = 'Satterthwaite';
    PARAMS15.Methods = METHODS;
    PARAMS15.Method1_FOOOF_Note = 'PRIMARY. Python fooof package, aperiodic_mode=fixed, peak_width_limits=[1 8], max_n_peaks=6, freq_range=[1 45]. Round-tripped via MATLAB system() call to X15_fit_fooof.py (NOT py.* interop). Hard error on nonzero exit status.';
    PARAMS15.Method2_RobustFit_Note = sprintf('SENSITIVITY/FALLBACK, MATLAB-native. robustfit(log10(freq),log10(power)) over the same 1-45 Hz range; Exponent = -slope, Offset = intercept, matching FOOOF''s power ~ 1/f^exponent convention. %s', ...
        local_yn_str(USE_ROBUSTFIT, 'robustfit() was available in this session.', robustFitFallbackNote));
    PARAMS15.Variants = VARIANTS;
    PARAMS15.MedianPSDBasedAvailable = hasMedianVariant;
    PARAMS15.MetricTypes = METRIC_TYPES;
    PARAMS15.Levels = LEVELS;
    PARAMS15.ROINames = ROI_NAMES;
    PARAMS15.ROIAggregationRule = 'Arithmetic mean of already-fit per-channel AperiodicExponent/AperiodicOffset across each ROI''s member channels (X07 locked ROI_DEFINITION; P9/P10 EXCLUDED, matching X07/X08). NOT a re-fit on an ROI-averaged PSD curve.';
    PARAMS15.GlobalAggregationRule = 'Arithmetic mean of already-fit per-channel AperiodicExponent/AperiodicOffset across all 62 reference channels (P9/P10 INCLUDED, matching X05''s whole-scalp convention). NOT a re-fit on a scalp-averaged PSD curve.';
    PARAMS15.FDRAlpha = FDR_ALPHA;
    PARAMS15.FDRFamilyScheme = 'Separate BH-FDR family per (Level x Method x Variant x MetricType), across that Level''s Units (Electrode family size=62, ROI family size=5, Global family size=1, reported as-is per X05''s Primary-tier convention). Never pooled with each other, and never pooled with X05/X06/X08/X09/X10''s FDR families.';
    PARAMS15.PairwiseContrastMethod = 'Holm-Bonferroni across the 3 pairwise contrasts (E-N, E-R, N-R -- matching X05/X06''s convention), computed only for Level x Unit x Method x Variant x MetricType combinations whose omnibus test survives its own family''s BH-FDR.';
    PARAMS15.CSVInterchangeFormat = 'WIDE (one row per Subject x TaskType x BlockNumber x Condition x Channel spectrum, one column per 1-45 Hz frequency bin, ~177 columns) -- chosen over long format because it keeps each row a single, directly FOOOF-fittable spectrum (no groupby/pivot needed on the Python side) and keeps per-Task row counts small (~nRecordings x 62 channels, a few thousand), while a long format at this scale (frequency bins x rows) would be far larger on disk with no benefit for a per-row/per-spectrum fitting workload.';
    PARAMS15.BaselineExclusionNote = sprintf('%d Baseline recording(s) excluded before any modeling; Baseline is processed upstream (X03) but never enters this or any other inferential model in this pipeline.', nBaselineExcluded);
    PARAMS15.SourcePhase2File = psdMatFile;
    PARAMS15.SourcePhase5DFile = roiDefFile;
    PARAMS15.PythonScriptPath = PYTHON_SCRIPT_PATH;
    PARAMS15.DesignVerification = struct('nRecordings', nRecordings, 'nSubjects', nSubjects, ...
        'PairTypeCounts', {pairTypeCount.keys}, 'IntegrityIssues', {integrityIssues});
    PARAMS15.RunDateTime = datestr(now); %#ok<TNOW1,DATST>

    save(fullfile(STATS_DIR, 'phase_aperiodic_parameters.mat'), 'PARAMS15');
    save(fullfile(STATS_DIR, 'phase_aperiodic_feature_table.mat'), 'FeatureT');
    save(fullfile(STATS_DIR, 'phase_aperiodic_omnibus_results.mat'), 'OmnibusT');
    save(fullfile(STATS_DIR, 'phase_aperiodic_pairwise_contrasts.mat'), 'PairwiseT');

    writetable(FeatureT, fullfile(STATS_DIR, 'phase_aperiodic_feature_table.xlsx'));
    writetable(FeatureT, fullfile(STATS_DIR, 'phase_aperiodic_feature_table.csv'));
    writetable(OmnibusT, fullfile(STATS_DIR, 'phase_aperiodic_omnibus_results.xlsx'));
    writetable(OmnibusT, fullfile(STATS_DIR, 'phase_aperiodic_omnibus_results.csv'));
    if ~isempty(PairwiseT)
        writetable(PairwiseT, fullfile(STATS_DIR, 'phase_aperiodic_pairwise_contrasts.xlsx'));
        writetable(PairwiseT, fullfile(STATS_DIR, 'phase_aperiodic_pairwise_contrasts.csv'));
    end

    fprintf('\nSaved output to: %s\n', STATS_DIR);

    %% --------------------------------------------------------------
    %  SECTION 14: REPORT
    % --------------------------------------------------------------
    reportLines = {};
    reportLines{end+1} = repmat('=', 1, 70);
    reportLines{end+1} = sprintf('PHASE "APERIODIC 1/F FAMILY" REPORT (Task=%s) -- X15, exploratory extension', TASK_LABEL);
    reportLines{end+1} = repmat('=', 1, 70);
    reportLines{end+1} = '';
    reportLines{end+1} = 'Exploratory extension (X11-X15 family). Never';
    reportLines{end+1} = 'reads or modifies any X05/X06/X08/X09/X10 statistical output; own FDR families only.';
    reportLines{end+1} = '';
    reportLines{end+1} = sprintf('Baseline recordings excluded before any modeling: %d', nBaselineExcluded);
    reportLines{end+1} = sprintf('Recordings verified: %d | Subjects: %d', nRecordings, nSubjects);
    reportLines{end+1} = sprintf('Condition-pair coverage: %s', strjoin(cellfun(@(k) sprintf('%s=%d', strrep(k,'_',' & '), pairTypeCount(k)), pairKeys, 'UniformOutput', false), ', '));
    reportLines{end+1} = sprintf('Fit range: %.4f-%.4f Hz (continuous 1-45 Hz retained portion only; 55-100 Hz never used).', freqVec45(1), freqVec45(end));
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Methods ---';
    reportLines{end+1} = 'Method 1 (PRIMARY): FOOOF (Python fooof package), aperiodic_mode=fixed, peak_width_limits=[1 8], max_n_peaks=6.';
    reportLines{end+1} = sprintf('Method 2 (sensitivity/fallback, MATLAB-native): %s', local_yn_str(USE_ROBUSTFIT, 'robustfit() (robust log-log linear regression).', ['polyfit() OLS fallback -- ' robustFitFallbackNote]));
    reportLines{end+1} = '';
    reportLines{end+1} = sprintf('MedianPSDBased variant: %s', local_yn_str(hasMedianVariant, 'AVAILABLE and computed alongside MeanPSDBased (psd_data carried PSD_WindowLevel).', 'SKIPPED cleanly (psd_data did not carry a per-window full-PSD field -- PARAMS.SaveFullWindowLevelPSD was false in X03, the default). Only MeanPSDBased was computed. This is NOT an error.'));
    reportLines{end+1} = '';
    reportLines{end+1} = '--- FOOOF fit success/failure counts (per Variant) ---';
    for vI = 1:numel(VARIANTS)
        variantName = VARIANTS{vI};
        outT = FOOOF_RESULTS.(variantName);
        nOK = sum(strcmp(local_to_cellstr(outT.FitStatus), 'OK'));
        reportLines{end+1} = sprintf('  %-16s: %d/%d OK, %d/%d FAILED', variantName, nOK, height(outT), height(outT)-nOK, height(outT)); %#ok<AGROW>
    end
    reportLines{end+1} = '';
    reportLines{end+1} = '--- FOOOF vs. RobustFit agreement (sanity/consistency check, Electrode level) ---';
    for i = 1:numel(agreementLines)
        reportLines{end+1} = ['  ' agreementLines{i}]; %#ok<AGROW>
    end
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Omnibus results (all Levels/Methods/Variants/MetricTypes; own FDR family per Level x Method x Variant x MetricType) ---';
    for r = 1:height(OmnibusT)
        reportLines{end+1} = sprintf('  [%-9s|%-9s] %-9s %-14s %-19s F(%g,%.1f)=%.3f p=%.4f pFDR=%.4f family_n=%d %s', ...
            OmnibusT.Level{r}, OmnibusT.Unit{r}, OmnibusT.Method{r}, OmnibusT.Variant{r}, OmnibusT.MetricType{r}, ...
            OmnibusT.DF1(r), OmnibusT.DF2(r), OmnibusT.FStat(r), OmnibusT.pValue(r), OmnibusT.pFDR(r), OmnibusT.FamilySize(r), ...
            local_yn(OmnibusT.SignificantFDR(r))); %#ok<AGROW>
    end
    reportLines{end+1} = '';
    reportLines{end+1} = '--- Pairwise contrasts (Level x Unit x Method x Variant x MetricType significant under its own FDR family; Holm-corrected, E-N/E-R/N-R convention) ---';
    if isempty(PairwiseT)
        reportLines{end+1} = '  (none -- nothing survived its family''s FDR correction)';
    else
        for r = 1:height(PairwiseT)
            reportLines{end+1} = sprintf('  [%-9s|%-9s] %-9s %-14s %-19s %-8s est=%.4f (95%% CI %.4f to %.4f) t=%.3f df=%.1f p=%.4f pHolm=%.4f', ...
                PairwiseT.Level{r}, PairwiseT.Unit{r}, PairwiseT.Method{r}, PairwiseT.Variant{r}, PairwiseT.MetricType{r}, PairwiseT.Contrast{r}, ...
                PairwiseT.Estimate(r), PairwiseT.CI_Low(r), PairwiseT.CI_High(r), PairwiseT.tStat(r), PairwiseT.DF(r), ...
                PairwiseT.pValue(r), PairwiseT.pHolm(r)); %#ok<AGROW>
        end
    end
    reportLines{end+1} = '';
    if ~isempty(integrityIssues)
        reportLines{end+1} = '--- INTEGRITY ISSUES (reported, not auto-corrected) ---';
        for i = 1:numel(integrityIssues)
            reportLines{end+1} = ['  - ' integrityIssues{i}]; %#ok<AGROW>
        end
    else
        reportLines{end+1} = 'Integrity issues: none.';
    end
    if ~isempty(convergenceWarnings)
        reportLines{end+1} = '';
        reportLines{end+1} = '--- Model convergence warnings ---';
        for i = 1:numel(convergenceWarnings)
            reportLines{end+1} = ['  - ' convergenceWarnings{i}]; %#ok<AGROW>
        end
    end
    reportLines{end+1} = '';
    reportLines{end+1} = 'Output files written to:';
    reportLines{end+1} = ['  ' STATS_DIR];
    reportLines{end+1} = repmat('=', 1, 70);

    reportText = strjoin(reportLines, '\n');
    fprintf('%s\n', reportText);
    fid = fopen(fullfile(REPORT_DIR, sprintf('phase_aperiodic_1f_family_report_%s_%s.txt', TASK_LABEL, RUN_MODE_LABEL)), 'w');
    fprintf(fid, '%s\n', reportText);
    fclose(fid);

end % TASK_TYPES loop

fprintf('\n=== X15 complete for all three Tasks. Do not proceed until every Task''s output has been reviewed. ===\n');

%% ======================= LOCAL FUNCTIONS (all after top-level code) =======================
% local_subject_numeric, local_fdr_bh, local_holm_bonferroni,
% local_condition_contrast, local_build_table, local_yn, local_fmt below
% are identical to those in X05_primary_statistics.m.
% local_fit_robust_1f, local_write_wide_psd_csv and local_yn_str are
% specific to this script.

function n = local_subject_numeric(subjVal)
    if iscategorical(subjVal)
        subjVal = cellstr(subjVal);
    end
    if iscell(subjVal)
        n = cellfun(@(s) str2double(regexprep(char(s), '\D', '')), subjVal);
    elseif isnumeric(subjVal)
        n = double(subjVal);
    else
        n = str2double(regexprep(char(subjVal), '\D', ''));
    end
end

function [adjP, reject] = local_fdr_bh(pvals, alphaFDR)
    pvals = pvals(:);
    m = numel(pvals);
    [sortedP, sortIdx] = sort(pvals);
    adjSorted = nan(m, 1);
    runningMin = 1;
    for i = m:-1:1
        val = sortedP(i) * m / i;
        runningMin = min(runningMin, val);
        adjSorted(i) = runningMin;
    end
    adjP = nan(m, 1);
    adjP(sortIdx) = adjSorted;
    adjP = min(adjP, 1);
    reject = adjP <= alphaFDR;
end

function [adjP, reject] = local_holm_bonferroni(pvals, alphaFDR)
    pvals = pvals(:);
    m = numel(pvals);
    [sortedP, sortIdx] = sort(pvals);
    adjSorted = nan(m, 1);
    runningMax = 0;
    for i = 1:m
        val = sortedP(i) * (m - i + 1);
        runningMax = max(runningMax, val);
        adjSorted(i) = runningMax;
    end
    adjP = nan(m, 1);
    adjP(sortIdx) = adjSorted;
    adjP = min(adjP, 1);
    reject = adjP <= alphaFDR;
end

function cRow = local_condition_contrast(coeffNames, varName, refLevel, level1, level2)
    cRow = zeros(1, numel(coeffNames));
    if ~strcmp(level1, refLevel)
        name1 = sprintf('%s_%s', varName, level1);
        idx1 = find(strcmp(coeffNames, name1));
        if isempty(idx1)
            error('Phase_Aperiodic:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name1);
        end
        cRow(idx1) = cRow(idx1) - 1;
    end
    if ~strcmp(level2, refLevel)
        name2 = sprintf('%s_%s', varName, level2);
        idx2 = find(strcmp(coeffNames, name2));
        if isempty(idx2)
            error('Phase_Aperiodic:ContrastNameNotFound', 'Coefficient name "%s" not found in model.', name2);
        end
        cRow(idx2) = cRow(idx2) + 1;
    end
end

function s = local_yn(tf)
    if tf
        s = 'SIGNIFICANT';
    else
        s = 'n.s.';
    end
end

function T = local_build_table(rows, textColIdx, numColIdx, varNames)
    nCols = numel(varNames);
    if size(rows, 2) ~= nCols
        error('Phase_Aperiodic:BuildTableColumnMismatch', ...
            'local_build_table: rows has %d columns but %d variable names were given.', size(rows, 2), nCols);
    end
    cols = cell(1, nCols);
    for j = 1:nCols
        if ismember(j, textColIdx)
            cols{j} = rows(:, j);
        elseif ismember(j, numColIdx)
            if isempty(rows)
                cols{j} = zeros(0, 1);
            else
                cols{j} = cell2mat(rows(:, j));
            end
        else
            error('Phase_Aperiodic:BuildTableColumnMismatch', 'Column %d not classified as text or numeric.', j);
        end
    end
    T = table(cols{:}, 'VariableNames', varNames);
end

function s = local_fmt(v)
    if isnan(v)
        s = 'NaN';
    else
        s = sprintf('%.4f', v);
    end
end

% ------------------------------------------------------------------
% NEW LOCAL FUNCTIONS (specific to X15)
% ------------------------------------------------------------------

function [slope, offsetVal, status] = local_fit_robust_1f(freqHz, powerVec, useRobustfit)
% Method 2: log-log linear fit of power ~ 1/f^exponent, i.e.
%   log10(power) = offset - exponent*log10(freq)  =>  slope = -exponent.
% Caller converts: exponent = -slope, offset = offsetVal. Uses
% robustfit() (outlier-downweighted) when available/successful for THIS
% row; falls back to plain polyfit() (ordinary least squares) otherwise,
% per-row, and marks the row's status accordingly ('OK' vs 'OK_LSQ').
% Never fabricates a fit from too few valid frequency bins -- returns
% NaN/'FAILED' instead.
    freqHz = freqHz(:);
    powerVec = powerVec(:);
    valid = isfinite(freqHz) & isfinite(powerVec) & powerVec > 0;
    if sum(valid) < 2
        slope = NaN; offsetVal = NaN; status = 'FAILED';
        return;
    end
    x = log10(freqHz(valid));
    y = log10(powerVec(valid));

    if useRobustfit
        try
            b = robustfit(x, y);
            offsetVal = b(1);
            slope = b(2);
            status = 'OK';
            return;
        catch
            % Fall through to the OLS fallback for THIS row only -- a
            % single row's robustfit failure (e.g. rank-deficient
            % weighting) must not be treated as toolbox-unavailable.
        end
    end
    try
        p = polyfit(x, y, 1);
        slope = p(1);
        offsetVal = p(2);
        status = 'OK_LSQ';
    catch
        slope = NaN; offsetVal = NaN; status = 'FAILED';
    end
end

function local_write_wide_psd_csv(csvPath, SubjectCol, TaskTypeCol, BlockNumberCol, ConditionCol, ChannelCol, RowIDCol, freqVec, powerMatrix)
    nRows = numel(SubjectCol);
    nFreq = numel(freqVec);
    if size(powerMatrix, 1) ~= nRows || size(powerMatrix, 2) ~= nFreq
        error('Phase_Aperiodic:CSVWriteDimMismatch', ...
            'local_write_wide_psd_csv: powerMatrix is %dx%d, expected %dx%d.', ...
            size(powerMatrix,1), size(powerMatrix,2), nRows, nFreq);
    end

    fid = fopen(csvPath, 'w');
    if fid < 0
        error('Phase_Aperiodic:CSVWriteFailed', 'Could not open %s for writing.', csvPath);
    end
    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    freqColNames = arrayfun(@(f) sprintf('F_%.4f', f), freqVec, 'UniformOutput', false);
    fprintf(fid, 'Subject,TaskType,BlockNumber,Condition,Channel,RowID,%s\n', strjoin(freqColNames, ','));

    rowFmt = ['%s,%s,%d,%s,%s,%d' repmat(',%.8g', 1, nFreq) '\n'];
    for r = 1:nRows
        fprintf(fid, rowFmt, SubjectCol{r}, TaskTypeCol{r}, BlockNumberCol(r), ConditionCol{r}, ChannelCol{r}, RowIDCol(r), powerMatrix(r, :));
    end
end

function s = local_yn_str(tf, sTrue, sFalse)
    if tf
        s = sTrue;
    else
        s = sFalse;
    end
end

function c = local_to_cellstr(x)
% Robustly converts char/string/categorical/cellstr input to a cell
% array of character vectors -- unlike MATLAB's own cellstr(), this
% does NOT error when x is already a cell array of char (which every
% text column this script builds itself, e.g. ConditionCol, is), and
% it also handles a table text column read back by readtable(), which
% may come back as char, string, categorical, or already cellstr
% depending on MATLAB version/import settings. Used everywhere this
% script needs a guaranteed cell array of char for strcmp()/categorical().
    if iscell(x)
        c = x;
    elseif iscategorical(x)
        c = cellstr(x);
    elseif isstring(x)
        c = cellstr(x);
    elseif ischar(x)
        c = cellstr(x);
    else
        c = cellstr(string(x));
    end
end
