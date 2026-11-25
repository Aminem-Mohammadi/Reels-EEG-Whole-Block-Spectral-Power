%% SART behavioral analysis
% -------------------------------------------------------------------------
% Computes, for every subject x block SART file:
%   - Commission errors (response to 3) and omission errors (no response to Go)
%   - Signal-detection indices d' and c (log-linear correction)
%   - Correct-Go RT: mean, median, SD, CV (= SD / mean)
%   - Commission (false-alarm) RT: mean, median, and Go-minus-Commission RT
%   - Ex-Gaussian parameters of correct-Go RT (mu, sigma, tau), MLE fit
%   - Pre-error speeding (RT in the 4 trials before commission vs before
%     correct withhold) and post-error slowing (traditional and robust)
%   - Time-on-task (quartiles) and font-size breakdowns
%
% Preprocessing: records < 100 ms are assigned to the previous trial or
% treated as anticipatory responses (trial excluded); see Settings.
% Responses without a measurable RT count for accuracy indices only.
%
% Output folder: <OUT_DIR>/results_sart
% Requirements: MATLAB R2016b or newer. No toolboxes needed.
% -------------------------------------------------------------------------
clear; clc;

%% ----------------------------- Settings ---------------------------------
config;
cfg.dataDir      = fullfile(DATA_DIR, 'sart_behavioral');
cfg.outDir       = fullfile(OUT_DIR, 'results_sart');
cfg.recursive    = true;     % also search sub-folders of dataDir
% N-back folder: the actual SART/N-back order in each block is taken from
% the file time stamps (N-back starting after SART = SN, before = NS).
% If a block's N-back file is missing, the planned order from the design
% table is used and a note is written to the log.
cfg.nbackDir     = fullfile(DATA_DIR, 'n_back_behavioral');
cfg.targetDigit  = 3;        % No-Go digit
cfg.SOA          = 1.150;    % digit 250 ms + mask 900 ms (seconds)
cfg.anticipThr   = 0.100;    % records faster than this are not a response to the current digit (s)
% ---- Records faster than anticipThr ----
% Rules for a record < anticipThr in trial k (applied in trial order):
%   - trial k-1 had no response                 -> response of k-1, RT = SOA + record
%   - trial k-1's record was itself reassigned  -> response of k-1
%                                                  (RT = SOA + record if record >= heldThr,
%                                                   otherwise no RT)
%   - trial k-1 responded with RT >= lateRT and the record is < heldThr
%                                               -> continuation of k-1's response, not counted
%   - otherwise                                 -> anticipatory response, trial k excluded
% Trial k then takes its response from trial k+1 if that record is
% reassigned; otherwise trial k is handled according to openPolicy.
cfg.lateRT     = 0.700;   % previous-trial RT threshold for the continuation rule (s)
cfg.heldThr    = 0.010;   % record threshold for the continuation rule (s)
cfg.openPolicy = 'exclude';   % 'exclude' or 'noresponse' for trials without an own response record
% If responses without a measurable RT exceed this share of a block's
% responses, the RT indices of that block are set to NaN (accuracy is kept).
cfg.maxUnknownRTProp = 0.10;
cfg.minCommission = 3;       % min commission errors for commission-RT indices
cfg.minGoForExG   = 50;      % min correct-Go RTs for Ex-Gaussian fit
cfg.preErrorWin   = 4;       % trials before a No-Go used for pre-error RT
cfg.nBins         = 4;       % time-on-task bins (4 = quartiles)
cfg.makeFigures   = true;
% Block-quality rules (blocks failing them get BlockValid = 0)
cfg.maxAnticipProp      = 0.20;  % > 20% anticipatory presses (not driven by the digit) -> invalid
cfg.maxExcludedProp     = 0.30;  % > 30% unusable trials in total (anticipatory + unobserved) -> invalid
cfg.requirePositiveDp   = true;  % d' <= 0 (no better than chance) -> invalid
cfg.excludeInvalidBlocks = true; % drop invalid blocks from change scores & summaries
% Blocks retained regardless of the rules above, as [Subject Block] rows.
% S3 blocks 2-3: retained after removal of anticipatory trials.
cfg.forceValidBlocks = [3 2; 3 3];
% Ex-Gaussian fits with sigma or tau below this (s) are degenerate -> NaN
cfg.exgMinParam   = 0.005;

%% ------------------------- Experimental design ---------------------------
% Subject, Gender, T1, Int1, T2, Int2, T3
design = {
 'S1','Female','NS','N','NS','R','SN'
 'S2','Female','NS','N','SN','R','NS'
 'S3','Male','NS','N','SN','R','SN'
 'S4','Female','SN','E','NS','R','SN'
 'S5','Female','SN','E','SN','N','NS'
 'S6','Female','NS','E','SN','N','NS'
 'S7','Female','NS','E','SN','R','NS'
 'S8','Female','SN','R','NS','N','SN'
 'S9','Male','NS','E','NS','R','SN'
 'S10','Male','SN','N','NS','R','NS'
 'S11','Female','SN','N','SN','E','NS'
 'S12','Male','SN','R','NS','N','NS'
 'S13','Female','NS','R','SN','N','NS'
 'S14','Female','SN','N','NS','E','SN'
 'S15','Female','SN','E','SN','R','NS'
 'S16','Male','SN','R','NS','E','NS'
 'S17','Male','SN','E','NS','R','NS'
 'S18','Female','NS','R','SN','E','NS'
 'S19','Male','SN','R','SN','N','NS'
 'S20','Male','NS','R','SN','E','SN'
 'S21','Female','NS','N','SN','E','NS'
 'S22','Male','NS','R','NS','E','SN'
 'S23','Male','SN','N','SN','R','NS'
 'S24','Female','NS','R','NS','N','SN'
 'S25','Male','NS','E','NS','N','SN'
 'S26','Male','NS','N','NS','E','SN'
 'S27','Male','SN','N','NS','E','NS'
 'S28','Male','SN','E','NS','N','NS'
 'S29','Female','SN','R','SN','E','NS'
 'S30','Male','NS','E','SN','N','SN'
};
designIDs = cellfun(@(s) str2double(s(2:end)), design(:,1));

%% ----------------------------- Find files --------------------------------
if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end
if cfg.recursive
    files = dir(fullfile(cfg.dataDir, '**', '*.mat'));
else
    files = dir(fullfile(cfg.dataDir, '*.mat'));
end
files = files(~contains({files.folder}, cfg.outDir));   % skip results folder

logLines = {};
logLines{end+1} = sprintf('SART analysis run: %s', datestr(now));
logLines{end+1} = sprintf('Data folder: %s', cfg.dataDir);

fileInfo = struct('path',{},'name',{},'subj',{},'block',{},'stamp',{});
for k = 1:numel(files)
    nm  = files(k).name;
    tok = regexpi(nm, '^S(\d+)_SART.*?block(\d+)', 'tokens', 'once');
    if isempty(tok)
        logLines{end+1} = sprintf('SKIPPED (name not recognised): %s', nm); %#ok<SAGROW>
        continue
    end
    fi.path  = fullfile(files(k).folder, nm);
    fi.name  = nm;
    fi.subj  = str2double(tok{1});
    fi.block = str2double(tok{2});
    fi.stamp = parseStamp(nm);
    fileInfo(end+1) = fi; %#ok<SAGROW>
end
fprintf('%d SART files found.\n', numel(fileInfo));
if isempty(fileInfo), error('No SART files found in %s', cfg.dataDir); end

%% ------------------- N-back time stamps (actual order) -------------------
nbStamp = containers.Map('KeyType', 'char', 'ValueType', 'double');
if exist(cfg.nbackDir, 'dir')
    nbFiles = dir(fullfile(cfg.nbackDir, '*.mat'));
    for k = 1:numel(nbFiles)
        tok = regexpi(nbFiles(k).name, '^s(\d+)_Nback.*?pre(\d+)', 'tokens', 'once');
        if isempty(tok), continue; end
        st = parseStamp(nbFiles(k).name);
        if isnat(st), continue; end
        nbStamp(sprintf('%d_%d', str2double(tok{1}), str2double(tok{2}))) = datenum(st);
    end
    fprintf('%d N-back files found (for task order).\n', nbStamp.Count);
else
    logLines{end+1} = sprintf('N-back folder not found (%s): planned task order used.', cfg.nbackDir);
end
orderRows = {};

%% ------------------------------ Process ----------------------------------
blockRows = {};  todRows = {};  fontRows = {};  trialTabs = {};  checkRows = {};

for k = 1:numel(fileInfo)
    fi = fileInfo(k);
    fprintf('Processing %s\n', fi.name);
    try
        S = load(fi.path);
        if ~isfield(S, 'data'), error('variable "data" not found'); end
        D = S.data;
        digit    = double(D.digit(:));
        fontSize = double(D.fontSize(:));
        rt       = double(D.RT(:));
        resp     = double(D.response(:));
    catch ME
        logLines{end+1} = sprintf('ERROR reading %s: %s', fi.name, ME.message); %#ok<SAGROW>
        continue
    end

    % ---- design info ----
    row = find(designIDs == fi.subj, 1);
    meta = struct();
    meta.Subject = fi.subj;
    meta.Block   = fi.block;
    if ~isempty(row) && fi.block >= 1 && fi.block <= 3
        taskCols = [3 5 7];
        meta.Gender    = design{row, 2};
        meta.TaskOrder = design{row, taskCols(fi.block)};
        if strcmp(meta.TaskOrder, 'SN'), meta.SARTPosition = 'First';
        else,                            meta.SARTPosition = 'AfterNback'; end
        switch fi.block
            case 1, meta.Intervention = 'Baseline'; meta.InterventionPosition = 0;
            case 2, meta.Intervention = design{row, 4}; meta.InterventionPosition = 1;
            case 3, meta.Intervention = design{row, 6}; meta.InterventionPosition = 2;
        end
    else
        meta.Gender = 'NA'; meta.TaskOrder = 'NA'; meta.SARTPosition = 'NA';
        meta.Intervention = 'NA'; meta.InterventionPosition = NaN;
        logLines{end+1} = sprintf('WARNING: S%d block %d not in design table (%s)', ...
            fi.subj, fi.block, fi.name); %#ok<SAGROW>
    end
    % ---- actual task order from time stamps ----
    meta.TaskOrderPlanned   = meta.TaskOrder;
    meta.NbackMinusSART_min = NaN;
    key = sprintf('%d_%d', fi.subj, fi.block);
    if isKey(nbStamp, key) && ~isnat(fi.stamp)
        dmin = (nbStamp(key) - datenum(fi.stamp)) * 24 * 60;
        meta.NbackMinusSART_min = dmin;
        if dmin > 0, meta.TaskOrder = 'SN'; else, meta.TaskOrder = 'NS'; end
        if abs(dmin) > 20
            logLines{end+1} = sprintf('WARNING: unusual SART/N-back gap for S%d block %d (%.1f min)', ...
                fi.subj, fi.block, dmin); %#ok<SAGROW>
        end
        if ~strcmp(meta.TaskOrder, meta.TaskOrderPlanned)
            logLines{end+1} = sprintf('ORDER CORRECTED: S%d block %d planned %s, actual %s (N-back started %.1f min %s SART)', ...
                fi.subj, fi.block, meta.TaskOrderPlanned, meta.TaskOrder, abs(dmin), ...
                ternary(dmin > 0, 'after', 'before')); %#ok<SAGROW>
        end
    else
        logLines{end+1} = sprintf('ORDER: no N-back file for S%d block %d, planned order used', ...
            fi.subj, fi.block); %#ok<SAGROW>
    end
    switch meta.TaskOrder
        case 'SN', meta.SARTPosition = 'First';
        case 'NS', meta.SARTPosition = 'AfterNback';
        otherwise, meta.SARTPosition = 'NA';
    end
    orderRows{end+1} = struct('Subject', fi.subj, 'Block', fi.block, ...
        'Planned', meta.TaskOrderPlanned, 'Actual', meta.TaskOrder, ...
        'NbackMinusSART_min', meta.NbackMinusSART_min, ...
        'Mismatch', double(~strcmp(meta.TaskOrder, meta.TaskOrderPlanned))); %#ok<SAGROW>

    if isnat(fi.stamp), meta.FileTime = '';
    else, meta.FileTime = datestr(datenum(fi.stamp), 'yyyy-mm-dd HH:MM:SS'); end
    meta.File = fi.name;

    % ---- preprocessing + metrics ----
    [tr, M, tod, fz] = sartCompute(digit, fontSize, rt, resp, cfg);

    blockRows{end+1} = catStruct(meta, M); %#ok<SAGROW>

    for b = 1:numel(tod)
        todRows{end+1} = catStruct(pickMeta(meta), tod(b)); %#ok<SAGROW>
    end
    for b = 1:numel(fz)
        fontRows{end+1} = catStruct(pickMeta(meta), fz(b)); %#ok<SAGROW>
    end

    n  = numel(digit);
    TT = struct2table(tr);
    TT = [table(repmat(fi.subj,n,1), repmat(fi.block,n,1), ...
                repmat({meta.Intervention},n,1), repmat({meta.TaskOrder},n,1), ...
                'VariableNames', {'Subject','Block','Intervention','TaskOrder'}), TT]; %#ok<AGROW>
    trialTabs{end+1} = TT; %#ok<SAGROW>

    % ---- data checks ----
    dc = struct();
    dc.Subject = fi.subj; dc.Block = fi.block; dc.File = fi.name; dc.FileTime = meta.FileTime;
    dc.nTrials = n;
    cnt = arrayfun(@(d) sum(digit == d), 1:9);
    dc.DigitCountsOK = double(n == 225 && all(cnt == 25));
    dc.nRespRTMismatch = M.nRespRTMismatch;
    dc.nAnticipatory   = M.nAnticipatory;
    dc.nMovedToPrev    = M.nMovedToPrev;
    dc.nHeldDuplicate  = M.nHeldDuplicate;
    dc.nAnticipExcl    = M.nAnticipExcl;
    dc.nOwnRespUnknown = M.nOwnRespUnknown;
    dc.nUnknownRT      = M.nUnknownRT;
    dc.nExcluded       = M.nExcluded;
    checkRows{end+1} = dc; %#ok<SAGROW>
    if M.nHeldDuplicate + M.nMovedToPrev >= 10
        logLines{end+1} = sprintf(['SPILL-OVER: S%d block %d: %d records moved to the previous trial, ' ...
            '%d held-key duplicates, %d anticipatory, %d trials with unobserved response, %d responses with unknown RT%s'], ...
            fi.subj, fi.block, M.nMovedToPrev, M.nHeldDuplicate, M.nAnticipExcl, M.nOwnRespUnknown, M.nUnknownRT, ...
            ternary(M.RTValid == 1, '', ' -> RT indices set to NaN')); %#ok<SAGROW>
    end
    if ~dc.DigitCountsOK
        logLines{end+1} = sprintf('WARNING: %s has %d trials / unexpected digit counts [%s]', ...
            fi.name, n, num2str(cnt)); %#ok<SAGROW>
    end
end

blockT = struct2table([blockRows{:}]);
todT   = struct2table([todRows{:}]);
fontT  = struct2table([fontRows{:}]);
trialT = vertcat(trialTabs{:});
checkT = struct2table([checkRows{:}]);
orderT = sortrows(struct2table([orderRows{:}]), {'Subject','Block'});
blockT = sortrows(blockT, {'Subject','Block'});
checkT = sortrows(checkT, {'Subject','Block'});

%% ---------------------- Duplicates / missing files -----------------------
[~, ~, g] = unique(blockT.Subject * 100 + blockT.Block);
dupCount = accumarray(g, 1);
blockT.IsDuplicate = double(dupCount(g) > 1);
blockT.UseForChange = ones(height(blockT), 1);
for gi = find(dupCount > 1)'
    idx = find(g == gi);
    logLines{end+1} = sprintf('WARNING: duplicate files for S%d block %d: %s', ...
        blockT.Subject(idx(1)), blockT.Block(idx(1)), strjoin(blockT.File(idx)', ' | ')); %#ok<SAGROW>
    % keep the latest file (by file time) for change scores
    [~, ord] = sort(blockT.FileTime(idx));
    blockT.UseForChange(idx) = 0;
    blockT.UseForChange(idx(ord(end))) = 1;
end
% ---- block quality ----
blockT.PropExcluded = blockT.nExcluded ./ (blockT.nValid + blockT.nExcluded);
blockT.PropAnticip  = blockT.nAnticipExcl ./ (blockT.nValid + blockT.nExcluded);
blockT.BlockValid = double(blockT.PropAnticip <= cfg.maxAnticipProp & ...
                           blockT.PropExcluded <= cfg.maxExcludedProp);
if cfg.requirePositiveDp
    blockT.BlockValid = blockT.BlockValid .* double(blockT.dprime > 0);
end
for r = 1:size(cfg.forceValidBlocks, 1)
    i = find(blockT.Subject == cfg.forceValidBlocks(r,1) & blockT.Block == cfg.forceValidBlocks(r,2));
    if ~isempty(i) && blockT.BlockValid(i(1)) == 0
        blockT.BlockValid(i) = 1;
        logLines{end+1} = sprintf('FORCED VALID: S%d block %d (anticipatory = %.0f%%, d'' = %.2f)', ...
            blockT.Subject(i(1)), blockT.Block(i(1)), 100 * blockT.PropAnticip(i(1)), blockT.dprime(i(1))); %#ok<SAGROW>
    end
end
for i = find(blockT.BlockValid == 0)'
    logLines{end+1} = sprintf('INVALID BLOCK: S%d block %d (anticipatory = %.0f%%, all excluded = %.0f%%, d'' = %.2f)', ...
        blockT.Subject(i), blockT.Block(i), 100 * blockT.PropAnticip(i), 100 * blockT.PropExcluded(i), blockT.dprime(i)); %#ok<SAGROW>
end
for i = find(blockT.ExG_valid == 0)'
    logLines{end+1} = sprintf('ExG not usable: S%d block %d', blockT.Subject(i), blockT.Block(i)); %#ok<SAGROW>
end
if cfg.excludeInvalidBlocks
    useRows = blockT.UseForChange == 1 & blockT.BlockValid == 1;
else
    useRows = blockT.UseForChange == 1;
end

for s = designIDs'
    for b = 1:3
        if ~any(blockT.Subject == s & blockT.Block == b)
            logLines{end+1} = sprintf('MISSING: S%d block %d', s, b); %#ok<SAGROW>
        end
    end
end

%% ------------------- Change scores (post - pre block) --------------------
metricList = {'CommissionRate','OmissionRate','dprime','c', ...
              'GoRT_mean','GoRT_median','GoRT_SD','GoRT_CV', ...
              'CommRT_mean','CommRT_median','GoMinusCommRT', ...
              'ExG_mu','ExG_sigma','ExG_tau', ...
              'PreErrorSpeeding','PES_traditional','PES_robust'};
B = blockT(useRows, :);
chRows = {};
for i = 1:height(B)
    if B.Block(i) < 2, continue; end
    j = find(B.Subject == B.Subject(i) & B.Block == B.Block(i) - 1, 1);
    if isempty(j), continue; end
    r = struct();
    r.Subject = B.Subject(i);
    r.Gender  = B.Gender{i};
    r.Intervention = B.Intervention{i};
    r.InterventionPosition = B.InterventionPosition(i);
    r.PreBlock = B.Block(j); r.PostBlock = B.Block(i);
    r.TaskOrderPre  = B.TaskOrder{j};
    r.TaskOrderPost = B.TaskOrder{i};
    for m = 1:numel(metricList)
        mn = metricList{m};
        r.(['pre_'  mn]) = B.(mn)(j);
        r.(['post_' mn]) = B.(mn)(i);
        r.(['delta_' mn]) = B.(mn)(i) - B.(mn)(j);
    end
    chRows{end+1} = r; %#ok<SAGROW>
end
if ~isempty(chRows)
    changeT = sortrows(struct2table([chRows{:}]), {'Subject','PostBlock'});
else
    changeT = table();
end

%% -------------------------- Descriptive summaries ------------------------
condLevels = {'Baseline','E','N','R'};
sumBlock = describeBy(blockT(useRows, :), 'Intervention', condLevels, metricList, '');
if ~isempty(changeT)
    sumChange = describeBy(changeT, 'Intervention', {'E','N','R'}, metricList, 'delta_');
else
    sumChange = table();
end

%% -------------------------------- Save -----------------------------------
writetable(blockT,  fullfile(cfg.outDir, 'SART_block_metrics.csv'));
writetable(changeT, fullfile(cfg.outDir, 'SART_change_scores.csv'));
writetable(todT,    fullfile(cfg.outDir, 'SART_time_on_task.csv'));
writetable(fontT,   fullfile(cfg.outDir, 'SART_font_size.csv'));
writetable(trialT,  fullfile(cfg.outDir, 'SART_trial_level_clean.csv'));
writetable(checkT,  fullfile(cfg.outDir, 'SART_data_check.csv'));
writetable(orderT,  fullfile(cfg.outDir, 'SART_task_order_check.csv'));
writetable(sumBlock,  fullfile(cfg.outDir, 'SART_summary_by_condition.csv'));
writetable(sumChange, fullfile(cfg.outDir, 'SART_summary_change_by_intervention.csv'));

xlsxFile = fullfile(cfg.outDir, 'SART_results.xlsx');
if exist(xlsxFile, 'file'), delete(xlsxFile); end
try
    writetable(blockT,    xlsxFile, 'Sheet', 'BlockMetrics');
    writetable(changeT,   xlsxFile, 'Sheet', 'ChangeScores');
    writetable(sumBlock,  xlsxFile, 'Sheet', 'SummaryByCondition');
    writetable(sumChange, xlsxFile, 'Sheet', 'SummaryChange');
    writetable(todT,      xlsxFile, 'Sheet', 'TimeOnTask');
    writetable(fontT,     xlsxFile, 'Sheet', 'FontSize');
    writetable(checkT,    xlsxFile, 'Sheet', 'DataCheck');
    writetable(orderT,    xlsxFile, 'Sheet', 'TaskOrder');
catch ME
    logLines{end+1} = sprintf('Could not write xlsx: %s', ME.message);
end

save(fullfile(cfg.outDir, 'SART_results.mat'), 'blockT', 'changeT', 'todT', ...
     'fontT', 'trialT', 'checkT', 'orderT', 'sumBlock', 'sumChange', 'cfg');

fid = fopen(fullfile(cfg.outDir, 'SART_log.txt'), 'w');
fprintf(fid, '%s\n', logLines{:});
fclose(fid);
fprintf('%s\n', logLines{:});

%% ------------------------------- Figures ---------------------------------
if cfg.makeFigures && ~isempty(changeT)
    figDir = fullfile(cfg.outDir, 'figures');
    if ~exist(figDir, 'dir'), mkdir(figDir); end
    plotList = {'CommissionRate','OmissionRate','dprime','c','GoRT_mean', ...
                'GoRT_CV','CommRT_median','ExG_tau'};
    ints = {'E','N','R'};  intNames = {'Reels (E)','Neutral (N)','Rest (R)'};
    cols = [0.85 0.33 0.10; 0.30 0.45 0.69; 0.47 0.67 0.19];
    fig = figure('Color','w','Position',[100 100 1300 650],'Visible','off');
    for p = 1:numel(plotList)
        subplot(2, 4, p); hold on
        v = changeT.(['delta_' plotList{p}]);
        for q = 1:3
            y = v(strcmp(changeT.Intervention, ints{q}));
            y = y(~isnan(y));
            if isempty(y), continue; end
            bar(q, mean(y), 0.6, 'FaceColor', cols(q,:), 'FaceAlpha', 0.35, 'EdgeColor', 'none');
            errorbar(q, mean(y), std(y)/sqrt(numel(y)), 'k', 'LineWidth', 1.2, 'CapSize', 8);
            scatter(q + (rand(size(y)) - 0.5) * 0.3, y, 18, cols(q,:), 'filled');
        end
        plot([0.4 3.6], [0 0], 'k:');
        set(gca, 'XTick', 1:3, 'XTickLabel', intNames, 'Box', 'off');
        title(['\Delta ' strrep(plotList{p}, '_', ' ')]);
        xlim([0.4 3.6]);
    end
    if exist('sgtitle', 'file') || exist('sgtitle', 'builtin')
        sgtitle('Change after intervention (post block - pre block), mean \pm SE');
    end
    print(fig, fullfile(figDir, 'SART_change_by_intervention.png'), '-dpng', '-r150');
    close(fig);
end

fprintf('\nDone. Results saved in:\n%s\n', cfg.outDir);


%% =========================================================================
%%                              Local functions
%% =========================================================================
function [tr, M, tod, fz] = sartCompute(digit, fontSize, rt, resp, cfg)
% Preprocess one SART block and compute all indices.
n = numel(digit);
digit = digit(:); fontSize = fontSize(:); rt = rt(:); resp = resp(:);
resp(isnan(resp)) = 0;
rawRT = rt; rawResp = resp;

M.nRespRTMismatch = sum((resp == 1 & isnan(rt)) | (resp == 0 & ~isnan(rt)));
resp(resp == 1 & isnan(rt)) = 0;      % a "response" without RT is unusable
rt(resp == 0) = NaN;

r0 = resp; t0 = rt;      % records after the consistency check
valid    = true(n, 1);
% action: 1 = record assigned to previous trial, 2 = anticipatory (trial excluded),
%         3 = continuation of the previous trial's response (record not counted)
action   = zeros(n, 1);
gotLate  = false(n, 1);   % trial received its response from the next trial's record
openT    = false(n, 1);   % own response not available (record reassigned)
isFast   = r0 == 1 & t0 < cfg.anticipThr;
M.nAnticipatory = sum(isFast);     % all records < anticipThr (before classification)
for k = find(isFast)'
    rk = t0(k);
    resp(k) = 0; rt(k) = NaN;
    if k == 1
        valid(k) = false; action(k) = 2; continue
    end
    p = k - 1;
    if r0(p) == 0                                % no response in k-1: late response to k-1
        resp(p) = 1; rt(p) = cfg.SOA + rk; gotLate(p) = true;
        action(k) = 1; openT(k) = true;
    elseif isFast(p) && openT(p)                 % k-1's record was itself reassigned
        resp(p) = 1; gotLate(p) = true; openT(p) = false;
        if rk >= cfg.heldThr, rt(p) = cfg.SOA + rk; else, rt(p) = NaN; end
        action(k) = 1; openT(k) = true;
    elseif ~isFast(p) && t0(p) >= cfg.lateRT && rk < cfg.heldThr
        action(k) = 3; openT(k) = true;          % continuation of k-1's response
    else                                         % anticipatory press for trial k
        valid(k) = false; action(k) = 2;
    end
end
% trials whose own response was never observed
nOpen = sum(openT);
if strcmpi(cfg.openPolicy, 'exclude')
    valid(openT) = false;
end
rtUnknown = resp == 1 & isnan(rt);           % response present, RT not measurable
M.nMovedToPrev   = sum(action == 1);
M.nHeldDuplicate = sum(action == 3);
M.nAnticipExcl   = sum(action == 2);
M.nOwnRespUnknown = nOpen;
M.nExcluded    = sum(~valid);
M.nValid       = sum(valid);

isT = digit == cfg.targetDigit;
cGo  = valid & ~isT & resp == 1;     % correct Go
om   = valid & ~isT & resp == 0;     % omission
com  = valid &  isT & resp == 1;     % commission
cw   = valid &  isT & resp == 0;     % correct withhold
cGoK = cGo & ~isnan(rt);             % correct Go with measurable RT
comK = com & ~isnan(rt);
% share of responses without a measurable RT
nRespValid = sum(valid & resp == 1);
M.nUnknownRT   = sum(valid & rtUnknown);
M.PropUnknownRT = M.nUnknownRT / max(nRespValid, 1);
M.RTValid = double(M.PropUnknownRT <= cfg.maxUnknownRTProp);
if ~M.RTValid, cGoK(:) = false; comK(:) = false; end

% ---- accuracy & SDT ----
nGo = sum(valid & ~isT);  nNoGo = sum(valid & isT);
M.nGo = nGo;  M.nNoGo = nNoGo;
M.nCommission = sum(com);  M.nOmission = sum(om);
M.CommissionRate = sum(com) / nNoGo;
M.OmissionRate   = sum(om)  / nGo;
M.HitRate = sum(cGo) / nGo;
M.FARate  = M.CommissionRate;
Hc = (sum(cGo) + 0.5) / (nGo + 1);       % log-linear correction
Fc = (sum(com) + 0.5) / (nNoGo + 1);
M.dprime = z(Hc) - z(Fc);
M.c      = -(z(Hc) + z(Fc)) / 2;

% ---- correct-Go RT ----
g = rt(cGoK);
M.GoRT_n      = numel(g);
M.GoRT_mean   = mean(g);
if isempty(g), M.GoRT_median = NaN; else, M.GoRT_median = median(g); end
M.GoRT_SD     = std(g);
M.GoRT_CV     = M.GoRT_SD / M.GoRT_mean;

% ---- commission RT ----
cr = rt(comK);
M.CommRT_n = numel(cr);
if numel(cr) >= cfg.minCommission
    M.CommRT_mean   = mean(cr);
    M.CommRT_median = median(cr);
    M.GoMinusCommRT = M.GoRT_median - M.CommRT_median;
else
    M.CommRT_mean = NaN; M.CommRT_median = NaN; M.GoMinusCommRT = NaN;
end

% ---- Ex-Gaussian ----
if numel(g) >= cfg.minGoForExG
    [M.ExG_mu, M.ExG_sigma, M.ExG_tau, M.ExG_converged] = fitExGauss(g);
    % sigma or tau collapsing to ~0 means the RT distribution is not
    % ex-Gaussian shaped (e.g. symmetric/left-skewed): the fit is not usable
    M.ExG_valid = double(M.ExG_sigma >= cfg.exgMinParam && M.ExG_tau >= cfg.exgMinParam);
    if ~M.ExG_valid
        M.ExG_mu = NaN; M.ExG_sigma = NaN; M.ExG_tau = NaN;
    end
else
    M.ExG_mu = NaN; M.ExG_sigma = NaN; M.ExG_tau = NaN; M.ExG_converged = 0;
    M.ExG_valid = 0;
end

% ---- pre-error speeding ----
M.PreErrorRT_commission = preRT(find(com), rt, cGoK, cfg.preErrorWin); %#ok<FNDSB>
M.PreErrorRT_withhold   = preRT(find(cw),  rt, cGoK, cfg.preErrorWin); %#ok<FNDSB>
M.PreErrorSpeeding      = M.PreErrorRT_withhold - M.PreErrorRT_commission;

% ---- post-error slowing ----
iC = find(com); iC = iC(iC < n);  iW = find(cw); iW = iW(iW < n);
postC = rt(iC(cGoK(iC + 1)) + 1);
postW = rt(iW(cGoK(iW + 1)) + 1);
M.PostCommRT     = mean(postC);
M.PostWithholdRT = mean(postW);
M.PES_traditional = M.PostCommRT - M.PostWithholdRT;
iR = find(com); iR = iR(iR > 1 & iR < n);
iR = iR(cGoK(iR - 1) & cGoK(iR + 1));
if isempty(iR), M.PES_robust = NaN;
else, M.PES_robust = mean(rt(iR + 1) - rt(iR - 1)); end

% ---- time on task ----
bin = ceil((1:n)' / n * cfg.nBins);
tod = struct('Bin',{},'nGo',{},'nNoGo',{},'CommissionRate',{},'OmissionRate',{}, ...
             'GoRT_mean',{},'GoRT_median',{},'GoRT_CV',{});
for b = 1:cfg.nBins
    s = bin == b;
    t.Bin = b;
    t.nGo   = sum(s & valid & ~isT);
    t.nNoGo = sum(s & valid & isT);
    if t.nNoGo > 0, t.CommissionRate = sum(s & com) / t.nNoGo; else, t.CommissionRate = NaN; end
    t.OmissionRate   = sum(s & om)  / max(t.nGo, 1);
    gb = rt(s & cGoK);
    if isempty(gb), gb = NaN; end
    t.GoRT_mean = mean(gb); t.GoRT_median = median(gb); t.GoRT_CV = std(gb) / mean(gb);
    tod(end+1) = t; %#ok<AGROW>
end

% ---- font size ----
fz = struct('FontSize',{},'nGo',{},'nNoGo',{},'CommissionRate',{},'GoRT_mean',{},'GoRT_median',{});
for f = unique(fontSize(:))'
    s = fontSize == f;
    u.FontSize = f;
    u.nGo   = sum(s & valid & ~isT);
    u.nNoGo = sum(s & valid & isT);
    if u.nNoGo > 0, u.CommissionRate = sum(s & com) / u.nNoGo; else, u.CommissionRate = NaN; end
    gf = rt(s & cGoK);
    if isempty(gf), gf = NaN; end
    u.GoRT_mean = mean(gf); u.GoRT_median = median(gf);
    fz(end+1) = u; %#ok<AGROW>
end

% ---- trial-level output ----
type = repmat({'Excluded'}, n, 1);
type(cGo) = {'CorrectGo'}; type(om) = {'Omission'};
type(com) = {'Commission'}; type(cw) = {'CorrectWithhold'};
tr.Trial    = (1:n)';
tr.Digit    = digit;
tr.FontSize = fontSize;
tr.RT_raw   = rawRT;
tr.Resp_raw = rawResp;
tr.RT       = rt;
tr.Resp     = resp;
tr.Valid    = double(valid);
tr.Action   = action;          % 1 assigned to previous, 2 anticipatory (excluded), 3 continuation of previous response
tr.OwnRespUnknown = double(openT);   % own response never observed
tr.RTUnknown = double(rtUnknown);    % response present but RT not measurable
tr.LateRespReceived = double(gotLate);
tr.TrialType = type;
tr.TimeBin  = bin;
end

function v = preRT(events, rt, cGo, W)
% Mean correct-Go RT in the W trials preceding each event, averaged over events.
vals = [];
for i = events(:)'
    j = max(1, i - W):(i - 1);
    j = j(cGo(j));
    if ~isempty(j), vals(end+1) = mean(rt(j)); end %#ok<AGROW>
end
if isempty(vals), v = NaN; else, v = mean(vals); end
end

function [mu, sigma, tau, ok] = fitExGauss(x)
% Maximum-likelihood Ex-Gaussian fit (no toolbox needed).
x = x(~isnan(x)); x = x(:);
m = mean(x); s = std(x);
sk = mean((x - m).^3) / s^3;
tau0   = min(s * (max(sk, 0.1) / 2)^(1/3), 0.9 * s);
sigma0 = sqrt(s^2 - tau0^2);
mu0    = m - tau0;
nll  = @(th) -sum(exgLogPdf(x, th(1), exp(th(2)), exp(th(3))));
opts = optimset('MaxFunEvals', 1e4, 'MaxIter', 1e4, 'TolX', 1e-8, 'TolFun', 1e-8, 'Display', 'off');
[th, ~, flag] = fminsearch(nll, [mu0, log(sigma0), log(tau0)], opts);
mu = th(1); sigma = exp(th(2)); tau = exp(th(3)); ok = double(flag == 1);
end

function lp = exgLogPdf(x, mu, sigma, tau)
zz = (x - mu) ./ sigma - sigma ./ tau;
t  = -zz ./ sqrt(2);
logPhi = zeros(size(t));
pos = t > 0;
logPhi(pos)  = log(0.5) + log(erfcx(t(pos))) - t(pos).^2;   % numerically stable tail
logPhi(~pos) = log(0.5 * erfc(t(~pos)));
lp = -log(tau) + (mu - x) ./ tau + sigma.^2 ./ (2 * tau.^2) + logPhi;
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function y = z(p)
% Inverse standard normal CDF (same as norminv, no toolbox needed).
y = -sqrt(2) * erfcinv(2 * p);
end

function dt = parseStamp(nm)
% Extracts e.g. "Sat-Feb-14-2026-14-33-47" from the file name.
dt = NaT;
tok = regexp(nm, '[A-Za-z]{3}-([A-Za-z]{3})-(\d{1,2})-(\d{4})-(\d{1,2})-(\d{1,2})-(\d{1,2})', 'tokens', 'once');
if isempty(tok), return; end
try
    dt = datetime(sprintf('%s %s %s %s:%s:%s', tok{:}), ...
        'InputFormat', 'MMM d yyyy H:m:s', 'Locale', 'en_US');
catch
    dt = NaT;
end
end

function s = catStruct(a, b)
s = a;
f = fieldnames(b);
for i = 1:numel(f), s.(f{i}) = b.(f{i}); end
end

function m = pickMeta(meta)
m.Subject = meta.Subject; m.Block = meta.Block;
m.Intervention = meta.Intervention; m.TaskOrder = meta.TaskOrder;
end

function T = describeBy(D, groupVar, levels, metrics, prefix)
% Mean, SD, SE and n of each metric per group level.
rows = {};
for L = 1:numel(levels)
    sel = strcmp(D.(groupVar), levels{L});
    for m = 1:numel(metrics)
        v = D.([prefix metrics{m}])(sel);
        v = v(~isnan(v));
        r.Condition = levels{L};
        r.Metric = metrics{m};
        r.n = numel(v);
        r.Mean = mean(v);
        r.SD = std(v);
        r.SE = std(v) / sqrt(max(numel(v), 1));
        rows{end+1} = r; %#ok<AGROW>
    end
end
T = struct2table([rows{:}]);
end
