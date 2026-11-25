%% N-back behavioral analysis
% -------------------------------------------------------------------------
% Task: 0 -> 1 -> 2 -> 3-back, repeated 3 times (12 sub-blocks x 20 trials)
% in each testing block (file s<N>_Nback-data-pre<B>-...mat).
%
% The same set of indices is computed at three levels of aggregation:
%   1) sub-block   : each of the 12 sub-blocks (load x repetition) -> time course
%   2) load level  : the 3 sub-blocks of each load combined (0-, 1-, 2-, 3-back)
%   3) overall     : all 12 sub-blocks combined
% Indices:
%   hits, misses, false alarms (FA), correct rejections; hit rate, FA rate,
%   accuracy, d', c;
%   hit RT (mean, median, SD, CV) and additional hit RT and FA indices;
%   FA RT (median).
% Block-level composites: working-memory composite (mean over 1-3-back) and
% load slopes (change per +1 n-back level).
%
% RT = responseTime - stimOnset (responses scored within the trial).
% Task order (SN / NS) comes from the SART and N-back file time stamps.
% Output folder: <OUT_DIR>/results_nback
% Requirements: MATLAB R2016b or newer. No toolboxes needed.
% -------------------------------------------------------------------------
clear; clc;

%% ----------------------------- Settings ---------------------------------
config;
cfg.dataDir   = fullfile(DATA_DIR, 'n_back_behavioral');
cfg.sartDir   = fullfile(DATA_DIR, 'sart_behavioral');
cfg.outDir    = fullfile(OUT_DIR, 'results_nback');
cfg.levels    = 0:3;
cfg.nRep      = 3;        % repetitions of the 0-1-2-3 cycle analysed per block
                          % (additional repetitions, if present, are not used)
cfg.stimDur   = 0.500;    % stimulus duration (s)
cfg.minRT     = 0.100;    % RTs below this are not used in RT indices (counted)
cfg.minRTn    = 3;        % min responses for RT mean/median
cfg.minRTnSD  = 5;        % min responses for RT SD / CV
cfg.min0backAcc = 0.80;   % 0-back accuracy below this -> block invalid
cfg.excludeInvalidBlocks = true;
cfg.makeFigures = true;

%% ------------------------- Experimental design ---------------------------
% Subject, Gender, Int1, Int2  (task order comes from the file time stamps)
design = {
 'S1','Female','N','R';  'S2','Female','N','R';  'S3','Male','N','R'
 'S4','Female','E','R';  'S5','Female','E','N';  'S6','Female','E','N'
 'S7','Female','E','R';  'S8','Female','R','N';  'S9','Male','E','R'
 'S10','Male','N','R';   'S11','Female','N','E'; 'S12','Male','R','N'
 'S13','Female','R','N'; 'S14','Female','N','E'; 'S15','Female','E','R'
 'S16','Male','R','E';   'S17','Male','E','R';   'S18','Female','R','E'
 'S19','Male','R','N';   'S20','Male','R','E';   'S21','Female','N','E'
 'S22','Male','R','E';   'S23','Male','N','R';   'S24','Female','R','N'
 'S25','Male','E','N';   'S26','Male','N','E';   'S27','Male','N','E'
 'S28','Male','E','N';   'S29','Female','R','E'; 'S30','Male','E','N'
};
designIDs = cellfun(@(s) str2double(s(2:end)), design(:,1));

%% ----------------------------- Find files --------------------------------
if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end
logLines = {sprintf('N-back analysis run: %s', datestr(now)), ...
            sprintf('Data folder: %s', cfg.dataDir)};

files = dir(fullfile(cfg.dataDir, '*.mat'));
fileInfo = struct('path',{},'name',{},'subj',{},'block',{},'stamp',{});
for k = 1:numel(files)
    nm  = files(k).name;
    tok = regexpi(nm, '^s(\d+)_Nback.*?pre(\d+)', 'tokens', 'once');
    if isempty(tok), continue; end
    fi.path = fullfile(files(k).folder, nm); fi.name = nm;
    fi.subj = str2double(tok{1}); fi.block = str2double(tok{2});
    fi.stamp = parseStamp(nm);
    fileInfo(end+1) = fi; %#ok<SAGROW>
end
fprintf('%d N-back files found.\n', numel(fileInfo));
if isempty(fileInfo), error('No N-back files found in %s', cfg.dataDir); end

sartStamp = containers.Map('KeyType', 'char', 'ValueType', 'double');
if exist(cfg.sartDir, 'dir')
    sf = dir(fullfile(cfg.sartDir, '**', '*.mat'));
    for k = 1:numel(sf)
        tok = regexpi(sf(k).name, '^S(\d+)_SART.*?block(\d+)', 'tokens', 'once');
        if isempty(tok), continue; end
        st = parseStamp(sf(k).name);
        if isnat(st), continue; end
        sartStamp(sprintf('%d_%d', str2double(tok{1}), str2double(tok{2}))) = datenum(st);
    end
else
    logLines{end+1} = 'SART folder not found: task order not available.';
end

%% ------------------------------ Process ----------------------------------
sbRows = {}; lvRows = {}; blRows = {}; trialTabs = {};

for k = 1:numel(fileInfo)
    fi = fileInfo(k);
    fprintf('Processing %s\n', fi.name);
    try
        S = load(fi.path);
        if ~isfield(S, 'data'), error('variable "data" not found'); end
        D = S.data;
        blk  = getf(D, 'block');    nb    = getf(D, 'nBack');
        trl  = getf(D, 'trial');    digit = getf(D, 'digit');
        isT  = getf(D, 'isTarget'); resp  = getf(D, 'responded');
        rt   = getf(D, 'responseTime') - getf(D, 'stimOnset');
        corrS = getf(D, 'correct');
    catch ME
        logLines{end+1} = sprintf('ERROR reading %s: %s', fi.name, ME.message); %#ok<SAGROW>
        continue
    end
    resp(isnan(resp)) = 0;

    % ---- meta ----
    meta = struct('Subject', fi.subj, 'Block', fi.block);
    row = find(designIDs == fi.subj, 1);
    if ~isempty(row) && fi.block >= 1 && fi.block <= 3
        meta.Gender = design{row, 2};
        ivs = {'Baseline', design{row, 3}, design{row, 4}};
        meta.Intervention = ivs{fi.block};
        meta.InterventionPosition = fi.block - 1;
    else
        meta.Gender = 'NA'; meta.Intervention = 'NA'; meta.InterventionPosition = NaN;
        logLines{end+1} = sprintf('WARNING: S%d block %d not in design table', fi.subj, fi.block); %#ok<SAGROW>
    end
    meta.TaskOrder = 'NA'; meta.NbackPosition = 'NA';
    key = sprintf('%d_%d', fi.subj, fi.block);
    dmin = NaN;
    if isKey(sartStamp, key) && ~isnat(fi.stamp)
        dmin = (datenum(fi.stamp) - sartStamp(key)) * 24 * 60;
        if dmin > 0, meta.TaskOrder = 'SN'; meta.NbackPosition = 'AfterSART';
        else,        meta.TaskOrder = 'NS'; meta.NbackPosition = 'First'; end
    else
        logLines{end+1} = sprintf('ORDER: no SART file for S%d block %d', fi.subj, fi.block); %#ok<SAGROW>
    end

    % ---- repetition index (0-1-2-3 cycle) and checks ----
    rep = zeros(size(blk));
    for L = cfg.levels
        bl = unique(blk(nb == L))';
        for r = 1:numel(bl), rep(blk == bl(r)) = r; end
    end
    [chk, expT] = checkFile(blk, nb, digit, isT, resp, rt, corrS);
    nRepFile = max(rep);
    if nRepFile ~= cfg.nRep
        logLines{end+1} = sprintf('NOTE %s: %d repetitions (%d trials); only repetitions 1-%d are analysed', ...
            fi.name, nRepFile, chk.nTrials, cfg.nRep); %#ok<SAGROW>
    end
    if chk.nTargetMismatch > 0 || chk.nCorrectMismatch > 0 || chk.nRespNoRT > 0
        logLines{end+1} = sprintf('CHECK %s: target-coding mismatches=%d, correct-field mismatches=%d, responses without RT=%d', ...
            fi.name, chk.nTargetMismatch, chk.nCorrectMismatch, chk.nRespNoRT); %#ok<SAGROW>
    end
    use = rep >= 1 & rep <= cfg.nRep;

    pm = pickMeta(meta);
    % ---- 1) sub-blocks ----
    sbs = unique(blk(use))';
    for b = sbs
        s = use & blk == b;
        M = nbMetrics(isT, resp, rt, s, cfg);
        sb = struct('SubBlock', find(sbs == b), 'Load', nb(find(s, 1)), 'Repetition', rep(find(s, 1)));
        sbRows{end+1} = catStruct(catStruct(pm, sb), M); %#ok<SAGROW>
    end
    % ---- 2) load levels ----
    Lv = struct([]);
    for L = cfg.levels
        M = nbMetrics(isT, resp, rt, use & nb == L, cfg);
        Lv = [Lv, M]; %#ok<AGROW>
        lvRows{end+1} = catStruct(catStruct(pm, struct('Load', L)), M); %#ok<SAGROW>
    end
    % ---- 3) overall + composites ----
    B = catStruct(meta, nbMetrics(isT, resp, rt, use, cfg));
    hi = ismember(cfg.levels, 1:3);
    B.dprime_WM       = mean([Lv(hi).dprime]);
    B.Accuracy_WM     = mean([Lv(hi).Accuracy]);
    B.HitRT_median_WM = mean([Lv(hi).HitRT_median], 'omitnan');
    B.ISIHitProp_WM   = mean([Lv(hi).ISIHitProp], 'omitnan');
    B.dprime_0back    = Lv(1).dprime;
    B.Accuracy_0back  = Lv(1).Accuracy;
    B.LoadSlope_dprime       = slope(cfg.levels, [Lv.dprime]);
    B.LoadSlope_Accuracy     = slope(cfg.levels, [Lv.Accuracy]);
    B.LoadSlope_HitRT_median = slope(cfg.levels, [Lv.HitRT_median]);
    B.LoadSlope_ISIHitProp   = slope(cfg.levels, [Lv.ISIHitProp]);
    B.NbackMinusSART_min = dmin;
    if isnat(fi.stamp), B.FileTime = '';
    else, B.FileTime = datestr(datenum(fi.stamp), 'yyyy-mm-dd HH:MM:SS'); end
    B.File = fi.name;
    B.nTrialsFile = chk.nTrials; B.nRepFile = nRepFile;
    B.nTargetMismatch = chk.nTargetMismatch; B.nCorrectMismatch = chk.nCorrectMismatch;
    blRows{end+1} = B; %#ok<SAGROW>

    % ---- trial level ----
    n = numel(blk);
    tr = struct('Subject', repmat(fi.subj, n, 1), 'Block', repmat(fi.block, n, 1), ...
        'SubBlock', blk(:), 'Load', nb(:), 'Repetition', rep(:), 'Used', double(use(:)), ...
        'Trial', trl(:), 'Digit', digit(:), 'isTarget', isT(:), 'isTarget_expected', expT(:), ...
        'Responded', resp(:), 'RT', rt(:), ...
        'Window', {repmat({''}, n, 1)});
    w = tr.Window; w(resp == 1 & rt <= cfg.stimDur) = {'Stim'}; w(resp == 1 & rt > cfg.stimDur) = {'ISI'};
    tr.Window = w;
    trialTabs{end+1} = struct2table(tr); %#ok<SAGROW>
end

sbT = sortrows(struct2table([sbRows{:}]), {'Subject','Block','SubBlock'});
lvT = sortrows(struct2table([lvRows{:}]), {'Subject','Block','Load'});
blT = sortrows(struct2table([blRows{:}]), {'Subject','Block'});
trialT = vertcat(trialTabs{:});

%% --------------------------- Validity / missing --------------------------
blT.BlockValid = double(blT.Accuracy_0back >= cfg.min0backAcc & blT.dprime > 0);
for i = find(blT.BlockValid == 0)'
    logLines{end+1} = sprintf('INVALID BLOCK: S%d block %d (0-back accuracy = %.2f, overall d'' = %.2f)', ...
        blT.Subject(i), blT.Block(i), blT.Accuracy_0back(i), blT.dprime(i)); %#ok<SAGROW>
end
[~, ~, g] = unique(blT.Subject * 100 + blT.Block);
dup = accumarray(g, 1);
for gi = find(dup > 1)'
    idx = find(g == gi);
    logLines{end+1} = sprintf('WARNING: duplicate files S%d block %d', blT.Subject(idx(1)), blT.Block(idx(1))); %#ok<SAGROW>
end
for s = designIDs'
    for b = 1:3
        if ~any(blT.Subject == s & blT.Block == b)
            logLines{end+1} = sprintf('MISSING: S%d block %d', s, b); %#ok<SAGROW>
        end
    end
end
sbT.BlockValid = lookupValid(sbT, blT);
lvT.BlockValid = lookupValid(lvT, blT);

%% ------------------------------ Change scores ----------------------------
mets = {'HitRate','FARate','Accuracy','dprime','c', ...
        'HitRT_mean','HitRT_median','HitRT_SD','HitRT_CV', ...
        'HitRT_stim_mean','HitRT_stim_median','HitRT_ISI_mean','HitRT_ISI_median', ...
        'ISIHitProp','FART_median'};
blkMets = [mets, {'dprime_WM','Accuracy_WM','HitRT_median_WM','ISIHitProp_WM', ...
        'LoadSlope_dprime','LoadSlope_Accuracy','LoadSlope_HitRT_median','LoadSlope_ISIHitProp'}];
if cfg.excludeInvalidBlocks
    LV = lvT(lvT.BlockValid == 1, :); BV = blT(blT.BlockValid == 1, :); SV = sbT(sbT.BlockValid == 1, :);
else
    LV = lvT; BV = blT; SV = sbT;
end
chLvT = changeScores(LV, mets, {'Load'});
chBlT = changeScores(BV, blkMets, {});
chSbT = changeScores(SV, mets, {'Load','Repetition'});

%% ------------------------------ Summaries --------------------------------
cond = {'Baseline','E','N','R'};
sumLvT = summarise(LV, {'Intervention','Load'}, cond, mets);
sumBlT = summarise(BV, {'Intervention'}, cond, blkMets);
sumSbT = summarise(SV, {'Intervention','Load','Repetition'}, cond, ...
    {'HitRate','FARate','Accuracy','dprime','HitRT_median','ISIHitProp'});

%% -------------------------------- Save -----------------------------------
writetable(sbT, fullfile(cfg.outDir, 'Nback_subblock_metrics.csv'));
writetable(lvT, fullfile(cfg.outDir, 'Nback_level_metrics.csv'));
writetable(blT, fullfile(cfg.outDir, 'Nback_block_metrics.csv'));
writetable(chSbT, fullfile(cfg.outDir, 'Nback_change_scores_subblock.csv'));
writetable(chLvT, fullfile(cfg.outDir, 'Nback_change_scores_level.csv'));
writetable(chBlT, fullfile(cfg.outDir, 'Nback_change_scores_block.csv'));
writetable(sumSbT, fullfile(cfg.outDir, 'Nback_summary_subblock.csv'));
writetable(sumLvT, fullfile(cfg.outDir, 'Nback_summary_level.csv'));
writetable(sumBlT, fullfile(cfg.outDir, 'Nback_summary_block.csv'));
writetable(trialT, fullfile(cfg.outDir, 'Nback_trial_level.csv'));
xf = fullfile(cfg.outDir, 'Nback_results.xlsx');
if exist(xf, 'file'), delete(xf); end
try
    writetable(blT, xf, 'Sheet', 'Overall');
    writetable(lvT, xf, 'Sheet', 'ByLoad');
    writetable(sbT, xf, 'Sheet', 'BySubBlock');
    writetable(chBlT, xf, 'Sheet', 'Change_Overall');
    writetable(chLvT, xf, 'Sheet', 'Change_ByLoad');
    writetable(sumBlT, xf, 'Sheet', 'Summary_Overall');
    writetable(sumLvT, xf, 'Sheet', 'Summary_ByLoad');
    writetable(sumSbT, xf, 'Sheet', 'Summary_BySubBlock');
catch ME
    logLines{end+1} = sprintf('Could not write xlsx: %s', ME.message);
end
save(fullfile(cfg.outDir, 'Nback_results.mat'), 'sbT', 'lvT', 'blT', 'chSbT', 'chLvT', ...
     'chBlT', 'sumSbT', 'sumLvT', 'sumBlT', 'trialT', 'cfg');
fid = fopen(fullfile(cfg.outDir, 'Nback_log.txt'), 'w');
fprintf(fid, '%s\n', logLines{:}); fclose(fid);
fprintf('%s\n', logLines{:});

%% ------------------------------- Figures ---------------------------------
if cfg.makeFigures
    figDir = fullfile(cfg.outDir, 'figures');
    if ~exist(figDir, 'dir'), mkdir(figDir); end
    cols  = [0.45 0.45 0.45; 0.85 0.33 0.10; 0.30 0.45 0.69; 0.47 0.67 0.19];
    names = {'Baseline','Reels (E)','Neutral (N)','Rest (R)'};
    % (a) by load
    pmet = {'dprime','Accuracy','HitRT_median','ISIHitProp'};
    plab = {'d''','Accuracy','Hit RT median (s)','Hits during ISI (prop.)'};
    fig = figure('Color', 'w', 'Position', [100 100 1200 800], 'Visible', 'off');
    for p = 1:4
        subplot(2, 2, p); hold on
        for c = 1:4
            [mu, se] = cellMeans(LV, pmet{p}, cond{c}, cfg.levels, []);
            errorbar(cfg.levels + (c - 2.5) * 0.06, mu, se, '-o', 'Color', cols(c,:), ...
                'MarkerFaceColor', cols(c,:), 'LineWidth', 1.4);
        end
        set(gca, 'XTick', cfg.levels, 'XTickLabel', {'0-back','1-back','2-back','3-back'}, 'Box', 'off');
        ylabel(plab{p}); xlim([-0.4 3.4]);
        if p == 1, legend(names, 'Location', 'southwest', 'Box', 'off'); end
    end
    print(fig, fullfile(figDir, 'Nback_by_load.png'), '-dpng', '-r150'); close(fig);
    % (b) time course: 12 sub-blocks
    fig = figure('Color', 'w', 'Position', [100 100 1300 700], 'Visible', 'off');
    pmet2 = {'dprime','HitRT_median'}; plab2 = {'d''','Hit RT median (s)'};
    for p = 1:2
        subplot(2, 1, p); hold on
        for c = 1:4
            mu = nan(1, 12); se = mu;
            for sbn = 1:12
                v = SV.(pmet2{p})(strcmp(SV.Intervention, cond{c}) & SV.SubBlock == sbn);
                v = v(~isnan(v)); mu(sbn) = mean(v); se(sbn) = std(v) / sqrt(max(numel(v), 1));
            end
            errorbar((1:12) + (c - 2.5) * 0.08, mu, se, '-o', 'Color', cols(c,:), ...
                'MarkerFaceColor', cols(c,:), 'LineWidth', 1.2);
        end
        lab = arrayfun(@(i) sprintf('%d-back r%d', mod(i-1,4), ceil(i/4)), 1:12, 'UniformOutput', false);
        set(gca, 'XTick', 1:12, 'XTickLabel', lab, 'Box', 'off');
        xtickangle(30); ylabel(plab2{p}); xlim([0.5 12.5]);
        for x = [4.5 8.5], plot([x x], ylim, 'k:'); end
        if p == 1, legend(names, 'Location', 'southwest', 'Box', 'off'); end
    end
    print(fig, fullfile(figDir, 'Nback_time_course.png'), '-dpng', '-r150'); close(fig);
end
fprintf('\nDone. Results saved in:\n%s\n', cfg.outDir);


%% =========================================================================
%%                              Local functions
%% =========================================================================
function v = getf(D, f)
v = nan(numel(D), 1);
for i = 1:numel(D)
    x = D(i).(f);
    if ~isempty(x), v(i) = double(x(1)); end
end
end

function [chk, expT] = checkFile(blk, nb, digit, isT, resp, rt, corrS)
chk.nTrials = numel(blk);
expT = zeros(size(isT));
for b = unique(blk(:))'
    idx = find(blk == b); n = nb(idx(1)); d = digit(idx);
    for i = 1:numel(idx)
        if n == 0, expT(idx(i)) = d(i) == 0;
        else,      expT(idx(i)) = i > n && d(i) == d(i - n); end
    end
end
chk.nTargetMismatch  = sum(expT ~= isT);
corrExp = double((isT == 1 & resp == 1) | (isT == 0 & resp == 0));
chk.nCorrectMismatch = sum(corrExp ~= corrS);
chk.nRespNoRT = sum(resp == 1 & isnan(rt));
end

function M = nbMetrics(isT, resp, rt, s, cfg)
T = s & isT == 1; N = s & isT == 0;
H = T & resp == 1; F = N & resp == 1;
M.nTarget = sum(T); M.nNonTarget = sum(N);
M.nHit = sum(H); M.nMiss = sum(T & resp == 0);
M.nFA  = sum(F); M.nCR   = sum(N & resp == 0);
M.HitRate  = M.nHit / M.nTarget;
M.FARate   = M.nFA / M.nNonTarget;
M.Accuracy = (M.nHit + M.nCR) / (M.nTarget + M.nNonTarget);
Hc = (M.nHit + 0.5) / (M.nTarget + 1);
Fc = (M.nFA + 0.5) / (M.nNonTarget + 1);
M.dprime = z(Hc) - z(Fc);
M.c      = -(z(Hc) + z(Fc)) / 2;
ok = rt >= cfg.minRT;
inStim = rt <= cfg.stimDur; inISI = rt > cfg.stimDur;
% hit RT
h = rt(H & ok);
M.HitRT_n = numel(h);
[M.HitRT_mean, M.HitRT_median] = mm(h, cfg.minRTn);
if numel(h) >= cfg.minRTnSD, M.HitRT_SD = std(h); M.HitRT_CV = std(h) / mean(h);
else, M.HitRT_SD = NaN; M.HitRT_CV = NaN; end
% additional hit RT indices
hs = rt(H & ok & inStim); hi = rt(H & inISI);
M.nHit_stim = numel(hs); M.nHit_ISI = numel(hi);
[M.HitRT_stim_mean, M.HitRT_stim_median] = mm(hs, cfg.minRTn);
[M.HitRT_ISI_mean,  M.HitRT_ISI_median]  = mm(hi, cfg.minRTn);
if M.nHit > 0, M.ISIHitProp = sum(H & inISI) / M.nHit; else, M.ISIHitProp = NaN; end
% false alarms
M.nFA_stim = sum(F & inStim); M.nFA_ISI = sum(F & inISI);
f = rt(F & ok);
[~, M.FART_median] = mm(f, cfg.minRTn);
M.nFastResp = sum(s & resp == 1 & ~ok);
end

function [a, b] = mm(x, nmin)
if numel(x) >= nmin, a = mean(x); b = median(x); else, a = NaN; b = NaN; end
end

function b = slope(x, y)
ok = ~isnan(y);
if sum(ok) < 3, b = NaN; return; end
p = polyfit(x(ok), y(ok), 1); b = p(1);
end

function v = lookupValid(T, blT)
v = zeros(height(T), 1);
for i = 1:height(T)
    j = find(blT.Subject == T.Subject(i) & blT.Block == T.Block(i), 1);
    v(i) = blT.BlockValid(j);
end
end

function C = changeScores(V, metrics, keys)
rows = {};
for i = 1:height(V)
    if V.Block(i) < 2, continue; end
    j = V.Subject == V.Subject(i) & V.Block == V.Block(i) - 1;
    for q = 1:numel(keys), j = j & V.(keys{q}) == V.(keys{q})(i); end
    j = find(j, 1);
    if isempty(j), continue; end
    r = struct('Subject', V.Subject(i), 'Gender', V.Gender{i}, ...
        'Intervention', V.Intervention{i}, 'InterventionPosition', V.InterventionPosition(i), ...
        'PreBlock', V.Block(j), 'PostBlock', V.Block(i), ...
        'TaskOrderPre', V.TaskOrder{j}, 'TaskOrderPost', V.TaskOrder{i});
    for q = 1:numel(keys), r.(keys{q}) = V.(keys{q})(i); end
    for m = 1:numel(metrics)
        mn = metrics{m};
        r.(['pre_' mn])   = V.(mn)(j);
        r.(['post_' mn])  = V.(mn)(i);
        r.(['delta_' mn]) = V.(mn)(i) - V.(mn)(j);
    end
    rows{end+1} = r; %#ok<AGROW>
end
if isempty(rows), C = table(); return; end
C = sortrows(struct2table([rows{:}]), [{'Subject','PostBlock'}, keys]);
end

function S = summarise(V, by, cond, metrics)
% mean / SD / SE / n of each metric for every combination of the grouping vars
rows = {};
if numel(by) == 1, nk = 1; keys = table();
else, keys = unique(V(:, by(2:end)), 'rows'); nk = height(keys); end
for c = 1:numel(cond)
    for kk = 1:nk
        s = strcmp(V.Intervention, cond{c});
        r0 = struct('Condition', cond{c});
        for q = 2:numel(by)
            s = s & V.(by{q}) == keys.(by{q})(kk);
            r0.(by{q}) = keys.(by{q})(kk);
        end
        for m = 1:numel(metrics)
            v = V.(metrics{m})(s); v = v(~isnan(v));
            r = r0; r.Metric = metrics{m}; r.n = numel(v);
            r.Mean = mean(v); r.SD = std(v); r.SE = std(v) / sqrt(max(numel(v), 1));
            rows{end+1} = r; %#ok<AGROW>
        end
    end
end
S = struct2table([rows{:}]);
end

function [mu, se] = cellMeans(V, met, cond, levels, ~)
mu = nan(size(levels)); se = mu;
for L = levels
    v = V.(met)(strcmp(V.Intervention, cond) & V.Load == L); v = v(~isnan(v));
    mu(L+1) = mean(v); se(L+1) = std(v) / sqrt(max(numel(v), 1));
end
end

function y = z(p)
y = -sqrt(2) * erfcinv(2 * p);
end

function dt = parseStamp(nm)
dt = NaT;
tok = regexp(nm, '[A-Za-z]{3}-([A-Za-z]{3})-(\d{1,2})-(\d{4})-(\d{1,2})-(\d{1,2})-(\d{1,2})', 'tokens', 'once');
if isempty(tok), return; end
try
    dt = datetime(sprintf('%s %s %s %s:%s:%s', tok{:}), 'InputFormat', 'MMM d yyyy H:m:s', 'Locale', 'en_US');
catch
    dt = NaT;
end
end

function s = catStruct(a, b)
s = a; f = fieldnames(b);
for i = 1:numel(f), s.(f{i}) = b.(f{i}); end
end

function m = pickMeta(meta)
m = struct('Subject', meta.Subject, 'Block', meta.Block, 'Gender', meta.Gender, ...
    'Intervention', meta.Intervention, 'InterventionPosition', meta.InterventionPosition, ...
    'TaskOrder', meta.TaskOrder, 'NbackPosition', meta.NbackPosition);
end
