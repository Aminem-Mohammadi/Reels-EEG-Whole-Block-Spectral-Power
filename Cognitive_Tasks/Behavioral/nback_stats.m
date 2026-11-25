%% N-back statistics: effect of the three interventions (E, N, R)
% -------------------------------------------------------------------------
% Input : results of nback_analysis_full.m (results_nback/*.csv)
% Output: results_nback/stats/
%
% All intervention analyses use change scores (post block - pre block)
% in linear mixed models with Satterthwaite df, controlling for
%   Position  : intervention 1st or 2nd           (centred)
%   OrderPre  : N-back after SART in pre block    (centred)
%   OrderPost : N-back after SART in post block   (centred)
% The models are fitted without intercept on cell dummies, so every
% coefficient is the adjusted mean change of one cell; all tests are
% contrasts on these cell means.
%
% A) Overall / composite indices (one change score per subject x post block)
%    delta ~ cells(E,N,R) + covariates + (1|Subject)
%    tests: E = N = R, pairwise (Holm), each intervention vs 0
% B) Per load level (0-3-back)
%    delta ~ cells(Intervention x Load) + covariates + (1|Subject) + (1|Subject:PostBlock)
%    tests: Intervention main effect, Load main effect, Intervention x Load,
%           simple effects of Intervention at each load (+ pairwise, Holm),
%           each cell vs 0
% C) Time course within the block (sub-blocks, repetitions 1-3; loads in cfg.tcLoads)
%    delta ~ cells(Intervention x Repetition) + Load + covariates
%            + (1|Subject) + (1|Subject:PostBlock)
%    tests: Intervention, Repetition, Intervention x Repetition,
%           simple effects of Intervention at each repetition
% D) Manipulation check: load effect on raw scores
%    y ~ Load + Block + (1|Subject)
% FDR (Benjamini-Hochberg) is applied to each family of main tests within
% the primary and secondary sets of indices. RT indices in ms.
% Requirements: MATLAB R2018b+ with Statistics and Machine Learning Toolbox.
% -------------------------------------------------------------------------
clear; clc;

%% ----------------------------- Settings ---------------------------------
config;
cfg.resDir = fullfile(OUT_DIR, 'results_nback');
cfg.outDir = fullfile(cfg.resDir, 'stats');
cfg.blkPrimary   = {'dprime_WM','Accuracy_WM','HitRT_median_WM', ...
                    'LoadSlope_dprime','LoadSlope_HitRT_median'};
cfg.blkSecondary = {'dprime','c','Accuracy','HitRT_median','HitRT_CV','LoadSlope_Accuracy'};
cfg.lvPrimary    = {'dprime','Accuracy','HitRT_median','ISIHitProp'};
cfg.lvSecondary  = {'HitRate','FARate','c','HitRT_stim_median','HitRT_ISI_median','HitRT_CV'};
cfg.tcMetrics    = {'dprime','Accuracy','HitRT_median','ISIHitProp'};
cfg.tcLoads      = 1:3;     % loads used for the time-course model (0-back is at ceiling)
cfg.minCellN     = 3;       % cells with fewer observations are left out
cfg.makeFigures  = true;
ints  = {'E','N','R'};
loads = 0:3; reps = 1:3;

if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end
CB = readtable(fullfile(cfg.resDir, 'Nback_change_scores_block.csv'));
CL = readtable(fullfile(cfg.resDir, 'Nback_change_scores_level.csv'));
CS = readtable(fullfile(cfg.resDir, 'Nback_change_scores_subblock.csv'));
LV = readtable(fullfile(cfg.resDir, 'Nback_level_metrics.csv'));
fprintf('Change scores: %d (block), %d (load), %d (sub-block); %d subjects.\n', ...
    height(CB), height(CL), height(CS), numel(unique(CB.Subject)));

%% ============================ A) overall ================================
A = struct('omni', {{}}, 'pair', {{}}, 'mean', {{}}, 'cov', {{}});
blkAll = [cfg.blkPrimary, cfg.blkSecondary];
blkFam = [repmat({'Primary'}, 1, numel(cfg.blkPrimary)), repmat({'Secondary'}, 1, numel(cfg.blkSecondary))];
for k = 1:numel(blkAll)
    m = blkAll{k}; [sc, unit] = scaleOf(m);
    D = prep(CB, CB.(['delta_' m]) * sc);
    D.cell = string(D.Intervention);
    [lme, nm, D] = fitCells(D, cellist(ints), {'Position','OrderPre','OrderPost'}, false, cfg.minCellN);
    if isempty(lme), continue; end
    H = [w(nm, {'E','N'}, [1 -1]); w(nm, {'N','R'}, [1 -1])];
    o = ctest(lme, H);
    A.omni{end+1} = struct('Metric', m, 'Family', blkFam{k}, 'Unit', unit, 'nObs', height(D), ...
        'nSubj', numel(unique(D.Subject)), 'F', o.F, 'df1', o.df1, 'df2', o.df2, 'p', o.p, ...
        'p_FDR', NaN, 'eta2_partial', o.eta2); %#ok<SAGROW>
    sdT = sdTotal(lme);
    pr = {'E','N'; 'E','R'; 'N','R'}; pp = nan(3,1); rr = cell(3,1);
    for q = 1:3
        c = ctest(lme, w(nm, pr(q,:), [1 -1]));
        pp(q) = c.p;
        rr{q} = struct('Metric', m, 'Family', blkFam{k}, 'Unit', unit, 'Contrast', [pr{q,1} '-' pr{q,2}], ...
            'Estimate', c.est, 'SE', c.se, 'CI_low', c.ci(1), 'CI_high', c.ci(2), 'df', c.df2, ...
            'p', c.p, 'p_Holm', NaN, 'd', c.est / sdT);
    end
    ph = holm(pp);
    for q = 1:3, rr{q}.p_Holm = ph(q); A.pair{end+1} = rr{q}; end
    for q = 1:3
        c = ctest(lme, w(nm, ints(q), 1));
        A.mean{end+1} = struct('Metric', m, 'Family', blkFam{k}, 'Unit', unit, 'Intervention', ints{q}, ...
            'n', sum(D.cell == ints{q}), 'AdjMeanChange', c.est, 'SE', c.se, 'CI_low', c.ci(1), ...
            'CI_high', c.ci(2), 'p_vs0', c.p); %#ok<SAGROW>
    end
    for cv = {'Position','OrderPre','OrderPost'}
        c = ctest(lme, w(nm, cv, 1));
        A.cov{end+1} = struct('Metric', m, 'Unit', unit, 'Covariate', cv{1}, 'Estimate', c.est, ...
            'SE', c.se, 'df', c.df2, 'p', c.p); %#ok<SAGROW>
    end
end
A_omni = fdrBy(struct2table([A.omni{:}]), 'p', 'p_FDR', 'Family');
A_pair = struct2table([A.pair{:}]);
A_mean = struct2table([A.mean{:}]);
A_cov  = struct2table([A.cov{:}]);

%% ======================= B) per load level ==============================
B = struct('main', {{}}, 'simple', {{}}, 'pair', {{}}, 'mean', {{}});
lvAll = [cfg.lvPrimary, cfg.lvSecondary];
lvFam = [repmat({'Primary'}, 1, numel(cfg.lvPrimary)), repmat({'Secondary'}, 1, numel(cfg.lvSecondary))];
cellsB = {};
for i = 1:3, for L = loads, cellsB{end+1} = sprintf('%s_L%d', ints{i}, L); end, end %#ok<SAGROW>
for k = 1:numel(lvAll)
    m = lvAll{k}; [sc, unit] = scaleOf(m);
    D = prep(CL, CL.(['delta_' m]) * sc);
    D.cell = string(D.Intervention) + "_L" + string(D.Load);
    [lme, nm, D] = fitCells(D, cellsB, {'Position','OrderPre','OrderPost'}, true, cfg.minCellN);
    if isempty(lme), continue; end
    cn = @(i, L) sprintf('%s_L%d', i, L);
    % main effect of intervention (average over loads)
    H = [avgDiff(nm, cn, 'E', 'N', loads); avgDiff(nm, cn, 'N', 'R', loads)];
    tI = ctest(lme, H);
    % main effect of load (average over interventions)
    H = [];
    for L = loads(2:end), H = [H; loadDiff(nm, cn, ints, L, loads(1))]; end %#ok<AGROW>
    tL = ctest(lme, H);
    % interaction
    H = [];
    for L = loads(2:end)
        H = [H; w(nm, {cn('E',L), cn('N',L), cn('E',loads(1)), cn('N',loads(1))}, [1 -1 -1 1]); ...
                w(nm, {cn('N',L), cn('R',L), cn('N',loads(1)), cn('R',loads(1))}, [1 -1 -1 1])]; %#ok<AGROW>
    end
    tX = ctest(lme, H);
    B.main{end+1} = struct('Metric', m, 'Family', lvFam{k}, 'Unit', unit, 'nObs', height(D), ...
        'nSubj', numel(unique(D.Subject)), ...
        'F_Intervention', tI.F, 'df_Intervention', tI.df1, 'df2_Intervention', tI.df2, 'p_Intervention', tI.p, 'pFDR_Intervention', NaN, ...
        'F_Load', tL.F, 'df_Load', tL.df1, 'df2_Load', tL.df2, 'p_Load', tL.p, ...
        'F_IxLoad', tX.F, 'df_IxLoad', tX.df1, 'df2_IxLoad', tX.df2, 'p_IxLoad', tX.p, 'pFDR_IxLoad', NaN); %#ok<SAGROW>
    sdT = sdTotal(lme);
    % simple effects at each load
    for L = loads
        H = [w(nm, {cn('E',L), cn('N',L)}, [1 -1]); w(nm, {cn('N',L), cn('R',L)}, [1 -1])];
        s = ctest(lme, H);
        B.simple{end+1} = struct('Metric', m, 'Family', lvFam{k}, 'Unit', unit, 'Load', L, ...
            'F', s.F, 'df1', s.df1, 'df2', s.df2, 'p', s.p, 'p_FDR_acrossLoads', NaN); %#ok<SAGROW>
        pr = {'E','N'; 'E','R'; 'N','R'}; pp = nan(3,1); rr = cell(3,1);
        for q = 1:3
            c = ctest(lme, w(nm, {cn(pr{q,1},L), cn(pr{q,2},L)}, [1 -1]));
            pp(q) = c.p;
            rr{q} = struct('Metric', m, 'Family', lvFam{k}, 'Unit', unit, 'Load', L, ...
                'Contrast', [pr{q,1} '-' pr{q,2}], 'Estimate', c.est, 'SE', c.se, ...
                'CI_low', c.ci(1), 'CI_high', c.ci(2), 'df', c.df2, 'p', c.p, 'p_Holm', NaN, 'd', c.est / sdT);
        end
        ph = holm(pp);
        for q = 1:3, rr{q}.p_Holm = ph(q); B.pair{end+1} = rr{q}; end
        for q = 1:3
            c = ctest(lme, w(nm, {cn(ints{q},L)}, 1));
            B.mean{end+1} = struct('Metric', m, 'Family', lvFam{k}, 'Unit', unit, 'Load', L, ...
                'Intervention', ints{q}, 'n', sum(D.cell == cn(ints{q},L)), 'AdjMeanChange', c.est, ...
                'SE', c.se, 'CI_low', c.ci(1), 'CI_high', c.ci(2), 'p_vs0', c.p); %#ok<SAGROW>
        end
    end
end
B_main = struct2table([B.main{:}]);
B_main = fdrBy(B_main, 'p_Intervention', 'pFDR_Intervention', 'Family');
B_main = fdrBy(B_main, 'p_IxLoad', 'pFDR_IxLoad', 'Family');
B_simple = fdrBy(struct2table([B.simple{:}]), 'p', 'p_FDR_acrossLoads', 'Metric');
B_pair = struct2table([B.pair{:}]);
B_mean = struct2table([B.mean{:}]);

%% ======================== C) time course ================================
C = struct('main', {{}}, 'simple', {{}}, 'mean', {{}});
cellsC = {};
for i = 1:3, for r = reps, cellsC{end+1} = sprintf('%s_R%d', ints{i}, r); end, end %#ok<SAGROW>
CSu = CS(ismember(CS.Load, cfg.tcLoads), :);
for k = 1:numel(cfg.tcMetrics)
    m = cfg.tcMetrics{k}; [sc, unit] = scaleOf(m);
    D = prep(CSu, CSu.(['delta_' m]) * sc);
    D.cell = string(D.Intervention) + "_R" + string(D.Repetition);
    covs = {'Position','OrderPre','OrderPost'};
    for L = cfg.tcLoads(2:end)
        vn = sprintf('Load%d', L);
        D.(vn) = double(D.Load == L); D.(vn) = D.(vn) - mean(D.(vn));
        covs{end+1} = vn; %#ok<SAGROW>
    end
    [lme, nm, D] = fitCells(D, cellsC, covs, true, cfg.minCellN);
    if isempty(lme), continue; end
    cn = @(i, r) sprintf('%s_R%d', i, r);
    tI = ctest(lme, [avgDiff(nm, cn, 'E', 'N', reps); avgDiff(nm, cn, 'N', 'R', reps)]);
    H = [];
    for r = reps(2:end), H = [H; loadDiff(nm, cn, ints, r, reps(1))]; end %#ok<AGROW>
    tR = ctest(lme, H);
    H = [];
    for r = reps(2:end)
        H = [H; w(nm, {cn('E',r), cn('N',r), cn('E',1), cn('N',1)}, [1 -1 -1 1]); ...
                w(nm, {cn('N',r), cn('R',r), cn('N',1), cn('R',1)}, [1 -1 -1 1])]; %#ok<AGROW>
    end
    tX = ctest(lme, H);
    C.main{end+1} = struct('Metric', m, 'Unit', unit, 'Loads', mat2str(cfg.tcLoads), 'nObs', height(D), ...
        'F_Intervention', tI.F, 'df2_Intervention', tI.df2, 'p_Intervention', tI.p, ...
        'F_Repetition', tR.F, 'df2_Repetition', tR.df2, 'p_Repetition', tR.p, ...
        'F_IxRep', tX.F, 'df2_IxRep', tX.df2, 'p_IxRep', tX.p, 'pFDR_IxRep', NaN); %#ok<SAGROW>
    for r = reps
        s = ctest(lme, [w(nm, {cn('E',r), cn('N',r)}, [1 -1]); w(nm, {cn('N',r), cn('R',r)}, [1 -1])]);
        C.simple{end+1} = struct('Metric', m, 'Unit', unit, 'Repetition', r, 'F', s.F, ...
            'df2', s.df2, 'p', s.p); %#ok<SAGROW>
        for q = 1:3
            c = ctest(lme, w(nm, {cn(ints{q},r)}, 1));
            C.mean{end+1} = struct('Metric', m, 'Unit', unit, 'Repetition', r, 'Intervention', ints{q}, ...
                'AdjMeanChange', c.est, 'SE', c.se, 'CI_low', c.ci(1), 'CI_high', c.ci(2), 'p_vs0', c.p); %#ok<SAGROW>
        end
    end
end
C_main = struct2table([C.main{:}]);
C_main.pFDR_IxRep = bh(C_main.p_IxRep);
C_simple = struct2table([C.simple{:}]);
C_mean = struct2table([C.mean{:}]);

%% ===================== D) manipulation check ============================
Dk = {};
LVv = LV(LV.BlockValid == 1, :);
LVv.SubjectCat = categorical(LVv.Subject);
LVv.LoadCat  = categorical(LVv.Load);
LVv.BlockCat = categorical(LVv.Block);
for m = {'dprime','Accuracy','HitRT_median','ISIHitProp'}
    [sc, unit] = scaleOf(m{1});
    T = LVv; T.y = T.(m{1}) * sc; T = T(~isnan(T.y), :);
    lme = fitlme(T, 'y ~ LoadCat + BlockCat + (1|SubjectCat)', 'FitMethod', 'REML');
    a = anova(lme, 'DFMethod', 'satterthwaite');
    for t = {'LoadCat','BlockCat'}
        i = find(strcmp(a.Term, t{1}));
        mu = arrayfun(@(L) mean(T.y(T.Load == L)), loads);
        Dk{end+1} = struct('Metric', m{1}, 'Unit', unit, 'Term', t{1}, 'F', a.FStat(i), ...
            'df1', a.DF1(i), 'df2', a.DF2(i), 'p', a.pValue(i), ...
            'Mean_0back', mu(1), 'Mean_1back', mu(2), 'Mean_2back', mu(3), 'Mean_3back', mu(4)); %#ok<SAGROW>
    end
end
D_check = struct2table([Dk{:}]);

%% -------------------------------- Save -----------------------------------
tabs = {'A_overall_omnibus', A_omni; 'A_overall_pairwise', A_pair; 'A_overall_means', A_mean; ...
        'A_overall_covariates', A_cov; 'B_load_main', B_main; 'B_load_simple', B_simple; ...
        'B_load_pairwise', B_pair; 'B_load_means', B_mean; 'C_timecourse_main', C_main; ...
        'C_timecourse_simple', C_simple; 'C_timecourse_means', C_mean; 'D_load_check', D_check};
xf = fullfile(cfg.outDir, 'Nback_stats.xlsx');
if exist(xf, 'file'), delete(xf); end
for t = 1:size(tabs, 1)
    writetable(tabs{t,2}, fullfile(cfg.outDir, [tabs{t,1} '.csv']));
    try, writetable(tabs{t,2}, xf, 'Sheet', tabs{t,1}(1:min(31, end))); catch, end
end
save(fullfile(cfg.outDir, 'Nback_stats.mat'), 'A_omni', 'A_pair', 'A_mean', 'A_cov', ...
     'B_main', 'B_simple', 'B_pair', 'B_mean', 'C_main', 'C_simple', 'C_mean', 'D_check', 'cfg');

%% ------------------------------- Report ----------------------------------
fid = fopen(fullfile(cfg.outDir, 'Nback_stats_report.txt'), 'w');
o = @(varargin) fprintf(fid, varargin{:});
o('N-back - intervention effects (LMM on change scores, Satterthwaite df)\n\n');
o('D) MANIPULATION CHECK (raw scores): y ~ Load + Block + (1|Subject)\n');
for i = 1:height(D_check)
    o('  %-14s %-9s F(%d, %.1f) = %.2f, p = %.3g   [means 0/1/2/3-back: %.3g / %.3g / %.3g / %.3g]\n', ...
        D_check.Metric{i}, D_check.Term{i}, D_check.df1(i), D_check.df2(i), D_check.F(i), D_check.p(i), ...
        D_check.Mean_0back(i), D_check.Mean_1back(i), D_check.Mean_2back(i), D_check.Mean_3back(i));
end
o('\nA) OVERALL / COMPOSITE INDICES: delta ~ E + N + R + covariates + (1|Subject)\n');
for i = 1:height(A_omni)
    mm = A_mean(strcmp(A_mean.Metric, A_omni.Metric{i}), :);
    o('  %-24s [%s] F(%d, %.1f) = %.2f, p = %.3f, p_FDR = %.3f | change E %+.3g, N %+.3g, R %+.3g %s\n', ...
        A_omni.Metric{i}, A_omni.Family{i}(1:3), A_omni.df1(i), A_omni.df2(i), A_omni.F(i), A_omni.p(i), ...
        A_omni.p_FDR(i), mm.AdjMeanChange(1), mm.AdjMeanChange(2), mm.AdjMeanChange(3), A_omni.Unit{i});
end
o('\nB) PER LOAD: delta ~ Intervention x Load + covariates + (1|Subject) + (1|Subject:PostBlock)\n');
for i = 1:height(B_main)
    o('  %-18s [%s] Intervention F = %.2f, p = %.3f (FDR %.3f) | Load F = %.2f, p = %.3g | IxLoad F(%d, %.1f) = %.2f, p = %.3f (FDR %.3f)\n', ...
        B_main.Metric{i}, B_main.Family{i}(1:3), B_main.F_Intervention(i), B_main.p_Intervention(i), ...
        B_main.pFDR_Intervention(i), B_main.F_Load(i), B_main.p_Load(i), B_main.df_IxLoad(i), ...
        B_main.df2_IxLoad(i), B_main.F_IxLoad(i), B_main.p_IxLoad(i), B_main.pFDR_IxLoad(i));
    ss = B_simple(strcmp(B_simple.Metric, B_main.Metric{i}), :);
    for j = 1:height(ss)
        mm = B_mean(strcmp(B_mean.Metric, B_main.Metric{i}) & B_mean.Load == ss.Load(j), :);
        o('      %d-back: F = %.2f, p = %.3f | change E %+.3g, N %+.3g, R %+.3g\n', ss.Load(j), ss.F(j), ss.p(j), ...
            mm.AdjMeanChange(1), mm.AdjMeanChange(2), mm.AdjMeanChange(3));
    end
end
o('\nC) TIME COURSE (loads %s): delta ~ Intervention x Repetition + Load + covariates + random\n', mat2str(cfg.tcLoads));
for i = 1:height(C_main)
    o('  %-14s Intervention p = %.3f | Repetition p = %.3f | IxRep F = %.2f, p = %.3f (FDR %.3f)\n', ...
        C_main.Metric{i}, C_main.p_Intervention(i), C_main.p_Repetition(i), C_main.F_IxRep(i), ...
        C_main.p_IxRep(i), C_main.pFDR_IxRep(i));
end
fclose(fid);
type(fullfile(cfg.outDir, 'Nback_stats_report.txt'));

%% ------------------------------- Figures ---------------------------------
if cfg.makeFigures
    cols = [0.85 0.33 0.10; 0.30 0.45 0.69; 0.47 0.67 0.19];
    names = {'Reels (E)','Neutral (N)','Rest (R)'};
    % (1) composites
    fig = figure('Color', 'w', 'Position', [80 80 1300 650], 'Visible', 'off');
    for p = 1:numel(cfg.blkPrimary)
        m = cfg.blkPrimary{p}; [sc, unit] = scaleOf(m);
        subplot(2, 3, p); hold on
        mm = A_mean(strcmp(A_mean.Metric, m), :);
        for q = 1:3
            y = CB.(['delta_' m])(strcmp(CB.Intervention, ints{q})) * sc; y = y(~isnan(y));
            scatter(q + (rand(size(y)) - 0.5) * 0.35, y, 14, cols(q,:), 'filled', 'MarkerFaceAlpha', 0.35);
            errorbar(q, mm.AdjMeanChange(q), mm.AdjMeanChange(q) - mm.CI_low(q), mm.CI_high(q) - mm.AdjMeanChange(q), ...
                'o', 'Color', cols(q,:) * 0.7, 'MarkerFaceColor', cols(q,:), 'LineWidth', 1.6, 'MarkerSize', 7);
        end
        plot([0.4 3.6], [0 0], 'k:');
        om = A_omni(strcmp(A_omni.Metric, m), :);
        lab = strrep(m, '_', ' '); if ~isempty(unit), lab = [lab ' (' unit ')']; end
        title(sprintf('\\Delta %s\np = %.3f (FDR %.3f)', lab, om.p, om.p_FDR), 'FontWeight', 'normal');
        set(gca, 'XTick', 1:3, 'XTickLabel', names, 'Box', 'off'); xlim([0.4 3.6]);
    end
    print(fig, fullfile(cfg.outDir, 'Nback_change_composites.png'), '-dpng', '-r150'); close(fig);
    % (2) by load
    fig = figure('Color', 'w', 'Position', [80 80 1300 650], 'Visible', 'off');
    for p = 1:numel(cfg.lvPrimary)
        m = cfg.lvPrimary{p}; [~, unit] = scaleOf(m);
        subplot(2, 2, p); hold on
        for q = 1:3
            mm = B_mean(strcmp(B_mean.Metric, m) & strcmp(B_mean.Intervention, ints{q}), :);
            mm = sortrows(mm, 'Load');
            errorbar(mm.Load + (q - 2) * 0.08, mm.AdjMeanChange, mm.AdjMeanChange - mm.CI_low, ...
                mm.CI_high - mm.AdjMeanChange, '-o', 'Color', cols(q,:), 'MarkerFaceColor', cols(q,:), 'LineWidth', 1.4);
        end
        plot([-0.4 3.4], [0 0], 'k:');
        bm = B_main(strcmp(B_main.Metric, m), :);
        lab = strrep(m, '_', ' '); if ~isempty(unit), lab = [lab ' (' unit ')']; end
        title(sprintf('\\Delta %s  |  Intervention p = %.3f, I x Load p = %.3f', lab, bm.p_Intervention, bm.p_IxLoad), ...
            'FontWeight', 'normal');
        set(gca, 'XTick', loads, 'XTickLabel', {'0-back','1-back','2-back','3-back'}, 'Box', 'off');
        xlim([-0.4 3.4]); if p == 1, legend(names, 'Location', 'best', 'Box', 'off'); end
    end
    print(fig, fullfile(cfg.outDir, 'Nback_change_by_load.png'), '-dpng', '-r150'); close(fig);
    % (3) time course
    fig = figure('Color', 'w', 'Position', [80 80 1300 650], 'Visible', 'off');
    for p = 1:numel(cfg.tcMetrics)
        m = cfg.tcMetrics{p}; [~, unit] = scaleOf(m);
        subplot(2, 2, p); hold on
        for q = 1:3
            mm = sortrows(C_mean(strcmp(C_mean.Metric, m) & strcmp(C_mean.Intervention, ints{q}), :), 'Repetition');
            errorbar(mm.Repetition + (q - 2) * 0.06, mm.AdjMeanChange, mm.AdjMeanChange - mm.CI_low, ...
                mm.CI_high - mm.AdjMeanChange, '-o', 'Color', cols(q,:), 'MarkerFaceColor', cols(q,:), 'LineWidth', 1.4);
        end
        plot([0.6 3.4], [0 0], 'k:');
        cm = C_main(strcmp(C_main.Metric, m), :);
        lab = strrep(m, '_', ' '); if ~isempty(unit), lab = [lab ' (' unit ')']; end
        title(sprintf('\\Delta %s  |  I x Repetition p = %.3f', lab, cm.p_IxRep), 'FontWeight', 'normal');
        set(gca, 'XTick', reps, 'XTickLabel', {'Cycle 1','Cycle 2','Cycle 3'}, 'Box', 'off');
        xlim([0.6 3.4]); if p == 1, legend(names, 'Location', 'best', 'Box', 'off'); end
    end
    print(fig, fullfile(cfg.outDir, 'Nback_change_timecourse.png'), '-dpng', '-r150'); close(fig);
end
fprintf('\nDone. Statistics saved in:\n%s\n', cfg.outDir);


%% =========================================================================
%%                              Local functions
%% =========================================================================
function [sc, unit] = scaleOf(m)
if contains(m, 'RT') && ~contains(m, 'CV'), sc = 1000; unit = 'ms'; else, sc = 1; unit = ''; end
end

function c = cellist(x), c = x; end

function D = prep(T, y)
% keep valid outcome rows; centred covariates; grouping variables
ok = ~isnan(y);
D = T(ok, :); D.y = y(ok);
D.SubjectCat = categorical(D.Subject);
D.UnitCat    = categorical(D.Subject * 10 + D.PostBlock);
D.Position  = D.InterventionPosition - mean(D.InterventionPosition);
D.OrderPre  = double(string(D.TaskOrderPre)  == "SN"); D.OrderPre  = D.OrderPre  - mean(D.OrderPre);
D.OrderPost = double(string(D.TaskOrderPost) == "SN"); D.OrderPost = D.OrderPost - mean(D.OrderPost);
end

function [lme, nm, D] = fitCells(D, cells, covs, unitRE, minN)
% Cell-means mixed model: one dummy per cell (cells with < minN obs dropped).
keep = false(height(D), 1); used = {};
for i = 1:numel(cells)
    s = D.cell == cells{i};
    if sum(s) >= minN
        used{end+1} = cells{i}; keep = keep | s; %#ok<AGROW>
    end
end
D = D(keep, :);
lme = []; nm = {};
if numel(used) < 2, return; end
terms = {};
for i = 1:numel(used)
    vn = ['c_' used{i}];
    D.(vn) = double(D.cell == used{i});
    terms{end+1} = vn; %#ok<AGROW>
end
f = ['y ~ -1 + ' strjoin([terms, covs], ' + ') ' + (1|SubjectCat)'];
if unitRE, f = [f ' + (1|UnitCat)']; end
lme = fitlme(D, f, 'FitMethod', 'REML');
nm  = strrep(lme.CoefficientNames, 'c_', '');
end

function h = w(nm, labels, weights)
% contrast row; NaN row if a needed cell is missing
h = zeros(1, numel(nm));
for i = 1:numel(labels)
    j = find(strcmp(nm, labels{i}));
    if isempty(j), h = nan(1, numel(nm)); return; end
    h(j) = h(j) + weights(i);
end
end

function h = avgDiff(nm, cn, a, b, levels)
% (mean over levels of a) - (mean over levels of b)
lab = {}; wt = [];
for L = levels
    lab = [lab, {cn(a, L), cn(b, L)}]; wt = [wt, 1 / numel(levels), -1 / numel(levels)]; %#ok<AGROW>
end
h = w(nm, lab, wt);
end

function h = loadDiff(nm, cn, ints, L, L0)
% (mean over interventions at L) - (mean over interventions at L0)
lab = {}; wt = [];
for i = 1:numel(ints)
    lab = [lab, {cn(ints{i}, L), cn(ints{i}, L0)}]; wt = [wt, 1 / numel(ints), -1 / numel(ints)]; %#ok<AGROW>
end
h = w(nm, lab, wt);
end

function r = ctest(lme, H)
% Satterthwaite F / t test of H*beta = 0 (rows with missing cells dropped)
H = H(~any(isnan(H), 2), :);
r = struct('F', NaN, 'df1', NaN, 'df2', NaN, 'p', NaN, 'eta2', NaN, 'est', NaN, 'se', NaN, 'ci', [NaN NaN]);
if isempty(H), return; end
[r.p, r.F, r.df1, r.df2] = coefTest(lme, H, zeros(size(H, 1), 1), 'DFMethod', 'satterthwaite');
r.eta2 = r.F * r.df1 / (r.F * r.df1 + r.df2);
if size(H, 1) == 1
    b = fixedEffects(lme); V = lme.CoefficientCovariance;
    r.est = H * b; r.se = sqrt(H * V * H');
    r.ci = r.est + [-1 1] * tinv(0.975, r.df2) * r.se;
end
end

function s = sdTotal(lme)
[psi, mse] = covarianceParameters(lme);
v = mse;
for i = 1:numel(psi), v = v + psi{i}(1); end
s = sqrt(v);
end

function T = fdrBy(T, pcol, qcol, grp)
g = T.(grp);
if ~iscell(g), g = cellstr(string(g)); end
for u = unique(g)'
    s = strcmp(g, u{1});
    T.(qcol)(s) = bh(T.(pcol)(s));
end
end

function q = bh(p)
p = p(:); q = nan(size(p));
v = find(~isnan(p));
[ps, o] = sort(p(v)); m = numel(ps);
qs = min(1, cummin(ps .* m ./ (1:m)', 'reverse'));
q(v(o)) = qs;
end

function q = holm(p)
p = p(:); q = nan(size(p));
v = find(~isnan(p));
[ps, o] = sort(p(v)); m = numel(ps);
qs = min(1, cummax((m - (1:m)' + 1) .* ps));
q(v(o)) = qs;
end
