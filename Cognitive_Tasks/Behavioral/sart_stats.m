%% SART statistics: effect of the three interventions (E, N, R)
% -------------------------------------------------------------------------
% Input : results of sart_analysis.m (results_sart/SART_change_scores.csv)
% Output: results_sart/stats/
%
% Primary model (one per index), linear mixed model on change scores:
%   delta = post block - pre block
%   delta ~ -1 + E + N + R + Position + OrderPre + OrderPost + (1 | Subject)
%     E, N, R    : intervention dummies (no intercept, so each coefficient is
%                  the adjusted mean change after that intervention)
%     Position   : intervention 1st or 2nd (centred)
%     OrderPre   : SART first (SN) in the pre block (centred)
%     OrderPost  : SART first (SN) in the post block (centred)
%   Tests (Satterthwaite df):
%     - omnibus: E = N = R
%     - pairwise: E-N, E-R, N-R  (Holm-corrected within each index)
%     - each intervention vs 0 (did the index change?)
%   Omnibus p-values are FDR-corrected (Benjamini-Hochberg) within the
%   primary and the secondary family of indices.
%
% Sensitivity model (primary indices): ANCOVA-type mixed model
%   post ~ -1 + E + N + R + Pre + Position + OrderPre + OrderPost + (1|Subject)
%
% RT-based indices are converted to milliseconds before modelling.
% Requirements: MATLAB R2018b+ with Statistics and Machine Learning Toolbox.
% -------------------------------------------------------------------------
clear; clc;

%% ----------------------------- Settings ---------------------------------
config;
cfg.resDir  = fullfile(OUT_DIR, 'results_sart');
cfg.outDir  = fullfile(cfg.resDir, 'stats');
cfg.primary = {'CommissionRate','OmissionRate','dprime','c', ...
               'GoRT_median','GoRT_CV','CommRT_median','ExG_tau'};
cfg.secondary = {'GoRT_mean','GoRT_SD','ExG_mu','ExG_sigma','GoMinusCommRT', ...
                 'PreErrorSpeeding','PES_traditional','PES_robust'};
cfg.msMetrics = {'GoRT_median','GoRT_mean','GoRT_SD','CommRT_median', ...
                 'GoMinusCommRT','ExG_mu','ExG_sigma','ExG_tau', ...
                 'PreErrorSpeeding','PES_traditional','PES_robust'};
cfg.runANCOVA   = true;
cfg.makeFigures = true;
ints = {'E','N','R'};

if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end
T = readtable(fullfile(cfg.resDir, 'SART_change_scores.csv'));
T.SubjectCat = categorical(T.Subject);
fprintf('%d change scores from %d subjects.\n', height(T), numel(unique(T.Subject)));

%% ------------------------------ Models -----------------------------------
allMetrics = [cfg.primary, cfg.secondary];
family     = [repmat({'Primary'}, 1, numel(cfg.primary)), ...
              repmat({'Secondary'}, 1, numel(cfg.secondary))];

omniRows = {}; pairRows = {}; meanRows = {}; covRows = {};
ancOmni = {}; ancPair = {};
fitted = struct();

for k = 1:numel(allMetrics)
    m = allMetrics{k};
    sc = 1; unit = '';
    if ismember(m, cfg.msMetrics), sc = 1000; unit = 'ms'; end

    % ---------- change-score model ----------
    D = prepData(T, T.(['delta_' m]) * sc, []);
    [lme, R] = fitAndTest(D, false, ints);
    fitted.(m) = lme;

    o = struct('Metric', m, 'Family', family{k}, 'Unit', unit, ...
        'nObs', height(D), 'nSubj', numel(unique(D.Subject)), ...
        'F', R.F, 'df1', R.df1, 'df2', R.df2, 'p', R.p, 'p_FDR', NaN, ...
        'eta2_partial', R.F * R.df1 / (R.F * R.df1 + R.df2));
    omniRows{end+1} = o; %#ok<SAGROW>

    for q = 1:3
        meanRows{end+1} = struct('Metric', m, 'Family', family{k}, 'Unit', unit, ...
            'Intervention', ints{q}, 'n', sum(D.(ints{q})), ...
            'AdjMeanChange', R.mean(q), 'SE', R.meanSE(q), ...
            'CI_low', R.meanCI(q,1), 'CI_high', R.meanCI(q,2), ...
            'df', R.meanDF(q), 'p_vs0', R.meanP(q), ...
            'RawMeanChange', mean(D.y(D.(ints{q}) == 1))); %#ok<SAGROW>
    end
    pH = holm(R.pairP);
    for q = 1:3
        pairRows{end+1} = struct('Metric', m, 'Family', family{k}, 'Unit', unit, ...
            'Contrast', R.pairName{q}, 'Estimate', R.pairEst(q), 'SE', R.pairSE(q), ...
            'CI_low', R.pairCI(q,1), 'CI_high', R.pairCI(q,2), 'df', R.pairDF(q), ...
            't', R.pairEst(q) / R.pairSE(q), 'p', R.pairP(q), 'p_Holm', pH(q), ...
            'd', R.pairEst(q) / R.sdTotal); %#ok<SAGROW>
    end
    for q = 1:numel(R.covName)
        covRows{end+1} = struct('Metric', m, 'Unit', unit, 'Covariate', R.covName{q}, ...
            'Estimate', R.covEst(q), 'SE', R.covSE(q), 'df', R.covDF(q), ...
            'p', R.covP(q)); %#ok<SAGROW>
    end

    % ---------- ANCOVA sensitivity (primary indices) ----------
    if cfg.runANCOVA && strcmp(family{k}, 'Primary')
        D2 = prepData(T, T.(['post_' m]) * sc, T.(['pre_' m]) * sc);
        [~, R2] = fitAndTest(D2, true, ints);
        ancOmni{end+1} = struct('Metric', m, 'Unit', unit, 'nObs', height(D2), ...
            'F', R2.F, 'df1', R2.df1, 'df2', R2.df2, 'p', R2.p, 'p_FDR', NaN); %#ok<SAGROW>
        pH2 = holm(R2.pairP);
        for q = 1:3
            ancPair{end+1} = struct('Metric', m, 'Unit', unit, 'Contrast', R2.pairName{q}, ...
                'Estimate', R2.pairEst(q), 'SE', R2.pairSE(q), ...
                'CI_low', R2.pairCI(q,1), 'CI_high', R2.pairCI(q,2), ...
                'p', R2.pairP(q), 'p_Holm', pH2(q)); %#ok<SAGROW>
        end
    end
end

omniT = struct2table([omniRows{:}]);
for f = {'Primary','Secondary'}
    s = strcmp(omniT.Family, f{1});
    omniT.p_FDR(s) = bh(omniT.p(s));
end
pairT = struct2table([pairRows{:}]);
meanT = struct2table([meanRows{:}]);
covT  = struct2table([covRows{:}]);
if ~isempty(ancOmni)
    ancOmniT = struct2table([ancOmni{:}]);
    ancOmniT.p_FDR = bh(ancOmniT.p);
    ancPairT = struct2table([ancPair{:}]);
else
    ancOmniT = table(); ancPairT = table();
end

%% -------------------------------- Save -----------------------------------
writetable(omniT, fullfile(cfg.outDir, 'LMM_omnibus.csv'));
writetable(pairT, fullfile(cfg.outDir, 'LMM_pairwise.csv'));
writetable(meanT, fullfile(cfg.outDir, 'LMM_adjusted_means.csv'));
writetable(covT,  fullfile(cfg.outDir, 'LMM_covariates.csv'));
writetable(ancOmniT, fullfile(cfg.outDir, 'ANCOVA_omnibus.csv'));
writetable(ancPairT, fullfile(cfg.outDir, 'ANCOVA_pairwise.csv'));
xf = fullfile(cfg.outDir, 'SART_stats.xlsx');
if exist(xf, 'file'), delete(xf); end
try
    writetable(omniT, xf, 'Sheet', 'Omnibus');
    writetable(pairT, xf, 'Sheet', 'Pairwise');
    writetable(meanT, xf, 'Sheet', 'AdjustedMeans');
    writetable(covT,  xf, 'Sheet', 'Covariates');
    if ~isempty(ancOmniT)
        writetable(ancOmniT, xf, 'Sheet', 'ANCOVA_Omnibus');
        writetable(ancPairT, xf, 'Sheet', 'ANCOVA_Pairwise');
    end
catch ME
    warning('xlsx not written: %s', ME.message);
end
save(fullfile(cfg.outDir, 'SART_stats.mat'), 'omniT', 'pairT', 'meanT', 'covT', ...
     'ancOmniT', 'ancPairT', 'fitted', 'cfg');

%% ------------------------------- Report ----------------------------------
fid = fopen(fullfile(cfg.outDir, 'SART_stats_report.txt'), 'w');
out = @(varargin) fprintf(fid, varargin{:});
out('SART - intervention effects (LMM on change scores)\n');
out('Model: delta ~ -1 + E + N + R + Position + OrderPre + OrderPost + (1|Subject)\n');
out('%d change scores, %d subjects. Satterthwaite df.\n\n', height(T), numel(unique(T.Subject)));
for k = 1:height(omniT)
    m = omniT.Metric{k};
    out('%-18s [%s] F(%d, %.1f) = %.2f, p = %.3f, p_FDR = %.3f, eta2p = %.3f\n', ...
        m, omniT.Family{k}, omniT.df1(k), omniT.df2(k), omniT.F(k), omniT.p(k), ...
        omniT.p_FDR(k), omniT.eta2_partial(k));
    s = strcmp(meanT.Metric, m);
    mm = meanT(s, :);
    for q = 1:height(mm)
        out('      %s: change = %+.4g %s [%.4g, %.4g], p vs 0 = %.3f\n', mm.Intervention{q}, ...
            mm.AdjMeanChange(q), mm.Unit{q}, mm.CI_low(q), mm.CI_high(q), mm.p_vs0(q));
    end
    s = strcmp(pairT.Metric, m);
    pp = pairT(s, :);
    for q = 1:height(pp)
        out('      %s: %+.4g, p = %.3f, p_Holm = %.3f, d = %.2f\n', pp.Contrast{q}, ...
            pp.Estimate(q), pp.p(q), pp.p_Holm(q), pp.d(q));
    end
end
fclose(fid);
type(fullfile(cfg.outDir, 'SART_stats_report.txt'));

%% ------------------------------- Figure ----------------------------------
if cfg.makeFigures
    cols = [0.85 0.33 0.10; 0.30 0.45 0.69; 0.47 0.67 0.19];
    names = {'Reels (E)','Neutral (N)','Rest (R)'};
    fig = figure('Color', 'w', 'Position', [80 80 1400 680], 'Visible', 'off');
    for p = 1:numel(cfg.primary)
        m = cfg.primary{p};
        subplot(2, 4, p); hold on
        sc = 1; if ismember(m, cfg.msMetrics), sc = 1000; end
        mm = meanT(strcmp(meanT.Metric, m), :);
        for q = 1:3
            y = T.(['delta_' m])(strcmp(T.Intervention, ints{q})) * sc;
            y = y(~isnan(y));
            scatter(q + (rand(size(y)) - 0.5) * 0.35, y, 14, cols(q,:), 'filled', ...
                'MarkerFaceAlpha', 0.35);
            errorbar(q, mm.AdjMeanChange(q), mm.AdjMeanChange(q) - mm.CI_low(q), ...
                mm.CI_high(q) - mm.AdjMeanChange(q), 'o', 'Color', cols(q,:) * 0.7, ...
                'MarkerFaceColor', cols(q,:), 'LineWidth', 1.6, 'MarkerSize', 7);
        end
        plot([0.4 3.6], [0 0], 'k:');
        o = omniT(strcmp(omniT.Metric, m), :);
        lab = strrep(m, '_', ' ');
        if ~isempty(mm.Unit{1}), lab = [lab ' (' mm.Unit{1} ')']; end
        title(sprintf('\\Delta %s\np = %.3f (FDR %.3f)', lab, o.p, o.p_FDR), 'FontWeight', 'normal');
        set(gca, 'XTick', 1:3, 'XTickLabel', names, 'Box', 'off');
        xlim([0.4 3.6]);
    end
    print(fig, fullfile(cfg.outDir, 'SART_LMM_adjusted_change.png'), '-dpng', '-r150');
    close(fig);
end
fprintf('\nDone. Statistics saved in:\n%s\n', cfg.outDir);


%% =========================================================================
%%                              Local functions
%% =========================================================================
function D = prepData(T, y, pre)
% Keeps rows with a valid outcome (and pre value), builds dummies and
% centred covariates within the analysed subset.
ok = ~isnan(y);
if ~isempty(pre), ok = ok & ~isnan(pre); end
D = T(ok, :);
D.y = y(ok);
iv = string(D.Intervention);
D.E = double(iv == "E"); D.N = double(iv == "N"); D.R = double(iv == "R");
D.Position  = D.InterventionPosition - mean(D.InterventionPosition);
D.OrderPre  = double(string(D.TaskOrderPre)  == "SN"); D.OrderPre  = D.OrderPre  - mean(D.OrderPre);
D.OrderPost = double(string(D.TaskOrderPost) == "SN"); D.OrderPost = D.OrderPost - mean(D.OrderPost);
if ~isempty(pre)
    D.Pre = pre(ok) - mean(pre(ok));
end
end

function [lme, R] = fitAndTest(D, ancova, ints)
if ancova
    f = 'y ~ -1 + E + N + R + Pre + Position + OrderPre + OrderPost + (1|SubjectCat)';
else
    f = 'y ~ -1 + E + N + R + Position + OrderPre + OrderPost + (1|SubjectCat)';
end
lme  = fitlme(D, f, 'FitMethod', 'REML');
b    = fixedEffects(lme);
C    = lme.CoefficientCovariance;
nm   = lme.CoefficientNames;
nC   = numel(b);
idx  = cellfun(@(s) find(strcmp(nm, s)), ints);

% omnibus E = N = R
H = zeros(2, nC); H(1, idx([1 2])) = [1 -1]; H(2, idx([2 3])) = [1 -1];
[R.p, R.F, R.df1, R.df2] = coefTest(lme, H, [0; 0], 'DFMethod', 'satterthwaite');

% adjusted means and test vs 0
for q = 1:3
    h = zeros(1, nC); h(idx(q)) = 1;
    [R.meanP(q), ~, ~, R.meanDF(q)] = coefTest(lme, h, 0, 'DFMethod', 'satterthwaite');
    R.mean(q)   = b(idx(q));
    R.meanSE(q) = sqrt(C(idx(q), idx(q)));
    tc = tinv(0.975, R.meanDF(q));
    R.meanCI(q, :) = R.mean(q) + [-1 1] * tc * R.meanSE(q);
end

% pairwise
pairs = [1 2; 1 3; 2 3];
for q = 1:3
    h = zeros(1, nC); h(idx(pairs(q,1))) = 1; h(idx(pairs(q,2))) = -1;
    [R.pairP(q), ~, ~, R.pairDF(q)] = coefTest(lme, h, 0, 'DFMethod', 'satterthwaite');
    R.pairEst(q) = h * b;
    R.pairSE(q)  = sqrt(h * C * h');
    tc = tinv(0.975, R.pairDF(q));
    R.pairCI(q, :) = R.pairEst(q) + [-1 1] * tc * R.pairSE(q);
    R.pairName{q} = sprintf('%s-%s', ints{pairs(q,1)}, ints{pairs(q,2)});
end

% covariates
covIdx = setdiff(1:nC, idx);
R.covName = nm(covIdx);
for q = 1:numel(covIdx)
    h = zeros(1, nC); h(covIdx(q)) = 1;
    [R.covP(q), ~, ~, R.covDF(q)] = coefTest(lme, h, 0, 'DFMethod', 'satterthwaite');
    R.covEst(q) = b(covIdx(q));
    R.covSE(q)  = sqrt(C(covIdx(q), covIdx(q)));
end

% total SD (random intercept + residual) for standardised effect sizes
[psi, mse] = covarianceParameters(lme);
R.sdTotal = sqrt(psi{1}(1) + mse);
end

function q = bh(p)
% Benjamini-Hochberg FDR
p = p(:); q = nan(size(p));
v = find(~isnan(p));
[ps, o] = sort(p(v)); m = numel(ps);
qs = ps .* m ./ (1:m)';
qs = min(1, cummin(qs, 'reverse'));
q(v(o)) = qs;
end

function q = holm(p)
% Holm step-down correction
p = p(:); q = nan(size(p));
v = find(~isnan(p));
[ps, o] = sort(p(v)); m = numel(ps);
qs = (m - (1:m)' + 1) .* ps;
qs = min(1, cummax(qs));
q(v(o)) = qs;
end
