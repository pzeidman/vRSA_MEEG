function [PEB,F] = spm_eeg_rsa_peb(RSAs,options)
%==========================================================================
% spm_eeg_rsa_peb.m
%==========================================================================
% This function takes a set of RSA (Representational Similarity Analysis)
% models, runs a PEB (Parametric Empirical Bayes) analysis on them,
% compares reduced vs. full models for each parameter, and accumulates
% evidence across basis functions and components.
%
% FORMAT:
%   [PEB, F] = spm_eeg_rsa_peb(RSAs, options)
%
% INPUTS:
%   RSAs    - Cell array of RSA models (size: [subjects x basisFunctions])
%             Each cell is an RSA struct with fields like .Ep, .Cp, .M, etc.
%
%   options - Struct containing:
%       .params    : (default='all') Parameter(s) to estimate in PEB
%       .FIR_bf    : (logical, default=false) If true, indicates an FIR basis
%       .t         : Time vector (optional)
%       .data      : Actual data for plotting/comparison (optional)
%
% OUTPUTS:
%   PEB - The resulting PEB models for each RSA
%   F   - Matrix of log-evidence comparisons for each PEB model
%
% STEPS:
%   1) Fit a PEB model to each RSA, specifying which parameters to include.
%   2) Compare each model with and without each parameter -> log BF
%   3) Summarize evidence across basis functions & contrasts
%   4) Plot relevant results, if requested
%
% REFERENCES:
%   - spm_dcm_peb, spm_softmax, spm_log_evidence_reduce
%   - RSA references in neuroimaging
%
%==========================================================================

%--------------------------------------------------------------------------
% Parse inputs & defaults
%--------------------------------------------------------------------------
arguments
    RSAs cell
    options.params = 'all';
    options.FIR_bf (1,1) logical = false
    options.t (:, 1) {mustBeNumeric} = []
    options.doplot (1, 1) logical = true
end

% Handle empty input
if isempty(options.t)
    options.t = 1:RSAs{1, 1}.M.ntimes;
end
if ~isfield(RSAs{1, 1}, 'cnames') | isempty(RSAs{1, 1}.cnames)
    ncomp = length(RSAs{1, 1}.M.Q);
    RSAs{1, 1}.cnames = arrayfun(@(x) sprintf('cond-%d', x), 1:ncomp, ...
        'UniformOutput', false);
end

%--------------------------------------------------------------------------
% (1) Run a PEB model for each RSA
%--------------------------------------------------------------------------
M    = struct();
PEB  = spm_dcm_peb(RSAs, M, options.params);

%--------------------------------------------------------------------------
% (2) Compare PEB models with / without each parameter
%--------------------------------------------------------------------------
F = [];
for i = 1:length(PEB)
    % conditions x basis functions
    F(:,i) = compare_models(PEB(i));
end

%--------------------------------------------------------------------------
% (3) Accumulate evidence over basis functions & contrasts
%--------------------------------------------------------------------------
F_bfs       = sum(F,1);   % evidence per BF (summed across conditions)
F_contrasts = sum(F,2);   % evidence per contrast/condition (summed across BFs)

Pp_bfs = spm_softmax([F_bfs(:)'; zeros(1,length(F_bfs))]);
Pp_bfs = Pp_bfs(1,:);

Pp_contrasts = spm_softmax([F_contrasts(:)'; zeros(1,length(F_contrasts))]);
Pp_contrasts = Pp_contrasts(1,:);

%--------------------------------------------------------------------------
% (4) Gather second-order parameters (BB) per basis function
%--------------------------------------------------------------------------
BB_per_bf = {};
for i = 1:size(RSAs,1)
    for j = 1:size(RSAs,2)
        if i == 1
            BB_per_bf{j} = RSAs{i,j}.BB;
        else
            BB_per_bf{j} = RSAs{i,j}.BB + RSAs{i,j}.BB;
        end
    end
end

%--------------------------------------------------------------------------
% (5) Observed and modelled 2nd order for "winning" or included BFs
%--------------------------------------------------------------------------
included_bfs = find(Pp_bfs > 0.9);
BB = [];
G  = [];
for i = 1:size(RSAs,1)
    for j = included_bfs
        if i == 1
            BB = RSAs{i,j}.BB;
            G  = RSAs{i,j}.G;
        else
            BB = BB + RSAs{i,j}.BB;
            G  = G + RSAs{i,j}.G;
        end
    end
end

%--------------------------------------------------------------------------
% (6) Extract parameters posterior expectation (Ep) and variance (Cp)
%--------------------------------------------------------------------------
Ep = full(horzcat(PEB.Ep));
Cp = cell2mat(arrayfun(@(x) diag(full(x.Cp)), PEB, ...
    'UniformOutput', false));

%--------------------------------------------------------------------------
% (6) Stop if no plotting is requested
%--------------------------------------------------------------------------
if ~options.doplot
    return
end

%==========================================================================
% PLOTTING SECTION
%==========================================================================

%--------------------------------------------------------------------------
% Extract basic info
%--------------------------------------------------------------------------
nXt    = RSAs{1, 1}.M.nXt;    % # basis functions
Xt     = RSAs{1, 1}.M.Xt;     % First-order basis design
Xt_names = RSAs{1, 1}.M.Xt_names;     % First-order basis design
Qnames = RSAs{1, 1}.cnames;   % Names of components
nQ     = size(Qnames, 2);

%--------------------------------------------------------------------------
% FIGURE 1: First and second order design
%--------------------------------------------------------------------------
fig = figure('Name', 'First and second order design', 'Position', [10, 10, 725, 725], 'Color',[1 1 1]);

% Create a 3-row tiled layout (top=1/3, bottom=2/3)
mainLayout = tiledlayout(fig, 3, 1, 'TileSpacing','compact', 'Padding','compact');

%------------------------------------
% (1) TOP: Plot first-order basis
%------------------------------------
nexttile(mainLayout, 1);
if options.FIR_bf
    imagesc(1:nXt, options.t, Xt);
    axis square;
    xlabel('FIR');
    ylabel('Time (s)');
    title('Basis functions');
    set(gca, 'FontSize',12);
    colormap gray;
else
    plot(options.t, Xt);
    xlabel('Time (ms)');
    title('Basis functions');
    set(gca, 'FontSize',12);
    xlim([options.t(1), options.t(end)]);
end

%------------------------------------
% (2) MIDDLE+BOTTOM: Second-order matrices
%------------------------------------
nGrid = ceil(sqrt(nQ));  % dimension of grid for Q-matrices

% Nested layout: tile #2, spanning two rows
layout2 = tiledlayout(mainLayout, nGrid, nGrid, 'TileSpacing','compact', 'Padding','compact');
layout2.Layout.Tile     = 2;
layout2.Layout.TileSpan = [2 1];
title(layout2, 'Model matrices');

for q = 1:nQ
    nexttile(layout2);
    imagesc(RSAs{1}.M.Q{q});
    colormap gray;
    axis square;
    title(Qnames{q});
    set(gca, 'FontSize',12);
end

%--------------------------------------------------------------------------
% FIGURE 2: Group evidence per basis function and component
%--------------------------------------------------------------------------
if options.FIR_bf
    fig = figure('Name', 'Time resolved Posterior and evidence', 'Position', [10, 10, 725, 725], 'Color',[1 1 1]);
    mainLayout = tiledlayout(3, 1, 'TileSpacing','compact', 'Padding','compact');

    % Downsample time vector according to FIR design
    t = spm_pinv(Xt) * options.t;

    % Approximate bin duration
    bin_dur = options.t(sum(Xt(:, 2))) - options.t(1);

    % (A) Plot the estimated parameters:
    nexttile(mainLayout)
    ctr = 1;
    for q = 1:nQ-1
        hold on
        spm_plot_ci(Ep(q,:), Cp(q,:), t);
        lbl{ctr}   = '';
        lbl{ctr+1} = Qnames{q};
        ctr        = ctr + 2; 
    end
    title('Estimated hyperparameters');
    ylabel('v');
    xlabel('Time');
    xlim([t(1), t(end)]);
    legend(lbl, 'Location','best');
    set(gca, 'FontSize',12);
    
    % (B) Plot the exponentiated estimated parameters:
    nexttile(mainLayout)
    for q = 1:nQ-1
        hold on
        spm_plot_ci(Ep(q,:), Cp(q,:), t, [], 'exp');
    end
    title('Estimated hyperparameters');
    ylabel('Weights (exp)');
    xlabel('Time');
    xlim([t(1), t(end)]);
    legend(lbl, 'Location','best');
    set(gca, 'FontSize',12);

    % (C) Plot log evidence per basis function
    nexttile(mainLayout)
    h = plot(t, F, '-o');
    set(h, {'MarkerFaceColor'}, get(h,'Color')); 
    title('Log evidence per basis function');
    ylabel('log BF');
    xlabel('Time');
    xlim([t(1), t(end)]);
    legend(RSAs{1}.cnames, 'Location','best');
    set(gca, 'FontSize',12);

    % (B) Covariance matrices at each time point, etc.
    fig = figure('Name', 'Group evidence BF', 'Position', [10, 10, 725, 725], 'Color',[1 1 1]);
    nsubrows = ceil(nXt/2);
    mainLayout = tiledlayout(nsubrows, 4, 'TileSpacing','compact','Padding','compact');
    title(mainLayout, 'Covariance matrices');
    ax_logbf = [];
    ctr = 1;
    for i = 1:(size(RSAs,2))
        nexttile(mainLayout)
        imagesc(BB_per_bf{i})
        colormap gray;
        set(gca,'YTickLabel',[], 'XTickLabel',[]);
        title(sprintf("%.2f-%.2f s", t(i) - bin_dur / 2, t(i) + bin_dur / 2))
        axis square;

        ax_logbf(ctr) = nexttile(mainLayout);
        bar(F(:, i), 'FaceColor', [0.3, 0.3, 0.3])
        hold on
        yline(3, 'Color', 'r', 'LineStyle', '--')
        set(gca,'XTickLabel',[]);
        ctr = ctr + 1;
    end
    set(gca,'XTickLabel',Qnames);
    ylabel('LogBF')
    set(gca,'FontSize',9);
    linkaxes(ax_logbf)

else
    fig = figure('Name', 'Group evidence BF', 'Position', [10, 10, 725, 725], 'Color',[1 1 1]);
    tiledlayout(size(RSAs,2), 2, 'TileSpacing','compact','Padding','compact');
    for i = 1:(size(RSAs,2))
        nexttile()
        imagesc(BB_per_bf{i})
        colormap gray;
        axis square;
        xlabel("Stimuli")
        ylabel("Stimuli")
        set(gca,'FontSize',9);
        title(sprintf("Beta 2nd order (%s)", Xt_names{i}), 'FontSize', 12)

        nexttile()
        bar(F(:, i), 'FaceColor', [0.3, 0.3, 0.3])
        axis square;
        ylabel('Log evidence (nats)')
        set(gca,'FontSize',9);
        set(gca,'XTickLabel',Qnames);
        title(sprintf("evidence %s", Xt_names{i}), 'FontSize', 12)
    end
end

%--------------------------------------------------------------------------
% FIGURE 3: Pooled evidence and second-order parameters
%--------------------------------------------------------------------------
fig = figure('Name', 'Group vRSA analysis', 'Position', [10, 10, 725, 725], 'Color',[1 1 1]);
mainLayout = tiledlayout(3, 2, 'TileSpacing','compact','Padding','compact');

% (A) Pooled evidence per component
nexttile(mainLayout)
bar(F_contrasts, 'FaceColor', [0.3, 0.3, 0.3]);
set(gca,'XTick',1:nQ,'XTickLabel',RSAs{1}.cnames);
ylabel('Log evidence (nats)');
set(gca,'FontSize',12);
title('Pooled evidence per component','FontSize',16);

nexttile(mainLayout)
bar(Pp_contrasts, 'FaceColor', [0.3, 0.3, 0.3]);
set(gca,'XTick',1:nQ,'XTickLabel',RSAs{1}.cnames);
ylabel('Probability');
set(gca,'FontSize',12);
title('Probability per component','FontSize',16);

% (B) Evidence per basis function
if options.FIR_bf
    subLayout = tiledlayout(mainLayout, 1, 1, 'TileSpacing','compact','Padding','compact');
    subLayout.Layout.Tile     = 3;
    subLayout.Layout.TileSpan = [1 2];
    nexttile(subLayout)

    yyaxis right
    hold on
    b = bar(t, Pp_bfs, 'FaceColor', [0.3, 0.3, 0.3]);
    b.FaceAlpha = 0.2;
    b.EdgeColor = 'none';
    ylabel('Probability');

    yyaxis left
    plot(t, F_bfs, 'LineWidth', 2);
    ylabel('log BF');
    xlabel('Time (sec)');
    set(gca,'FontSize',12);
    title('Log evidence over time','FontSize',16)
    ax = gca;
    ax.YAxis(2).Color = [0.3, 0.3, 0.3];
else
    nexttile(mainLayout)
    bar(F_bfs, 'FaceColor', [0.3, 0.3, 0.3]);
    ylabel('Log evidence (nats)');
    xlabel('Basis function');
    xticklabels(Xt_names)
    set(gca,'FontSize',12);
    title('Pooled evidence per BF','FontSize',16);

    nexttile(mainLayout)
    bar(Pp_bfs, 'FaceColor', [0.3, 0.3, 0.3]);
    ylabel('Probability');
    xticklabels(Xt_names)
    xlabel('Basis function');
    set(gca,'FontSize',12);
    title('Probability per BF','FontSize',16);
end

% (C) Visualize second-order Betas vs. vRSA
nexttile(mainLayout)
imagesc(BB);
colormap gray;
xlabel('Condition'); ylabel('Condition');
set(gca,'FontSize',12);
title('Betas (second order)','FontSize',16);
axis square;

nexttile(mainLayout)
imagesc(G);
xlabel('Condition'); ylabel('Condition');
set(gca,'FontSize',12);
title('vRSA','FontSize',16);
axis square;

%==========================================================================
% SUBFUNCTION: compare_models
%==========================================================================
function [F,P] = compare_models(PEB)
% Compare each parameter switched off vs. on in a PEB model
% and produce log BF values for each parameter.

n = length(PEB.Pnames);

F = zeros(1,n);
for i = 1:n
    F(i) = reduce_parameters(PEB,i);
end

% Flip log BF to evidence in favour of the full model
F = -F;

% Compute posterior probability that each parameter is needed
P = spm_softmax([F; zeros(1,n)]);
P = P(1,:);

%==========================================================================
% SUBFUNCTION: reduce_parameters
%==========================================================================
function logBF = reduce_parameters(PEB,q)
% Switch off the q-th parameter in the PEB model, then compute the log Bayes
% factor in favour of that reduced model (versus the full one).

Ep = PEB.Ep;
Cp = PEB.Cp;
pE = PEB.M.pE;
pC = PEB.M.pC;

% Remove prior variance from parameter q
rC = pC;
rC(q,q) = 0;

% Evaluate difference in log evidence for the reduced model
logBF = spm_log_evidence_reduce(Ep, Cp, pE, pC, pE, rC);