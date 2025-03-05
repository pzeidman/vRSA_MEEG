function spm_eeg_rsa_review(RSA, options)
%==========================================================================
% spm_eeg_rsa_review.m
%==========================================================================
% Review vRSA/MEEG results (for a single subject). This function plots the 
% first and second order models (i.e. design matrix and model trial by
% trial covariance matrices), the fitted first order responses alongside
% the posterior distribution of the weights for each model covariance
% matrices, and finally, the evidence for each model matrices and basis
% function. 
% If the first order responses were modelled using an FIR filter, the
% posteriors and evidence are plotted as time series, otherwise as barplots
% and heatmaps
%
% FORMAT:
%   spm_eeg_rsa_review(RSA, options)
%
% Inputs:
%   RSA - Cell array containing model-fitting results. Each cell is
%         a structure with fields:
%           .Ep      Posterior estimates (struct or numeric)
%           .Cp      Posterior covariance
%           .logBF   Log Bayesian Factor evidence
%           .M       Model structure with basis functions, etc.
%           .cnames  Names of model components/conditions
%           .BB      Second-order betas
%           .G       Summed vRSA
%
%   options - Struct with optional fields:
%       .FIR_bf  (default = false)
%                     If true, treat design as FIR/time-based, and plotting
%                     the results accordingly
%       .t            Time vector (default = 1:RSA{1}.M.ntimes)
%       .data         Data for comparison against model fits (optional)
%
% This function:
%   1) Plots the first-level (basis functions, second-order model matrices)
%   2) Plots posteriors for each basis function and component
%   3) Summarizes evidence (logBF) by condition and basis function
%==========================================================================

%--------------------------------------------------------------------------
% Parse input arguments & defaults
%--------------------------------------------------------------------------
arguments
    RSA cell
    options.FIR_bf (1,1) logical = false
    options.t (:, 1) {mustBeNumeric} = []
    options.data (:, :, :) = []
end

if isempty(options.t)
    % Default time vector covers all time bins in the model
    options.t = (1:RSA{1}.M.ntimes)';
end
if ~isfield(RSA{1, 1}, 'cnames') | isempty(RSA{1, 1}.cnames)
    ncomp = length(RSA{1, 1}.M.Q);
    RSA{1, 1}.cnames = arrayfun(@(x) sprintf('cond-%d', x), 1:ncomp, ...
        'UniformOutput', false);
end
if ~isempty(options.data)
    for m = 1:RSA{1}.M.nmodes
        for c = 1:RSA{1}.M.nconditions
            options.data(m,:,c) = options.data(m,:,c) - ...
                mean(options.data(m,:,c));
        end
    end
end
%--------------------------------------------------------------------------
% Basic model info & extraction
%--------------------------------------------------------------------------
nXt        = RSA{1}.M.nXt;                     % # basis functions
Xt         = RSA{1}.M.Xt;                      % First-order basis design
Ep         = cellfun(@(x) spm_vec(x.Ep),  RSA, 'UniformOutput', false);
Vp         = cellfun(@(x) diag(x.Cp),     RSA, 'UniformOutput', false);
logBF      = cellfun(@(x) spm_vec(x.logBF.cond), RSA, 'UniformOutput', false);
logBF      = cell2mat(logBF);                  % (components x basis functions)
Qnames     = RSA{1}.cnames;                    % Names of components
Xt_names   = RSA{1}.M.Xt_names;                % Names of the basis functions
nQ         = size(Qnames, 2);

%--------------------------------------------------------------------------
% FIGURE 1: First and second order design
%--------------------------------------------------------------------------
fig = figure('Name', 'First and second order design', 'Position', [10, 10, 725, 725], 'Color',[1 1 1]);
% Split vertically into 3 rows (top=1/3, bottom=2/3)
mainLayout = tiledlayout(fig, 3, 1, ...
    'TileSpacing','compact', 'Padding','compact');

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
% (2) MIDDLE+BOTTOM: Plot Q-matrices (Second-order)
%------------------------------------
nGrid = ceil(sqrt(nQ));  % dimension of grid for Q-matrices

% Create nested layout in tile #2, spanning 2 rows (tiles #2 & #3)
layout2 = tiledlayout(mainLayout, nGrid, nGrid, ...
    'TileSpacing','compact', 'Padding','compact');
layout2.Layout.Tile     = 2; 
layout2.Layout.TileSpan = [2 1];
title(layout2, 'Model matrices');

for q = 1:nQ
    nexttile(layout2);
    imagesc(RSA{1}.M.Q{q});
    colormap gray;
    axis square;
    title(Qnames{q});
    set(gca, 'FontSize',12);
end

%--------------------------------------------------------------------------
% FIGURE 2: Estimated parameters & fits
%--------------------------------------------------------------------------
fig = figure('Name', 'Estimated parameters', 'Position', [10, 10, 725, 725], 'Color',[1 1 1]);

% Layout for various subplots
mainLayout = tiledlayout(fig, 3, 1, ...
    'TileSpacing','compact', 'Padding','compact');
%------------------------------------
% (A) Plot first-order fit (for demonstration)
%------------------------------------
layout1 = tiledlayout(mainLayout,3,3, ...
    'TileSpacing','compact', 'Padding','compact');
layout1.Layout.Tile = 1;
% Draw 9 random channels:
if RSA{1}.M.nmodes > 8
    ch_ind = randsample([RSA{1}.M.nmodes], 9, 1);
else
    ch_ind = 1:RSA{1}.M.nmodes;
end
% Recreate the first order design matrix:
X = kron(eye(RSA{1}.M.nconditions), Xt);

for i = 1:length(ch_ind)
    nexttile(layout1)
    % Average the beta across trials for each basis function
    B = cell2mat(cellfun(@(x) x.B(:, ch_ind(i)),  RSA, ...
        'UniformOutput', false));
    % Reshape betas:
    B = reshape(B.', [], 1);

    % Calculated fitted:
    fitted = X * B;
    
    % Compute average across conditions/trials
    fitted2D = reshape(fitted, size(X, 1)/RSA{1}.M.nconditions, ...
        RSA{1}.M.nconditions);
    plot(options.t, fitted2D, 'Color', [0, 0, 0, 0.2], 'LineWidth', 0.5)
    hold on
    if ~isempty(options.data)
        % Mean correct each ERP over time
        Y = options.data - mean(options.data, 1);
        Y = squeeze(Y(ch_ind(i), :, :));
        for c = 1:size(Y, 2)
            Y(:,c) = Y(:,c) - mean(Y(:,:), 'all');
        end
        % Mean correct over channel:
        Y = Y - mean(Y);
        plot(options.t, Y, 'Color', [1, 0, 0, 0.2], 'LineWidth', 0.5)
    end
    xlim([options.t(1), options.t(end)]);
    xline(0);
    title(sprintf('Channel %d',ch_ind(i)));
    set(gca,'YTickLabel',[]);
end

%------------------------------------
% (B) Posterior parameters
%------------------------------------
if options.FIR_bf
    %----------------------------------------------------------------------
    % FIR/time-based model
    %----------------------------------------------------------------------
    t = spm_pinv(Xt) * options.t;  % Downsample time to match the FIR

    Ep_mat = cell2mat(Ep);
    Vp_mat = cell2mat(Vp);

    % (1) Posterior (Ep) over time
    nexttile(mainLayout, 2);
    lbl = cell(nQ*2,1);
    ctr = 1;
    for q = 1:nQ
        hold on
        spm_plot_ci(Ep_mat(q,:), Vp_mat(q,:), t);
        lbl{ctr}   = '';
        lbl{ctr+1} = Qnames{q};
        ctr        = ctr + 2;
    end
    title('Posterior over time');
    set(gca, 'FontSize',12);
    ylabel('v');
    xlabel('Time');
    xlim([t(1), t(end)]);
    legend(lbl);

    % (2) Exponential of posterior (weights)
    nexttile(mainLayout, 3);
    hold on
    for q = 1:nQ
        hold on
        spm_plot_ci(Ep_mat(q,:), Vp_mat(q,:), t, [], 'exp');
    end
    title('Weight over time');
    set(gca, 'FontSize',12);
    ylabel('$\lambda$', 'interpreter','latex');
    xlabel('Time');
    xlim([t(1), t(end)]);
    legend(lbl);

else
    %----------------------------------------------------------------------
    % Non-FIR basis
    %----------------------------------------------------------------------
    % (1) Posterior per component
    layout2 = tiledlayout(mainLayout,1,nXt, ...
        'TileSpacing','compact', 'Padding','compact');
    layout2.Layout.Tile = 2;
    for bf = 1:nXt
        ax = nexttile(layout2);
        spm_plot_ci(Ep{bf}, diag(Vp{bf}));
        axis square;
        title(Xt_names{bf});
        set(ax, 'FontSize',12);
        xticklabels(Qnames);
        if bf == 1
            ylabel('Posterior');
        end
    end

    % (2) Exponential of posterior
    layout3 = tiledlayout(mainLayout,1,nXt);
    layout3.Layout.Tile = 3; 
    for bf = 1:nXt
        ax = nexttile(layout3);
        spm_plot_ci(Ep{bf}, diag(Vp{bf}), [], [], 'exp');
        axis square;
        title(Xt_names{bf});
        set(ax, 'FontSize',12);
        xticklabels(Qnames);
        if bf == 1
            ylabel('Weight');
        end
    end
end

%--------------------------------------------------------------------------
% FIGURE 3: Summaries (2nd-order parameters & log evidences)
%--------------------------------------------------------------------------
fig = figure('Name', 'vRSA - matrices', 'Position', [10, 10, 725, 725], 'Color',[1 1 1]);

rows = 3; 
cols = 4;

% (1) Sum second-order betas (BB) across all BFs
BB = RSA{1}.BB;
for bf = 2:length(RSA)
    BB = BB + RSA{bf}.BB;
end
subplot(rows, cols, 1:2);
imagesc(BB);
colormap gray;
axis square;
title('Betas (2nd order)', 'FontSize',16);

% (2) Summation of vRSA (G)
G = RSA{1}.G;
for bf = 2:length(RSA)
    G = G + RSA{bf}.G;
end
subplot(rows, cols, 3:4);
imagesc(G);
axis square;
title('vRSA', 'FontSize',16);
colormap gray;

% (3) Log evidence per condition & basis function
subplot(rows, cols, 6:7);
if options.FIR_bf
    % For FIR, logBF is time-varying
    t = spm_pinv(Xt) * options.t;
    plot(t, logBF);
    title('Log evidence per basis function');
    ylabel('log BF');
    xlabel('Time');
    xlim([t(1), t(end)]);
    legend(RSA{1}.cnames, 'Location','best');
    set(gca, 'FontSize',12);
else
    % For non-FIR, logBF is a matrix (nQ x nXt)
    if nQ == 1
        bar(logBF, 'FaceColor','k');
        str = sprintf('Total: %2.2f nats', sum(logBF));
        title({'Log evidence', str}, 'FontSize',16);
        ylabel('F (with vs without)');
    else
        imagesc(logBF);
        colormap gray;
        set(gca,'YTick',1:nQ, 'YTickLabel',RSA{1}.cnames);
        set(gca,'XTick',1:nXt, 'XTickLabel',Xt_names);
        colorbar;
        title({'Log evidence','per condition & BF'}, 'FontSize',16);
    end
    xlabel('Basis function');
    set(gca, 'FontSize',12);
    axis square;
end

if nQ > 1
    % (4) Log evidence per condition
    F_cond = sum(logBF,2);
    subplot(rows, cols, 9:10);
    bar(F_cond, 'FaceColor','k');
    set(gca, 'XTick',1:nQ, 'XTickLabel',RSA{1}.cnames(1:nQ));
    colormap gray;
    set(gca, 'FontSize',12);
    title({'Log evidence','per condition'}, 'FontSize',16);
    axis square;

    % (5) Log evidence per basis function
    F_bf = sum(logBF,1);
    subplot(rows, cols, 11:12);
    bar(F_bf, 'FaceColor','k');
    set(gca, 'XTick',1:nQ, 'XTickLabel',Xt_names);
    xlabel('Basis function');
    colormap gray;
    set(gca, 'FontSize',12);
    title({'Log evidence','per BF'}, 'FontSize',16);
    axis square;
end

end % function spm_eeg_rsa_review
