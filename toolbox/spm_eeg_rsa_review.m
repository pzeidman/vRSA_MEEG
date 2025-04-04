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
%   RSA - Structure containing model-fitting results from spm_eeg_rsa_estimate. 
%       The structure must contain the following fields:
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
%       .t            Time vector (default = 1:RSA.M.ntimes)
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
    RSA struct
    options.FIR_bf (1,1) logical = false
    options.t (:, 1) {mustBeNumeric} = []
    options.data (:, :, :) = []
end

if isempty(options.t)
    % Default time vector covers all time bins in the model
    options.t = (1:RSA.M.ntimes)';
end
if ~isfield(RSA, 'cnames') | isempty(RSA.cnames)
    ncon = length(RSA.con);
    RSA.cnames = arrayfun(@(x) sprintf('contrast-%d', x), 1:ncon, ...
        'UniformOutput', false);
end
if ~isempty(options.data)
    for m = 1:RSA.M.nmodes
        for c = 1:RSA.M.nconditions
            options.data(m,:,c) = options.data(m,:,c) - ...
                mean(options.data(m,:,c));
        end
    end
end
%--------------------------------------------------------------------------
% Unpack info from RSA structure
%--------------------------------------------------------------------------
nmodes      = RSA.M.nmodes;                       % Number of modes/channels
nconditions = RSA.M.nconditions;                  % Number of conditions
B           = RSA.B;                              % First order betas
nXt         = RSA.M.nXt;                          % # basis functions
Xt          = RSA.M.Xt;                           % First-order basis design
con         = RSA.con;                            % Contrasts
ncon        = length(con);                        % Number of contrasts
nQ          = length(RSA.M.Q);                    % Number of components (nXt * ncon)
qcon        = RSA.qcon;                           % Index of each contrast
qbf         = RSA.qbf;                             % Index of each bf
Ep          = [RSA.Ep.cond, RSA.Ep.noise];        % Posterior mean
Vp          = RSA.Cp;                             % Posterior covariance
logBF       = reshape(RSA.logBF.cond, ncon, nXt); % (contrasts x basis functions)
cnames      = RSA.cnames;                         % Names of components
qnames      = RSA.qnames;                         % Names of components
Xt_names    = RSA.M.Xt_names;                     % Names of the basis functions
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
nGrid = ceil(sqrt(ncon));  % dimension of grid for Q-matrices

% Create nested layout in tile #2, spanning 2 rows (tiles #2 & #3)
layout2 = tiledlayout(mainLayout, nGrid, nGrid, ...
    'TileSpacing','compact', 'Padding','compact');
layout2.Layout.Tile     = 2; 
layout2.Layout.TileSpan = [2 1];
title(layout2, 'Model matrices');

for c = 1:ncon
    nexttile(layout2);
    imagesc(con{c} * con{c}');
    colormap gray;
    axis square;
    title(cnames{c});
    set(gca, 'FontSize',12);
end

%--------------------------------------------------------------------------
% FIGURE 2: First order fit:
%--------------------------------------------------------------------------
fig = figure('Name', 'First order fit', 'Position', [10, 10, 725, 725], 'Color',[1 1 1]);

% Layout for various subplots
mainLayout = tiledlayout(fig, 3, 3, ...
    'TileSpacing','compact', 'Padding','compact');

% Draw 9 random channels if more than 9 channels:
if nmodes > 8
    ch_ind = randsample(nmodes, 9, 1);
else
    ch_ind = 1:nmodes;
end
% Recreate the first order design matrix:
X = kron(eye(nconditions), Xt);

for i = 1:length(ch_ind)
    nexttile(mainLayout)
    % Average the beta across trials for each basis function
    Bi = B(:, ch_ind(i));
    % Reshape betas:
    Bi = reshape(Bi.', [], 1);

    % Calculated fitted:
    fitted = X * Bi;
    
    % Compute average across conditions/trials
    fitted2D = reshape(fitted, size(X, 1)/nconditions, ...
        nconditions);
    plot(options.t, fitted2D, 'Color', [0, 0, 0, 0.2], 'LineWidth', 0.5)
    hold on
    if ~isempty(options.data)
        % Mean correct each ERP over channel
        Y = options.data - mean(options.data, 1);
        Y = squeeze(Y(ch_ind(i), :, :));
        plot(options.t, Y, 'Color', [1, 0, 0, 0.2], 'LineWidth', 0.5)
    end
    xlim([options.t(1), options.t(end)]);
    xline(0);
    title(sprintf('Channel %d',ch_ind(i)));
    set(gca,'YTickLabel',[]);
end

%--------------------------------------------------------------------------
% FIGURE 3: Estimated parameter in log space
%--------------------------------------------------------------------------
if options.FIR_bf
    %----------------------------------------------------------------------
    % FIR/time-based model
    %----------------------------------------------------------------------
    fig = figure('Name', 'Estimated parameters', 'Position', [10, 10, 725, 725], 'Color',[1 1 1]);
    mainLayout = tiledlayout(fig, 3, 1, 'TileSpacing','compact', 'Padding','compact');
    
    t = spm_pinv(Xt) * options.t;  % Downsample time to match the FIR
    % (1) Posterior (Ep) over time
    nexttile(mainLayout, 1);
    lbl = cell(ncon*2,1);
    ctr = 1;
    for c = 1:ncon
        hold on
        spm_plot_ci(Ep(:,qcon == c), diag(Vp(qcon == c,qcon == c)), t);
        lbl{ctr}   = '';
        lbl{ctr+1} = cnames{c};
        ctr        = ctr + 2;
    end
    title('Posterior over time');
    set(gca, 'FontSize',12);
    ylabel('v');
    xlabel('Time');
    xlim([t(1), t(end)]);
    legend(lbl);

    % (2) Exponential of posterior (weights)
    nexttile(mainLayout, 2);
    hold on
    for c = 1:ncon
        hold on
        spm_plot_ci(Ep(:,qcon == c), diag(Vp(qcon == c,qcon == c)), t, [], 'exp');
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
    fig = figure('Name', 'Estimated parameters', 'Position', [10, 10, 725, 725], 'Color',[1 1 1]);
    nGrid = ceil(sqrt(nXt));  % dimension of grid for Q-matrices
    mainLayout = tiledlayout(fig, nGrid, nGrid, 'TileSpacing','compact', 'Padding','compact');
    
    for bf = 1:nXt
        ax = nexttile(mainLayout);
        spm_plot_ci(Ep(:, qbf == bf)', diag(Vp(qbf == bf, qbf == bf)));
        axis square;
        title(Xt_names{bf});
        set(ax, 'FontSize',12);
        if bf == nXt
            ylabel('Posterior');
            xticklabels(cnames);
        end
    end

    % (2) Exponential of posterior
    fig = figure('Name', 'Estimated weights', 'Position', [10, 10, 725, 725], 'Color',[1 1 1]);
    nGrid = ceil(sqrt(nXt));  % dimension of grid for Q-matrices
    mainLayout = tiledlayout(fig, nGrid, nGrid, 'TileSpacing','compact', 'Padding','compact');
    for bf = 1:nXt
        ax = nexttile(mainLayout);
        spm_plot_ci(Ep(:, qbf == bf)', diag(Vp(qbf == bf, qbf == bf)), [], [], 'exp');
        axis square;
        title(Xt_names{bf});
        set(ax, 'FontSize',12);
        if bf == nXt
            ylabel('Posterior');
            xticklabels(cnames);
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
BB = zeros(nconditions);
for bf = 1:nXt
    BB = BB + RSA.BB(bf:nXt:end, bf:nXt:end);
end
subplot(rows, cols, 1:2);
imagesc(BB);
colormap gray;
axis square;
title('Betas (2nd order)', 'FontSize',16);

% (2) Summation of vRSA (G)
G = zeros(nconditions);
for bf = 1:length(RSA)
    G = G + RSA.G(bf:nXt:end, bf:nXt:end);
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
    legend(cnames, 'Location','best');
    set(gca, 'FontSize',12);
else
    % For non-FIR, logBF is a matrix (ncon x nXt)
    if ncon == 1
        bar(logBF, 'FaceColor','k');
        str = sprintf('Total: %2.2f nats', sum(logBF));
        title({'Log evidence', str}, 'FontSize',16);
        ylabel('F (with vs without)');
    else
        imagesc(logBF);
        colormap gray;
        set(gca,'YTick',1:nQ, 'YTickLabel', cnames);
        set(gca,'XTick',1:nXt, 'XTickLabel',1:nXt);
        colorbar;
        title({'Log evidence','per condition & BF'}, 'FontSize',16);
    end
    xlabel('Basis function');
    set(gca, 'FontSize',12);
    axis square;
end

if ncon > 1
    % (4) Log evidence per condition
    F_cond = sum(logBF,2);
    subplot(rows, cols, 9:10);
    bar(F_cond, 'FaceColor','k');
    set(gca, 'XTick',1:ncon, 'XTickLabel', cnames);
    colormap gray;
    set(gca, 'FontSize',12);
    title({'Log evidence','per condition'}, 'FontSize',16);
    axis square;

    % (5) Log evidence per basis function
    F_bf = sum(logBF,1);
    subplot(rows, cols, 11:12);
    bar(F_bf, 'FaceColor','k');
    set(gca, 'XTick',1:nXt, 'XTickLabel',1:nXt);
    xlabel('Basis function');
    colormap gray;
    set(gca, 'FontSize',12);
    title({'Log evidence','per BF'}, 'FontSize',16);
    axis square;
end

end % function spm_eeg_rsa_review
