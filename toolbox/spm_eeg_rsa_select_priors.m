function [pE, pV] = spm_eeg_rsa_select_priors(c, bf, s, nmodes, nconditions, nsub,  options)
% Grid search for optimal priors in RSA-based EEG analysis.
%
% This function performs a grid search over prior expectations and variances 
% to identify the optimal prior parameters for a representational similarity 
% analysis (RSA) model applied to MEEG data. The goal is to find the prior 
% settings that maximize the sum Bayesian evidence in favor of effects
% known to be on and against effects known to be off
%
% The approach involves:
%   1) Simulating multivariate EEG data with known ground-truth contrasts.
%   2) Fitting RSA models across different prior expectations and variances.
%   3) Extracting Bayesian model evidence for "on" vs. "off" effects.
%   4) Selecting the prior parameters that optimize free energy differences.
%
% FORMAT:
%   [pV, pE] = spm_eeg_rsa_select_priors(D, c, bf, options)
%
% Inputs:
%   D  - EEG data matrix or cell array of matrices [modes x time x trials].
%   c  - Contrast matrix defining experimental conditions.
%   bf - Basis function matrix modeling expected neural responses.
%
% Options:
%   options.doplot  - Plot the evidence or not (default = []).
%
% Outputs:
%   pV - Optimal prior variance maximizing Bayesian evidence difference.
%   pE - Optimal prior expectation maximizing Bayesian evidence difference.
%
% This function:
%   1) Estimates noise variance directly from the empirical EEG data.
%   2) Simulates EEG responses using realistic noise levels and spatial 
%       covariance.
%   3) Fits RSA models across a predefined range of prior parameters.
%   4) Performs Bayesian Model Comparison (BMC) to assess model fit.
%   5) Selects the optimal prior expectation (pE) and variance (pV).
%   6) Optionally, visualizes the grid search results if plotting is 
%       enabled.
%==========================================================================

%--------------------------------------------------------------------------
% Parse input arguments & defaults
%--------------------------------------------------------------------------
arguments
    c                   % Contrast matrix (numeric or cell array)
    bf {mustBeNumeric}  % Basis function
    s {mustBeNumeric}   % Variance
    nmodes  {mustBeNumeric}   % Number of channels
    nconditions  {mustBeNumeric}   % Number of conditions to simulate
    nsub   {mustBeNumeric}   % Number of subjects to simulate
    options.doplot (1, 1) logical = true
end

%--------------------------------------------------------------------------
% Handle arguments;
%--------------------------------------------------------------------------
if iscell(c),    c = cell2mat(c); end

%--------------------------------------------------------------------------
% Define basic parameters
%--------------------------------------------------------------------------
nXt     = size(bf, 2);   % Number of basis functions
nC      = size(c, 2);    % Number of contrasts

% Define prior search space
pEs     = -16:2:-2;
pVs     = [2 4 8 16 32 64 128];

% Preallocate for on and off and face validity:
F_on          = zeros(length(pVs), length(pEs));
F_off         = zeros(length(pVs), length(pEs));
Eps_beta_corr = zeros(length(pVs), length(pEs), nsub);

%  Switch half of the effects on, the other half off:
CV = zeros(nXt, nC);
ind_on = randsample(nXt * nC, ceil(nXt * nC/2));
CV(ind_on) = 1; % Turn these effects on
%% 1) Simulate the data:
[Y, B_true] = spm_eeg_simulate_covariance(bf, c, s, nmodes, nsub, CV);

%% 2) Grid search for optimal priors:

% Loop through each Prior expectations:
for pE_i = 1:length(pEs)
    % Loop through each Prior variances:
    parfor pV_i = 1:length(pVs)
        S        = struct();
        S.Xt     = bf;
        S.con_c  = mat2cell(c, size(c, 1), ones(1, size(c, 2)));  % Between trial;
        S.pE = pEs(pE_i);
        S.pV = pVs(pV_i);

        %------------------------------------------------------------------
        % Fit RSA for each subject with current priors
        %------------------------------------------------------------------
        %RSA = cell(nsub,nXt);
        RSA = cell(nsub,1);
        for i_sub = 1:nsub
            % Convert the data to 3d:
            D = reshape(Y{i_sub}', [nmodes, size(Y{i_sub}, 1)/(nconditions), nconditions]);
            % Calculate RSA for current subject
            RSA(i_sub,:) = spm_eeg_rsa_specify(S,D);
            %RSA(i_sub,:) = spm_eeg_rsa_estimate(RSA(i_sub,:));
            RSA{i_sub} = spm_eeg_rsa_estimate(RSA{i_sub});
            
            % Reshape the ground truth beta to compare to the estimates
            % (the vRSA pipeline returns estimates in a different order)
            b = reshape(var(B_true{i_sub}, [], 2), [size(CV)])'; 
            % Compute the correlation between estimated parameters and
            % beta variance for each contrast
            Eps_beta_corr(pV_i, pE_i, i_sub) = corr(b(:), exp(RSA{1}.Ep.cond)');
        end

        if nsub > 1
            %------------------------------------------------------------------
            % Group level evidence with PEB
            %------------------------------------------------------------------
            try
                [~,F] = spm_eeg_rsa_peb(RSA, params=1:size(c, 2), FIR_bf=false, doplot=false);
            catch
                F_on(pV_i, pE_i) = nan;
                % Evidence against off effect (hence the -F)
                F_off(pV_i, pE_i) = nan;
                continue
            end

            %------------------------------------------------------------------
            % Gather evidence for "on" and "off" effects
            %------------------------------------------------------------------
            % Evidence in favor of on effect:
            F_on(pV_i, pE_i) = sum(F(find(CV' == 1, 1)));     % PZ: what does this do??
            % Evidence against off effect (hence the -F)
            F_off(pV_i, pE_i) = sum(F(find(~CV' == 1, 1)));
        else
            % Single subject analysis
            F_on(pV_i, pE_i)  = sum(RSA{1}.logBF.cond(spm_vec(CV')==1));
            F_off(pV_i, pE_i) = sum(RSA{1}.logBF.cond(spm_vec(CV')==0));
        end
    end
end

%% 3) Find the (pE, pV) combination that maximizes free energy for on effects while avoiding false positives
% Take the sum of the evidence for on and off effects:
sum_F = F_on - F_off;
[x,y]=find(sum_F==max(sum_F(F_off < 0), [], "all"));

% Optimized prior expectation & variance
pE = pEs(y);
pV = pVs(x);

if ~options.doplot
    return
end

%% Plot the results of the grid search:

% Convert free energies to probabilities
P_on = [];
P_off = [];
for pE_i = 1:length(pEs)    
    for pV_i = 1:length(pVs)
        F    = [F_on(pV_i,pE_i),0]';
        P    = spm_softmax(F);        
        P_on(pV_i,pE_i) = P(1);
        
        F    = [F_off(pV_i,pE_i),0]';
        P    = spm_softmax(F);        
        P_off(pV_i,pE_i) = P(1);        
    end
end

% Create x and y labels:
lbl_pE  = arrayfun(@(x) sprintf('%0.2f', x), pEs,...
    'UniformOutput', false); 
lbl_pV  = arrayfun(@(x) sprintf('%0.2f', x), pVs,...
    'UniformOutput', false); 

fig = figure('Name', 'Prior Selection', 'Position', [10, 10, 725, 725], 'Color',[1 1 1]);

% Layout for various subplots
tiledlayout(fig, 3, 2, 'TileSpacing','compact', 'Padding','compact');

% Plot the effects on (F):
nexttile()
imagesc(F_on); xticks(1:length(lbl_pE)); yticks(1:length(lbl_pV)); xticklabels(lbl_pE); yticklabels(lbl_pV)
xlabel("pE"); ylabel("pV"); title('Effect on (F)'); axis square;
colormap gray; colorbar; set(gca, 'FontSize',12);

% P:
nexttile()
imagesc(P_on); xticks(1:length(lbl_pE)); yticks(1:length(lbl_pV)); xticklabels(lbl_pE); yticklabels(lbl_pV)
xlabel("pE"); ylabel("pV"); title('Effect on (P)'); axis square;
colormap gray; colorbar; set(gca, 'FontSize',12);

% Plot the effects off (F):
nexttile()
imagesc(F_off); xticks(1:length(lbl_pE)); yticks(1:length(lbl_pV)); xticklabels(lbl_pE); yticklabels(lbl_pV)
xlabel("pE"); ylabel("pV"); title('Effect off (F)'); axis square;
colormap gray; colorbar; set(gca, 'FontSize',12);

% Plot the effects off (P):
nexttile()
imagesc(P_off); xticks(1:length(lbl_pE)); yticks(1:length(lbl_pV)); xticklabels(lbl_pE); yticklabels(lbl_pV)
xlabel("pE"); ylabel("pV"); title('Effect off (P)'); axis square;
colormap gray; colorbar; set(gca, 'FontSize',12);

% Correlation between Ep and beta variance as a function of pE and pE:
fig = figure('Name', 'Face validity', 'Position', [10, 10, 725, 725], 'Color',[1 1 1]);
imagesc(mean(Eps_beta_corr, 3)); xticks(1:length(lbl_pE)); yticks(1:length(lbl_pV)); xticklabels(lbl_pE); yticklabels(lbl_pV)
xlabel("pE"); ylabel("pV"); title('Correlation between Ep and Beta variance (r)'); axis square;
colormap gray; colorbar; set(gca, 'FontSize',12);

% Plot the distribution:
spm_plot_lognormal(pE, pV)
end