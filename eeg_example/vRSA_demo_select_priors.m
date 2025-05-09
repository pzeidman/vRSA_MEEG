%% vRSA Example Script: Select priors
%
% This scripts demonstrates how to identify empirical priors that maximize
% sensitivity (i.e. detecting present effects) and specificity (rejecting
% absent effects). To do so, we simulate data with the same experimental
% design and parameters as the data we aim to analyze, but with known
% covariance structures. The vRSA pipeline is then applied while varying
% prior Expectation and Variance. The optimal priors are those that yield
% maximal evidence for the effects present in the data, and maximum
% (negative) evidence for the effects known to be absent in the data.
% 
% Step by Step:
% 1. Compute observation noise from the data to analyze
%    - For each subject, fit the within trial design matrix to each trial
%    - Compute the residuals variance (relative to the variance of the
%    fitted data)
% 2. Simulate data with known effects:
%    - Use a multivariate GLM with beta parameters selected from a normal
%    distribution across channels for some effects (contrast at given time
%    point) while setting others to 0
%    - The Simulated data have the same dimensions as the data to be
%    analyzed: same number of subjects, channels, time points, trials
% 3. Fit the vRSA while varying prior Expectations and Variance
%    - Sample pEs from -16:2:-2; and pVs from [2 4 8 16 32 64 128];
%    - Fit the vRSA with each combination
%    - Take the sum of evidence for on effects, and the sum of - evidence
%    for off effects. 
% 4. Identify optimal priors
%    - Normalize the on and off priors between 0 and 1, to account for the
%    difference in range between the two
%    - Select priors maximizing the of the normalized evidence
% 
% Step 2, 3 and 4 have been combined in a standalone function, to
% facilitate reuse of the procedure by other users
% 
% Data for this example are from:
%   Kaneshiro, Blair, et al. "A representational similarity analysis of the 
%   dynamics of object processing using single-trial EEG classification." 
%   Plos one 10.8 (2015): e0135697.
%
% Data downloaded from:
% https://purl.stanford.edu/bq914sc3730
% 
% Author:      Peter Zeidman, Alex Lepauvre
% Date:        2025-02-14
% Version:     1.0

%% Settings
% subjects = {'S1','S2','S3','S4','S5','S6','S7', 'S8', 'S9','S10'};
subjects = {'S1'};
subject = 1;
D = spm_eeg_load(sprintf('subjects/RmD_%s.mat',subjects{subject}));
[nmodes,ntimes,nstimuli] = size(D(:,:,:));
rng(51); % Set seed

%% Create contrasts and basis set:

% Get binary vectors for which stimuli belong to each condition
is_human = double(contains(D.conditions,'Human'));
is_animal = double(contains(D.conditions,'Animal'));
is_face = double(contains(D.conditions,'Face'));
is_body = double(contains(D.conditions,'Body'));
is_natural = double(contains(D.conditions,'Natural'));
is_manmade = double(contains(D.conditions,'ManMade'));

% Define contrasts:
c = zeros(nstimuli, 5);
c(:, 1) = is_human + is_animal - (is_natural + is_manmade); % Animate vs. inanimate
c(:, 2) = is_human - is_animal; % Human vs. animal
c(:, 3) = is_face - is_body; % Body vs. face
c(:, 4) = is_natural - is_manmade; % Natural vs. manmade
c(:, 5) = c(:, 2) .* c(:, 3); % Interaction between body part and species

% Add conditions names:
cnames = {};
cnames{1} = 'Animate - Inanimate';
cnames{2} = 'Human - Animal';
cnames{3} = 'Face - Body';
cnames{4} = 'Natural - Manmade';
cnames{5} = 'Species * Body';

% Create FIR with-trial design matrix with bins of 50ms:
xBF = struct();
xBF.dt = 1 / D.fsample;  % Sampling rate
xBF.name = 'Finite Impulse Response';
xBF.length = size(D, 2) * xBF.dt;
xBF.order  = 16;
xBF = spm_get_bf(xBF);
Xt = xBF.bf(1:size(D, 2), 2:end);

%% Estimate variance from the data:
s = zeros(length(subjects), 1);
% Create the design matrix:
X = kron(eye(nstimuli),Xt);
% Loop through each subject:
for sub_i = 1:length(subjects)
    % Load the data:
    D = spm_eeg_load(sprintf('subjects/RmD_%s.mat',subjects{sub_i}));
    % Stack trials vertically:
    splitY = num2cell(permute(D(:, :, :),[2 1 3]), [1 2]);
    Y = vertcat(splitY{:});
    % Estimate beta with OLS:
    B  = spm_pinv(X) * Y;
    % Compute fitted response and the residuals:
    fitted = X * B; res = Y - fitted;
    % Compute the residual standard deviation (averaged over channels)
    s(sub_i) = mean(std(res));
end

% Average residual s.d. across subjects:
s = mean(s);

% Plot sanity check (average betas over stimuli)
Bhat = B;
Yhat = zeros(ntimes,1);
mode = 1;
ncov = size(Xt,2);
for i = 1:nstimuli
    Yhat = Yhat + Xt * Bhat(1:ncov,mode);
    Bhat = Bhat((ncov+1):end,:);
end
Yhat = Yhat ./ nstimuli;
figure;plot(Yhat);
xlabel('Time (measurements)');
title({'Average modelled ERP';'(when estimating residual variance)'});
%% Search for optimal priors
% Grid search on simulated data to select priors maximizing sensitivity and
% specificity:
tic
[pE, pV] = spm_eeg_rsa_select_priors(c, Xt, s, nmodes, nstimuli, ...
    length(subjects));

fprintf('\n Selected Priors: \n for v=%0.2f \n Prior Expectation=%0.2f, \n Prior variance=%0.2f\n', s, pE, pV)
toc

% Save the priors for later use:
priors = struct();
priors.pE = pE;
priors.pV = pV;
save("priors.mat", 'priors')

%% Saving the figures:
if ~exist('./figures', 'dir')
    mkdir('./figures')
end

fh = findobj( 'Type', 'Figure', 'Name', 'Prior Selection' );
print(fh, './figures/Figure5ABCD.svg', '-dsvg')
fh = findobj( 'Type', 'Figure', 'Name', 'Face validity' );
print(fh, './figures/Figure5E.svg', '-dsvg')
% Plot of the selected distribution (export to pdf as svg messes font up)
fh = findobj( 'Type', 'Figure', 'Name', 'Normal Log normal' );
print(fh, './figures/Figure5F.svg', '-dsvg');