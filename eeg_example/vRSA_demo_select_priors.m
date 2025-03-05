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
%    Miller, K.J., Schalk, G., Hermes, D., Ojemann, J.G. and Rao, R.P., 2016.
%    Spontaneous decoding of the timing and content of human object perception
%    from cortical surface recordings reveals complementary information in the
%    event-related potential and broadband spectral change.
%    PLoS Computational Biology, 12(1), p.e1004660.
% 
% Author:      Peter Zeidman, Alex Lepauvre
% Date:        2025-02-14
% Version:     1.0

%% Settings
subjects = {'S1','S2','S3','S4','S5','S6','S7', 'S8', 'S9','S10'};
subject = 1;
D = spm_eeg_load(sprintf('subjects/RmD_%s.mat',subjects{subject}));
[nmodes,ntimes,nstimuli] = size(D(:,:,:));

%% Create contrasts and basis set:
% Get binary vectors of each condition:
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

% Create FIR with bins of 50ms:
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
    % Compute the residuals variance relative to the fitted variance:
    s(sub_i) = var(res, [], 'all') / var(fitted, [], 'all');
end
% Average variance across subjects:
s = mean(s);

%% Search for optimal priors
% Grid search on simulated data to select priors maximizing sensitivity and
% specificity:
[pE, pV] = spm_eeg_rsa_select_priors(c, Xt, s, nmodes, nstimuli, ...
    length(subjects));

fprintf('\n Selected Priors: \n for v=%0.2f \n Prior Expectation=%0.2f, \n Prior variance=%0.2f\n', s, pE, pV)

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
print(fh, './figures/Figure5ABC.svg', '-dsvg')
% Plot of the selected distribution (export to pdf as svg messes font up)
fh = findobj( 'Type', 'Figure', 'Name', 'Normal Log normal' );
print(fh, './figures/Figure5D.svg', '-dsvg')