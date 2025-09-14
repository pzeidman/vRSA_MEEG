%% vRSA Example Script: vRSA with FIR basis set
%
% This script demonstrates how to perform 
% variational Representational Similarity Analysis (vRSA) using a fourier
% basis sete to characterize the contribution of experimental
% contrasts on the second order statistics of the data
% 
% *vRSA Overview*
% - Models the observed trial-by-trial covariance as a linear combination 
%   of model covariance components, each scaled by a log-space parameter 
%   (lambda).
% - Variational Laplace provides posterior estimates of these parameters 
%   and an approximation to the model evidence (free energy).
%
% *Fourier basis set*
% - For each trial/condition/channel, we first obtain beta parameters by 
%   projecting the data onto specified basis functions
% - In this example, the basis functions are fourier set
% - We then construct a trial-by-trial covariance matrix from those betas 
%   and decompose it into user-specified model covariance structures.
%
% *Single subject vRSA*
% - The trial by trial beta covariance matrices capture the trial by trial
%   covariance of the temporal dynamics encoded by each basis function in
%   the basis set. 
% - For each subject, each covariance matrix (i.e. for each basis function)
%   is decomposed as a linear mixture of the experimental contrasts using a
%   GLM. The lambda weights are estimated using restricted maximum
%   likelihood, yielding a posterior distribution for each 
% 
% *Group-Level Analysis*
% - Parametric Empirical Bayes (PEB) then combines these subject-level results 
%   into a group-level analysis.
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
subjects = {'S1','S2','S3','S4','S5','S6','S7', 'S8', 'S9','S10'};
subject = 1;
D = spm_eeg_load(sprintf('subjects/RmD_%s.mat',subjects{subject}));
[nmodes,ntimes,nstimuli] = size(D(:,:,:));  
rng(51); % Set seed

%% RSA settings
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

% Create fourier basis set with 16 basis functions
xBF = struct();
xBF.dt = 1 / D.fsample;  % Sampling rate
xBF.name = 'Fourier set';
xBF.length = ntimes * xBF.dt;
xBF.order  = 8;
xBF = spm_get_bf(xBF);    
Xt = xBF.bf(1:ntimes, 2:end);  % Remove intercept term

% RSA settings
S = struct();
S.Xt = Xt;
S.con_c = mat2cell(c, size(c, 1), ones(1, size(c, 2)));
S.con_c_names = cnames;
if exist('priors.mat', 'file')
    load('priors.mat')
    S.pE = priors.pE;
    S.pV = priors.pV;
else
    S.pE = -8;
    S.pV = 4;
end
%% Fit an examplar subject
RSA = spm_eeg_rsa_specify(S,D);
RSA = spm_eeg_rsa_estimate(RSA);
save(sprintf('subjects/RSA_s%d.mat',subject),'RSA','-v7.3');
% Computing the time axis of the FIR:
spm_eeg_rsa_review(RSA, FIR_bf=false, t=D.time', data=D(:, :, :));
%% Run first level analysis on group-level empirical data
RSAs = cell(length(subjects),1);
for s = 1:length(subjects)
    disp(subjects{s});
    % Load data    
    D = spm_eeg_load(sprintf('subjects/RmD_%s.mat',subjects{s}));
    % RSA settings
    settings = S;
    % Specify and estimate
    RSAs{s} = spm_eeg_rsa_specify(settings,D);    
    RSAs{s} = spm_eeg_rsa_estimate(RSAs{s});
end
save('subjects/RSAs.mat','RSAs','-v7.3');
%% Run second level analysis
load('subjects/RSAs.mat');

% Run PEB / BMC
[PEB,F] = spm_eeg_rsa_peb(RSAs, params=1:(length(RSAs{1}.qtypes)-1), FIR_bf=false, t=D.time');

save('subjects/PEB.mat','F','PEB');