%% vRSA Example Script: Simulation
%
% This scripts demonstrates how to test the face validity of the vRSA
% pipeline with the selected priors. Data are simulated with known effects
% (and dimensions matching the data to analyze) and the vRSA pipeline is
% run. As the ground truth effects are known, we can test the capability of
% the pipeline to retrieve the parameters.
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

% Create FIR with bins of 50ms:
xBF = struct();  
xBF.dt = 1 / D.fsample;  % Sampling rate
xBF.name = 'Finite Impulse Response';
xBF.length = size(D, 2) * xBF.dt;
xBF.order  = 16;
xBF = spm_get_bf(xBF);
Xt = xBF.bf(1:size(D, 2), 2:end);

% RSA settings:
S = struct();
S.Xt = Xt;
S.con_c = mat2cell(c, size(c, 1), ones(1, size(c, 2)));
S.con_c_names = cnames;
S.pE = -4;
S.pV = 2;

%% Run vRSA on simulated data

% Simulate the data:
CV = zeros(size(Xt, 2), size(c, 2)); % Controls which effects are on and which are off: time bin x contrast
CV(4, 5) = 1; % Turn the first contrast on in 3r to 5th time window
CV(6, 1) = 1; % Turn the second contrast on in the 6th to 8th time window
s = 0.14;
Y = spm_eeg_simulate_covariance(Xt, c, s, nmodes, length(subjects), CV);

nbases = size(Xt,2);
RSAs = cell(length(subjects),nbases);
for s = 1:length(subjects)
    % Load data    
    y =  reshape(Y{s}', [nmodes, size(Y{s}, 1)/(nstimuli), nstimuli]);
    % RSA settings
    settings = S;
    % Specify and estimate
    RSAs(s,:) = spm_eeg_rsa_specify(settings,y);    
    RSAs(s,:) = spm_eeg_rsa_estimate(RSAs(s,:));
end
% Save the results:
save('subjects/RSAs-sim.mat','RSAs','-v7.3');
% Run PEB / BMC
[PEB,F] = spm_eeg_rsa_peb(RSAs, params=[1, 2, 3, 4, 5], FIR_bf=true, t=D.time');

% Save the results
save('subjects/PEB-sim.mat','F','PEB');