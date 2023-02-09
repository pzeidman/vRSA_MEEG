% Runs an example vRSA analysis using the data from:
%
% 
% Miller, K.J., Schalk, G., Hermes, D., Ojemann, J.G. and Rao, R.P., 2016. 
% Spontaneous decoding of the timing and content of human object perception 
% from cortical surface recordings reveals complementary information in the 
% event-related potential and broadband spectral change. 
% PLoS computational biology, 12(1), p.e1004660.
%
% Peter Zeidman
%% Settings
subjects = {'ca','de','fp','ja','mv','wc','zt'};

% Load data from example subject
subject = 1;
D = spm_eeg_load(sprintf('subjects/CfdmD_%s.mat',subjects{subject}));

% Get sizes
[nmodes,ntimes,nstimuli] = size(D(:,:,:));  
%% RSA settings

% Define contrasts between conditions
is_face  = startsWith(D.conditions,'Face');
is_house = startsWith(D.conditions,'House');
c = zeros(nstimuli,1);
c(is_face) = 1;
c(is_house) = -1;

% Create and orthonormalise contrast vectors 
% (constant, faces vs houses, random)
rng(1);
C = [ones(nstimuli,1) c randn(nstimuli,1)];
C = spm_orth(C,'norm');

con_c  = {};
cnames = {};

% Genuine first condition
con_c{1}  = C(:,2);
cnames{1} = 'Face - House';    

% Random second condition
con_c{2}  = C(:,3);
cnames{2} = 'Emotional valence';

load('subjects/basis_set.mat');

% Within-trial basis set
%Xt = [ones(ntimes,1) bf];
Xt = bf;

% RSA settings
S = struct();
S.pE = -4;
S.pV = 8;
S.Xt = Xt;
S.con_c = con_c;
S.con_c_names = cnames;
%% Fit an examplar subject
RSA = spm_eeg_rsa_specify(S,D);
RSA = spm_eeg_rsa_estimate(RSA);
save(sprintf('subjects/RSA_s%d.mat',subject),'RSA','-v7.3');

spm_eeg_rsa_review(RSA);
%% Run first level analysis on group-level empirical data
nbases = size(Xt,2);
RSAs = cell(length(subjects),nbases);

for s = 1:length(subjects)
    disp(subjects{s});
    
    % Load data    
    D = spm_eeg_load(sprintf('subjects/CfdmD_%s.mat',subjects{s}));
    
    % RSA settings
    settings = S;
    
    % Specify and estimate
    RSAs(s,:) = spm_eeg_rsa_specify(settings,D);    
    RSAs(s,:) = spm_eeg_rsa_estimate(RSAs(s,:));
end
save('subjects/RSAs.mat','RSAs','-v7.3');
%% Run second level analysis
load('subjects/RSAs.mat');

% Run PEB / BMC
[PEB,F] = spm_eeg_rsa_peb(RSAs);

save('subjects/PEB.mat','F','PEB');