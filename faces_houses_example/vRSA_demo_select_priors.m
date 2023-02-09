% Performs vRSA analyses, varying the prior variance on the parameters, to
% select the prior variance that maximize the evidence for real effects and
% minimizes the evidence for fictious effects. This is done using
% simulations based on the paradigm of:
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

% Initialize SPM in EEG mode
spm('defaults','eeg');
spm_jobman('initcfg');

% Load data from example subject
subject = 1;
D = spm_eeg_load(sprintf('subjects/CfdmD_%s.mat',subjects{subject}));

% Get sizes
[nmodes,ntimes,nstimuli] = size(D(:,:,:));  
%% Sensitivity / specificity analysis

% Define contrasts between conditions
is_face  = startsWith(D.conditions,'Face');
is_house = startsWith(D.conditions,'House');
c = zeros(nstimuli,1);
c(is_face) = 1;
c(is_house) = -1;

load('subjects/basis_set.mat');

% Within-trial basis set
Xt = bf;

pVs = [1/64 1/32 1/16 1/8 1/4 1/2 1 2 4 8 16 32];
labels = {'1/64','1/32','1/16','1/8','1/4','1/2','1','2','4','8','16','32'};
n_it = 12;

rng(1);
js = 1:length(pVs);

Fs_on  = zeros(n_it,length(pVs));
Fs_off = zeros(n_it,length(pVs));

nbases = size(Xt,2);

for i = 1:n_it

    con_c     = {};
    cnames    = {};
    
    % Create and orthonormalise contrast vectors 
    % (constant, faces vs houses, random)
    C = [ones(nstimuli,1) c randn(nstimuli,1)];
    C = spm_orth(C,'norm');
    
    % Genuine first condition
    con_c{1}  = C(:,2);
    cnames{1} = 'Face - House';    
    
    % Random second condition
    con_c{2}  = C(:,3);
    cnames{2} = 'Emotional valence';
    
    % RSA settings
    S = struct();
    S.pE = -4;
    S.Xt = Xt;    
    S.con_c = con_c;
    S.con_c_names = cnames;
                    
    % Prior variance
    for j = js
        S.pV = pVs(j);

        RSA = cell(length(subjects),nbases);
        parfor s = 1:length(subjects)
            disp(subjects{s});

            % Load data    
            D = spm_eeg_load(sprintf('subjects/CfdmD_%s.mat',subjects{s}));

            % Specify and estimate
            RSA(s,:) = spm_eeg_rsa_specify(S,D);    
            RSA(s,:) = spm_eeg_rsa_estimate(RSA(s,:));
        end

        % Run PEB / BMC
        [~,F] = spm_eeg_rsa_peb(RSA, [1 2]);    

        % Evidence for face-house regressor
        Fs_on(i,j) = sum(F(1,:));
        
        % Evidence for random regressor
        Fs_off(i,j) = sum(F(2,:));        
    end
end

save('subjects/priors_selection.mat','Fs_on','Fs_off','labels');
