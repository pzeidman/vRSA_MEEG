% Generates an informed basis set for ERP analysis based on the data of:
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

addpath('../toolbox');

% Make output folder
out_dir = fullfile('subjects','DCM');
if ~exist(out_dir,'file')
    mkdir(out_dir);
end

%% Fit DCMs to each subject's data

% Estimate with 150ms onset
GCM = cell(length(subjects),1);
for s = 1:length(subjects)
    GCM{s} = specify_fit_dcm(subjects{s},150);
end
save('subjects/DCM/GCM_allsubjects_ons_150ms.mat','GCM');

% Average
BPA = spm_dcm_bpa(GCM,'true');
save('subjects/DCM/BPA_allsubjects_150ms.mat','BPA');

% Estimate with 300ms onset
GCM = cell(length(subjects),1);
for s = 1:length(subjects)
    GCM{s} = specify_fit_dcm(subjects{s},300);
end
save('subjects/DCM/GCM_allsubjects_ons_300ms.mat','GCM');

% Average
BPA = spm_dcm_bpa(GCM,'true');
save('subjects/DCM/BPA_allsubjects_300ms.mat','BPA');

%% Prepare within-trial design matrix (DCM-derived basis set)
BPAs = cell(2,1);

load('subjects/DCM/BPA_allsubjects_150ms');
BPAs{1} = BPA;

load('subjects/DCM/BPA_allsubjects_300ms');
BPAs{2} = BPA;

% We will just use the modelled voltage (not conductance)
% for each population. See spm_fx_cmc.m
% 1 = spiny stellate
% 3 = superficial pyramidal cells
% 5 = inhibitory interneurons
% 7 = deep pyramidal cells
states = [1 3 5 7];

[bf,pst] = spm_dcm_eeg2bf(BPAs,states);

save('subjects/basis_set.mat','bf','pst');