function DCM = specify_fit_dcm(subject, ons)
% Specifies and estimates a single-region DCM for the vRSA example
%
% subject - name of subject
% ons     - prior expectation for onset of stimulus (ms)
%
% Peter Zeidman

% Data and analysis directories
%--------------------------------------------------------------------------
Pbase     = fullfile(pwd,'subjects');

Pdata     = fullfile(Pbase, '.');   % data directory in Pbase
Panalysis = fullfile(Pbase, 'DCM'); % analysis directory in Pbase

% Data filename
%--------------------------------------------------------------------------
fn = sprintf('RCfdmD_%s_2conditions_1channel.mat',subject);
DCM.xY.Dfile = fullfile(Pdata,fn);

% Parameters and options used for setting up model
%--------------------------------------------------------------------------
DCM.options.analysis = 'ERP'; % analyze evoked responses
DCM.options.model    = 'CMC'; % CMC model
DCM.options.spatial  = 'LFP'; % spatial model
DCM.options.trials   = [1 2]; % index of ERPs within ERP/ERF file
DCM.options.Tdcm(1)  = -50;   % start of peri-stimulus time to be modelled
DCM.options.Tdcm(2)  = 600;   % end of peri-stimulus time to be modelled
DCM.options.Nmodes   = 1;     % nr of modes for data selection
DCM.options.h        = 0;     % nr of DCT components
DCM.options.onset    = ons;   % selection of onset (prior mean)
DCM.options.dur      = 8;     % prec
DCM.options.D        = 1;     % downsampling

%--------------------------------------------------------------------------
% Data and spatial model
%--------------------------------------------------------------------------
DCM  = spm_dcm_erp_data(DCM);

%--------------------------------------------------------------------------
% Location priors for dipoles
%--------------------------------------------------------------------------
DCM.Lpos  = [0; 0; 0];
DCM.Sname = {'ROI'};

%--------------------------------------------------------------------------
% Spatial model
%--------------------------------------------------------------------------
DCM = spm_dcm_erp_dipfit(DCM);

%--------------------------------------------------------------------------
% Specify connectivity model
%--------------------------------------------------------------------------
start_dir = pwd;
cd(Panalysis)

% Fwd (sp -> ss)
DCM.A{1} = 0;

% Bwd (sp -> dp)
DCM.A{2} = 0;

% Modulatory (dp -> sp)
DCM.A{3} = 1;

% B-effect
DCM.B{1} = 1;

% Input
DCM.C = 1;

%--------------------------------------------------------------------------
% Between trial effects
%--------------------------------------------------------------------------
DCM.xU.X = [1; -1];
DCM.xU.name = {'faces vs houses'};

%--------------------------------------------------------------------------
% Invert
%--------------------------------------------------------------------------
DCM.name = ['DCM_' subject];

DCM = spm_dcm_erp(DCM);

cd(start_dir);