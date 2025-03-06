% Example vRSA analysis of EEG data
% 
% In this script, we perform the entire pipeline of the vRSA toolbox. To
% investigate which experimental contrasts are experessed in the second
% order statistics of the data (i.e. which model matrices are are
% represented in the multivariate patterns of the data in representational
% analysis (RSA) linguo), we perform the following steps:
%
% 1. Preprocessing of the epoched data
%   - Average across trials in each condition (i.e. stimulus examplar)
%   - Dimensionality reduction with PCA across channels, 7 components
% 2. Select Priors
%   - Simulate data with known effects
%   - Run the vRSA pipeline to identify priors maximizing sensitivity and
%   specificity
% 3. Test vRSA on simulated data (face validity)
%   - Simulate data with effects of 2 experimental contrasts at 2 time
%   - Test that the expected effects are retrieved
%   points
% 4. Run VRSA with selected priors
%   - Compute vRSA on each subject
%   - Group level analysis using spm_peb
% 
% In addition, the script plots all the figures presented in the paper
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

% Initialize SPM in EEG mode
spm('defaults','eeg');
spm_jobman('initcfg');
%% Run analyses
% Pre-process EEG data
vRSA_demo_preprocess;

% Select the priors (Takes a while)
if ~exist('priors.mat', 'file')
    vRSA_demo_select_priors;
end
% Run simulation (to show face validity of the pipeline with selected
% priors)
vRSA_demo_run_vRSA_simulation;

% Run vRSA analysis with FIR
vRSA_demo_run_vRSA_FIR;

close all
%% Figures for the paper:
% From here on, it is only plotting of the figuresd of the paper:
if ~exist('./figures', 'dir')
    mkdir('./figures')
end

% Load one example data set:
D = spm_eeg_load(sprintf('subjects/RmD_%s.mat',subjects{subject}));
[nmodes,ntimes,nstimuli] = size(D(:,:,:));

%% Figure 1:
% Plot of the normal and log normal distribution with specified parame=ters
spm_plot_lognormal(-8, 16)
% Fetch the figure handle and save it:
fh = findobj( 'Type', 'Figure', 'Name', 'Normal Log normal' );
saveas(fh, './figures/Figure1.svg')

%% Figure 3:

% Figure 3A: plot the first components across trials
fig = figure;
h = plot(D.time, squeeze(D(1, :, :)), '-o', 'MarkerSize', 1.5);
set(h, {'MarkerFaceColor'}, get(h,'Color'));
xline(0:D.time(3) - D.time(1):0.5, 'LineWidth', 0.5)
lbl = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"};
t = spm_pinv(Xt) * D.time';
for i = 1:length(t)
    text(t(i), 0.9, lbl{i}, 'HorizontalAlignment', ...
        'center', 'FontSize', 12)

end
xlabel('Time (secs)')
ylabel('Amplitude')
set(gca, 'FontSize',12);
saveas(fig, './figures/Figure3A.svg')

% Figure 3B: plot the basis functions
fig = figure;
Xt_plot = Xt + repmat(0:size(Xt,2)-1, size(Xt,1), 1)*1.3;
plot(Xt_plot', D.time, 'k')
set(gca, 'YDir','reverse')
xlim([-0.5, 19.5])
ax1 = gca;
ax1.XAxis.Visible = 'off'; % remove x-axis
set(gca,'YAxisLocation','left', 'box','off')
ylabel("Time (secs)")
set(gca, 'FontSize',12);
saveas(fig, './figures/Figure3B.svg')

%% Figure 4:

% Figure 4A: plot the Within trial contrast vectors
fig = figure;
nXt = size(Xt, 2);
fig = figure;
t = tiledlayout(fig, 1, nXt,'TileSpacing', 'compact');
for i = 1:nXt
    nexttile(t)
    imagesc(Xt(:, i))
    colormap gray;
    if i == 1
        ylabel('Time')
    else
        yticks([])
        yticklabels([])
    end
    xticks([])
    xticklabels([])
    set(gca, 'FontSize',12);
end
saveas(fig, './figures/Figure4A.svg')

% Figure 4B: plot the between trials contrast vectors:
c = cell2mat(RSA{1}.con);
nC = size(c, 2);
fig = figure;
t = tiledlayout(fig, 1, nC,'TileSpacing', 'compact');
for i = 1:nC
    nexttile(t)
    imagesc(c(:, i))
    colormap gray;
    if i == 1
        ylabel('Stimulus image')
    else
        yticks([])
        yticklabels([])
    end
    xticks([])
    xticklabels([])
    set(gca, 'FontSize',12);
end
saveas(fig, './figures/Figure4B.svg')

%% Figure 6:

% Plot the results of the simulation to showcase face validity:
load('./subjects/RSAs-sim.mat')
spm_eeg_rsa_peb(RSAs, params=[1, 2, 3, 4, 5], FIR_bf=true, t=D.time');

% Save the time resolved posterior and evidence:
fh = findobj( 'Type', 'Figure', 'Name', 'Time resolved Posterior and evidence');
saveas(fh, './figures/Figure6ABC.svg')
% Save the pooled evidence:
fh = findobj( 'Type', 'Figure', 'Name', 'Group vRSA analysis');
saveas(fh, './figures/Figure6DEF.svg')
close all
%% Figure 6:
% Same but for the real data:
% Plot the results of the simulation to showcase face validity:
load('./subjects/RSAs.mat')
spm_eeg_rsa_peb(RSAs, params=[1, 2, 3, 4, 5], FIR_bf=true, t=D.time');

% Save the time resolved posterior and evidence:
fh = findobj( 'Type', 'Figure', 'Name', 'Time resolved Posterior and evidence');
saveas(fh, './figures/Figure7ABC.svg')
% Save the pooled evidence:
fh = findobj( 'Type', 'Figure', 'Name', 'Group vRSA analysis');
saveas(fh, './figures/Figure7DEF.svg')
