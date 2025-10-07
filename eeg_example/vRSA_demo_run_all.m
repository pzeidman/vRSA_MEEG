% Example vRSA analysis of EEG data
% 
% In this script, we perform the entire pipeline of the vRSA toolbox. To
% investigate which experimental contrasts are expressed in the second
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
% Data URL:
% https://purl.stanford.edu/bq914sc3730
%
% Author:      Peter Zeidman, Alex Lepauvre
% Date:        2025-02-14
% Version:     1.0

%% Settings
subjects = {'S1','S2','S3','S4','S5','S6','S7', 'S8', 'S9','S10'};
% subjects = {'S1'};

% Initialize SPM in EEG mode
spm('defaults','eeg');
spm_jobman('initcfg');
addpath('C:\Users\alexander.lepauvre\Documents\GitHub\vRSA_MEEG\toolbox')
%% Download the data
% NB if this download hangs, try using your browser to downlaod the data
% manually and save it in a folder called 'raw_data'

% Make output folder
out_dir = 'raw_data';
if ~exist(out_dir,'file')
    mkdir(out_dir);
end

disp('Downloading data...');
for i = 1:length(subjects)    
    fname = fullfile(out_dir,sprintf('%s.mat',subjects{i}));
    url   = sprintf('https://stacks.stanford.edu/file/bq914sc3730/%s.mat',subjects{i});
    websave(fname, url);    
end
disp('Done');
%% Run analyses
% Pre-process EEG data
vRSA_demo_preprocess(subjects);

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
% From here on, it is only plotting of the figures of the paper:
if ~exist('./figures', 'dir')
    mkdir('./figures')
end

% Load one example data set:
subject = 1;
D = spm_eeg_load(sprintf('subjects/RmD_%s.mat',subjects{subject}));
[nmodes,ntimes,nstimuli] = size(D(:,:,:));

% Load one example RSA model
load('subjects/RSA_s1.mat');
Xt = RSA.M.Xt;

%% Figures:
%% Figure 1:
% Plot of the normal and log normal distribution with specified parame=ters
mu = -8;
v  = 4;
x  = linspace(-30,30,1000);
L  = spm_Npdf(x,mu,v);
V  = spm_Npdf(exp(x),mu,v);

fh = figure;
tiledlayout(1,2);
nexttile
area(x,L,'FaceColor',[0.7 0.7 0.7]);
xlabel('lambda');ylabel('Probability density');axis square;
xline(0,':');
set(gca,'FontSize',12);
nexttile()
area(exp(x),V,'FaceColor',[0.7 0.7 0.7]);
xlabel('v');ylabel('Probability density');axis square;
xlim([-0.5 10]);
xline(0,':');
set(gca,'FontSize',12);
saveas(fh, './figures/Figure1.svg')

%% Figure 3:
% Figure 3A: plot the first components across trials
fig = figure('Position', [10, 10, 725, 500]);
tiledlayout(1,2);
nexttile
h = plot(D.time, squeeze(D(1, :, :)), '-o', 'MarkerSize', 1.5);
set(h, {'MarkerFaceColor'}, get(h,'Color'));
xline(0:D.time(3) - D.time(1):0.5, 'LineWidth', 0.5)
lbl = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"};
t = spm_pinv(Xt) * D.time';
for i = 1:length(t)
    text(t(i), 0.9, lbl{i}, 'HorizontalAlignment', ...
        'center', 'FontSize', 10)

end
xlabel('Time (secs)')
ylabel('Amplitude')
axis square
% Figure 3B: plot the basis functions
nexttile
imagesc(Xt);
ticks = get(gca,'YTick');
tick_labels = D.time(ticks);
set(gca,'YTick',ticks,'YTickLabel',compose('%.2f',tick_labels));
ylabel('Time (secs)')
xlabel('Basis function (covariate)');
colormap gray
set(gca, 'FontSize',12);
axis square
saveas(fig, './figures/Figure3.svg')

%% Figure 4
% Plot overall design matrix
fig=figure;imagesc(RSA.M.X)
colormap gray
ylim([0 120]);
xlim([0 50]);
xlabel('Covariate');
ylabel('Measurement');
set(gca, 'FontSize',12);
saveas(fig, './figures/Figure4.svg')

%% Figure 5: plot the between trials contrast vectors:
c = cell2mat(RSA.con);
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
saveas(fig, './figures/Figure5.svg')

%% Figure 6:
% Plot the results of the simulation to showcase face validity:
load('./subjects/RSAs-sim.mat')
[PEB,F,BB,G,F_bfs,Pp_bfs,F_contrasts,Pp_contrasts] = spm_eeg_rsa_peb(RSAs, FIR_bf=true, t=D.time');
close all
% Extract info:
ncon = length(RSAs{1}.con);
cnames = RSAs{1}.cnames;
nXt  = RSAs{1}.M.nXt;
qcon = RSAs{1}.qcon;
t = spm_pinv(Xt) * D.time';
bin_dur = D.time(sum(Xt(:, 2))) - D.time(1);
Ep = reshape(full(PEB.Ep), ncon, nXt);
Cp = PEB.Cp;

% Figure 6
fig = figure('Position', [10, 10, 1000, 725]);
tiledlayout(3, 4, 'TileSpacing','compact', 'Padding','compact');

% Figure 6A
% Plot the estimated parameters:
nexttile([1, 2])
ctr = 1;
for c = 1:ncon
    hold on
    spm_plot_ci(Ep(c,:), diag(Cp(c == qcon, c == qcon)), t);
    lbl{ctr}   = '';
    lbl{ctr+1} = cnames{c};
    ctr        = ctr + 2;
end
title('Estimated hyperparameters','FontSize',12);
ylabel('v');
xlabel('Time (sec)');
xlim([t(1), t(end)]);
set(gca, 'FontSize',12);

% Figure 6B
% Plot the exponentiated estimated parameters:
nexttile([1, 2])
for c = 1:ncon
    hold on
    spm_plot_ci(Ep(c,:), diag(Cp(c == qcon,c == qcon)), t, [], 'exp');
end
title('Estimated hyperparameters','FontSize',12);
ylabel('Weights (exp)');
xlabel('Time (sec)');
xlim([t(1), t(end)]);
set(gca, 'FontSize',12);

% Figure 6C
% Plot the log evidence per basis function and contrast
nexttile([1, 2])
h = plot(t, F, '-o');
set(h, {'MarkerFaceColor'}, get(h,'Color'));
title('Log evidence per basis function','FontSize',12);
ylabel('log BF');
xlabel('Time (sec)');
xlim([t(1), t(end)]);
legend(cnames, 'Location','best');
set(gca, 'FontSize',12);

% Figure 6D
% Plot log evidence for each basis function
nexttile([1, 2])
yyaxis right
hold on
b = bar(t, Pp_bfs, 'FaceColor', [0.3, 0.3, 0.3]);
b.FaceAlpha = 0.2;
b.EdgeColor = 'none';
ylabel('Probability');

yyaxis left
plot(t, F_bfs, 'LineWidth', 2);
ylabel('log BF');
xlabel('Time (sec)');
set(gca,'FontSize',12);
title('Log evidence over time','FontSize',12)
ax = gca;
ax.YAxis(2).Color = [0.3, 0.3, 0.3];

% Figure 6E
% Plot log evidence for each contrast
nexttile
bar(F_contrasts, 'FaceColor', [0.3, 0.3, 0.3]);
set(gca,'XTick',1:ncon,'XTickLabel',cnames);
ylabel('Log evidence (nats)');
set(gca,'FontSize',12);
title('Pooled evidence per component','FontSize',12);

nexttile
bar(Pp_contrasts, 'FaceColor', [0.3, 0.3, 0.3]);
set(gca,'XTick',1:ncon,'XTickLabel',cnames);
ylabel('Probability');
set(gca,'FontSize',12);
title('Probability per component','FontSize',12);

% Figure 6E
% Plot fitted second order model
nexttile
imagesc(BB);
colormap gray;
xlabel('Condition'); ylabel('Condition');
set(gca,'FontSize',12);
title('Betas (second order)','FontSize',12);
axis square;
nexttile
G(eye(size(G))==1) = nan;
imagesc(G);
xlabel('Condition'); ylabel('Condition');
set(gca,'FontSize',12);
title('vRSA','FontSize',12);
axis square;
saveas(fig, './figures/Figure6.svg')

%% Figure 7:
% Same but for the real data:
% Plot the results of the simulation to showcase face validity:
load('./subjects/RSAs.mat')
[PEB,F,BB,G,F_bfs,Pp_bfs,F_contrasts,Pp_contrasts] = spm_eeg_rsa_peb(RSAs, FIR_bf=true, t=D.time');
close all
% Extract info:
ncon = length(RSAs{1}.con);
cnames = RSAs{1}.cnames;
nXt  = RSAs{1}.M.nXt;
qcon = RSAs{1}.qcon;
t = spm_pinv(Xt) * D.time';
bin_dur = D.time(sum(Xt(:, 2))) - D.time(1);
Ep = reshape(full(PEB.Ep), ncon, nXt);
Cp = PEB.Cp;

% Figure 7
fig = figure('Position', [10, 10, 1000, 725]);
tiledlayout(3, 4, 'TileSpacing','compact', 'Padding','compact');

% Figure 7A
% Plot the estimated parameters:
nexttile([1, 2])
ctr = 1;
for c = 1:ncon
    hold on
    spm_plot_ci(Ep(c,:), diag(Cp(c == qcon, c == qcon)), t);
    lbl{ctr}   = '';
    lbl{ctr+1} = cnames{c};
    ctr        = ctr + 2;
end
title('Estimated hyperparameters','FontSize',12);
ylabel('v');
xlabel('Time (sec)');
xlim([t(1), t(end)]);
set(gca, 'FontSize',12);

% Figure 7B
% Plot the exponentiated estimated parameters:
nexttile([1, 2])
for c = 1:ncon
    hold on
    spm_plot_ci(Ep(c,:), diag(Cp(c == qcon,c == qcon)), t, [], 'exp');
end
title('Estimated hyperparameters','FontSize',12);
ylabel('Weights (exp)');
xlabel('Time (sec)');
xlim([t(1), t(end)]);
set(gca, 'FontSize',12);

% Figure 7C
% Plot the log evidence per basis function and contrast
nexttile([1, 2])
h = plot(t, F, '-o');
set(h, {'MarkerFaceColor'}, get(h,'Color'));
title('Log evidence per basis function','FontSize',12);
ylabel('log BF');
xlabel('Time (sec)');
xlim([t(1), t(end)]);
legend(cnames, 'Location','best');
set(gca, 'FontSize',12);

% Figure 7D
% Plot log evidence for each basis function
nexttile([1, 2])
yyaxis right
hold on
b = bar(t, Pp_bfs, 'FaceColor', [0.3, 0.3, 0.3]);
b.FaceAlpha = 0.2;
b.EdgeColor = 'none';
ylabel('Probability');

yyaxis left
plot(t, F_bfs, 'LineWidth', 2);
ylabel('log BF');
xlabel('Time (sec)');
set(gca,'FontSize',12);
title('Log evidence over time','FontSize',12)
ax = gca;
ax.YAxis(2).Color = [0.3, 0.3, 0.3];

% Figure 7E
% Plot log evidence for each contrast
nexttile
bar(F_contrasts, 'FaceColor', [0.3, 0.3, 0.3]);
set(gca,'XTick',1:ncon,'XTickLabel',cnames);
ylabel('Log evidence (nats)');
set(gca,'FontSize',12);
title('Pooled evidence per component','FontSize',12);

nexttile
bar(Pp_contrasts, 'FaceColor', [0.3, 0.3, 0.3]);
set(gca,'XTick',1:ncon,'XTickLabel',cnames);
ylabel('Probability');
set(gca,'FontSize',12);
title('Probability per component','FontSize',12);

% Figure 7E
% Plot fitted second order model
nexttile
imagesc(BB);
colormap gray;
xlabel('Condition'); ylabel('Condition');
set(gca,'FontSize',12);
title('Betas (second order)','FontSize',12);
axis square;
nexttile
G(eye(size(G))==1) = nan;
imagesc(G);
xlabel('Condition'); ylabel('Condition');
set(gca,'FontSize',12);
title('vRSA','FontSize',12);
axis square;
set(gca, 'FontSize',12);
saveas(fig, './figures/Figure7.svg')
