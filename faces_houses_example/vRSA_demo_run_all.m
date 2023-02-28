% Example vRSA analysis of EEG data from:
 
% Miller, K.J., Schalk, G., Hermes, D., Ojemann, J.G. and Rao, R.P., 2016. 
% Spontaneous decoding of the timing and content of human object perception 
% from cortical surface recordings reveals complementary information in the 
% event-related potential and broadband spectral change. 
% PLoS computational biology, 12(1), p.e1004660.
% 
% Data downloaded from:
% https://purl.stanford.edu/xd109qh3109
% 
% Peter Zeidman
%% Settings
subjects = {'ca','de','fp','ja','mv','wc','zt'};

% Initialize SPM in EEG mode
spm('defaults','eeg');
spm_jobman('initcfg');
%% Run analyses

% Pre-process ECoG data
vRSA_demo_preprocess;

% Generate an informed basis set derived from DCMs
vRSA_demo_create_basis_set;

% Run vRSA analyses with different priors to optimise the model (SLOW)
%vRSA_demo_select_priors;

% Run vRSA analysis
vRSA_demo_run_vRSA;
%% Load example subject's data. From here onwards it's all plotting.

% Load data from example subject
subject = 1;
D = spm_eeg_load(sprintf('subjects/CfdmD_%s.mat',subjects{subject}));

% Get sizes
[nmodes,ntimes,nstimuli] = size(D(:,:,:));  
%% Illustrate data
% Gray lines are individual channels, red lines are the first PC
is_house = startsWith(D.conditions,'House');

% Load reduced data (1st PC)
rD = spm_eeg_load(sprintf('subjects/RCfdmD_%s_2conditions_1channel.mat',subjects{subject}));

figure;
ax1=subplot(1,2,1);
y = D(:,:,:);
y = mean(y(:,:,is_house),3); % average over trials
ry = squeeze(rD(1,:,2));
plot(D.time,y','Color',[0.5 0.5 0.5]); hold on
plot(rD.time,ry,'r');
xline(0);
title('House stimuli');
xlabel('Time (secs)');
set(gca,'FontSize',12);
axis square

ax2=subplot(1,2,2);
y = D(:,:,:);
y = mean(y(:,:,~is_house),3); % average over trials
ry = squeeze(rD(1,:,1));
plot(D.time,y','Color',[0.5 0.5 0.5]); hold on
plot(rD.time,ry,'r');
xline(0);
title('Face stimuli');
xlabel('Time (secs)');
set(gca,'FontSize',12);
axis square

linkaxes([ax1 ax2]);
% figure;
% y = D(:,:,:);
% y = mean(y(:,:,~is_house),3); % average over trials
% for i = 1:D.nchannels
%     subplot(5,6,i);
%     plot(D.time,y(i,:));
% end

%% Illustrate the DCM fits
GCMs = {};

load('subjects/DCM/GCM_allsubjects_ons_150ms');
GCMs{1} = GCM;

load('subjects/DCM/GCM_allsubjects_ons_300ms');
GCMs{2} = GCM;

pst = GCMs{1}{1}.xY.pst;

conditions = {'face','house'};
onsets = {'150ms','300ms'};

figure;
plotcounter = 1;
for k = 1:2 % condition
    for j = 1:2 % onset time
        for i = 1:length(subjects)
            subplot(4,length(subjects),plotcounter);

            % DCM with onset time j subject i
            DCM = GCMs{j}{i};

            plot(pst,DCM.H{k},'LineWidth',2); hold on
            plot(pst,DCM.H{k}+DCM.R{k},'--','LineWidth',2,'Color',[0.5 0.5 0.5]);
            title(subjects{i});
            ylabel(sprintf('%s %s',conditions{k},onsets{j}),'FontSize',12);
            axis square
            xline(0);

            ylim([-4 4]);
            plotcounter = plotcounter + 1;
        end
    end
end
%% Ilustrate basis set
load('subjects/basis_set.mat');

%bf = [ones(size(bf,1),1) bf];

figure
plot(pst,bf);
xlabel('Time (ms)');
set(gca,'FontSize',12);
title('Basis functions','FontSize',16);

nbf = size(bf,2);
figure;
ax = [];
for i = 1:nbf
    ax(i)=subplot(nbf,1,i);
    plot(pst,bf(:,i));
    xline(0);
    set(gca,'FontSize',12);
    if i == 1
        title('Basis functions','FontSize',16);
    end
end
linkaxes(ax);
xlabel('Time (ms)');
%% Illustrate 2nd order parameters
load('subjects/basis_set.mat');

nbf = size(bf,2);

figure;
ax = [];
for i = 1:nbf
    BB = RSA{i}.BB;
    
    ax(i) = subplot(1,nbf,i);
    imagesc(BB);
    colormap gray;
    axis square;
end
linkaxes(ax);
%% Illustrate basis set fit to one subject's per-channel data
load('subjects/basis_set.mat');

is_house = startsWith(D.conditions,'House');

figure;
y = D(:,:,:);
y = mean(y(:,:,~is_house),3); % average over trials
y = y - mean(y,2);

h = [];
for i = 1:D.nchannels
    
    % Fit GLM
    X = full(bf);
    b = X \ y(i,:)';
    
    fitted = X * b;
    
    h(i) = subplot(4,7,i);
    plot(D.time*1000,y(i,:),'k--');
    hold on;
    plot(D.time*1000,fitted);
    xline(0);
    axis square;
    title(sprintf('Channel %d',i));
    set(gca,'YTickLabel',[]);
end
%linkaxes(h,'y');
%% Plot actual second order parameters and fitted, for a single subject
load('subjects/RSA_s1.mat');
nbases = length(RSA);
c = 1;

% Find min and max values of matrices excluding leading diagonal
minG = []; maxG = [];
for i = 1:nbases
    G = RSA{i}.G;
    G = G - diag(diag(G));
    maxG = [maxG max(G(:))];
    minG = [minG min(G(:))];
end
maxG = max(maxG);
minG = min(minG);
    
figure;
for i = 1:nbases
             
    % Prior and posterior variance
    pV = diag(RSA{i}.M.pC);
    Cp = diag(RSA{i}.Cp);
    
    % Priors and posteriors for grouped bar plot
    P = [spm_vec(RSA{i}.M.pE) spm_vec(RSA{i}.Ep)];
    V = [pV Cp];
     
    subplot(nbases,3,c);
    imagesc(RSA{i}.BB - diag(diag(RSA{i}.BB)));        
    caxis([minG maxG]);
    axis square; colormap gray; colorbar
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    c = c + 1;
    
    subplot(nbases,3,c);
    imagesc(RSA{i}.G - diag(diag(RSA{i}.G)));    
    caxis([minG maxG]);
    axis square; colormap gray; colorbar
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    c = c + 1;    
    
    subplot(nbases,3,c);
    spm_plot_ci(P,V);

    set(gca,'XTickLabel',{'Faces-Houses','Valence','Noise'});
    axis square    
    c = c + 1;        
end
%% Plot priors and posteriors from one subject
load('subjects/RSA_s1.mat');
figure;

subplot(1,2,1);
x = -3:0.1:30;
y1 = spm_Npdf(x, RSA{1}.M.pE.cond(1), RSA{1}.M.pC(1,1));
y2 = spm_Npdf(x, RSA{1}.Ep.cond(1), RSA{1}.Cp(1,1));
area(x,y1);
hold on
area(x,y2);
legend({'Prior','Posterior'});
xlabel('\lambda'); ylabel('P(\lambda)');
set(gca,'FontSize',12);
title('Faces-houses','FontSize',16);
axis square

subplot(1,2,2);
x = -3:0.1:30;
y1 = spm_Npdf(x, RSA{1}.M.pE.noise, RSA{1}.M.pC(2,2));
y2 = spm_Npdf(x, RSA{1}.Ep.noise, RSA{1}.Cp(2,2));
area(x,y1);
hold on
area(x,y2);
legend({'Prior','Posterior'});
xlabel('\lambda'); ylabel('P(\lambda)');
set(gca,'FontSize',12);
title('Noise','FontSize',16);
axis square
%% Plot within-subject model evidences
load('subjects/RSA_s1.mat');
spm_eeg_rsa_review(RSA);
%% Plot reconstructed ERP (single subject level)
load('subjects/basis_set.mat');
load('subjects/RSA_s1.mat');

% Load reduced data (1st PC)
rD = spm_eeg_load(sprintf('subjects/RCfdmD_%s_2conditions_1channel.mat',subjects{subject}));
y_faces  = squeeze(rD(:,:,1));
y_houses = squeeze(rD(:,:,2));

figure;
ax1=subplot(1,2,1);
h1=plot(pst,y_faces,'Color',[0.5 0.5 0.5]); hold on
h2=plot(pst,y_houses,'--','Color',[0.5 0.5 0.5]);
h3=plot(pst,y_faces-y_houses,'Color','r');
%xline(0); yline(0);
grid on
legend([h1,h2,h3],{'faces','houses','difference'});
xlabel('Time (ms)');
axis square;
set(gca,'FontSize',12);
title('ERPs (first PC)','FontSize',14);

ax2=subplot(1,2,2);
plot(pst,RSA{1}.y(:,1),'Color','b'); hold on
plot(pst,RSA{1}.y(:,2),'Color','k');
%xline(0); yline(0);
grid on
xlabel('Time (ms)');
axis square;
set(gca,'FontSize',12);
legend({'Faces-Houses','Random'});
title('Reconstructed','FontSize',14);

%% Plot results of search over different priors
load('subjects/priors_selection.mat');

figure;
subplot(2,1,1);
boxplot(Fs_on,'labels',labels);
yline(3);
xlabel('Prior variance');ylabel('dF');
set(gca,'FontSize',12);
title('Evidence for real effect','FontSize',16);

subplot(2,1,2);
boxplot(Fs_off,'labels',labels);
yline(3);
xlabel('Prior variance');
ylabel('dF');
set(gca,'FontSize',12);
title('Evidence for false effect','FontSize',16);