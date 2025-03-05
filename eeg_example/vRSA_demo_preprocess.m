function vRSA_demo_preprocess
% Imports EEG data from:
%
% Kaneshiro B, Perreau Guimaraes M, Kim HS, Norcia AM, and Suppes P (2015). A 
% Representational Similarity Analysis of the Dynamics of Object Processing 
% Using Single-Trial EEG Classification. PLoS ONE 10(8): e0135697. 
% doi:10.1371/journal.pone.0135697 
% 
% Data downloaded from:
% https://purl.stanford.edu/bq914sc3730
%
% Outputs:
%
% subjects/RmD_**.mat - channel-wise, stimulus-wise ERP data for vRSA
% 
% by Peter Zeidman & Alex Lepauvre
% _________________________________________________________________________

% Settings

% Path to folder with raw data
raw_dir = fullfile(pwd,'raw_data');

% Make output folder
out_dir = 'subjects';
if ~exist(out_dir,'file')
    mkdir(out_dir);
end

% Subject names
subjects = {'S1','S2','S3','S4','S5','S6','S7','S8', 'S9','S10'};
nsubjects = length(subjects);

% Name of each stimulus category in the right order:
categories = {'HumanBody', 'HumanFace', 'AnimalBody', 'AnimalFace', ...
    'Natural', 'ManMade'};

% Convert data to Fieldtrip and then to SPM format
% -------------------------------------------------------------------------
start_dir = pwd;

f = figure;
for s = 1:nsubjects 
    disp(subjects{s});
    
    % Prepare subject's data in Fieldtrip format
    % ---------------------------------------------------------------------
    
    % Load data (times x channels)
    load(fullfile(raw_dir, [subjects{s} '.mat']));    
    
    % Get the dimension:
    [nchannels, ntimes,ntrials] = size(X_3D);

    % The data are already epoched, we only need to generate a time vector:
    t = 0:1/Fs:0.5;
    
    % Set channel labels
    labels = cell(nchannels,1);
    for i = 1:nchannels
        labels{i} = num2str(i);
    end    
    ft = struct();
    ft.label = labels;

    % Create condition strings:
    trial_conditions = cell(1, ntrials);
    for trial_i = 1:ntrials
        trial_conditions{trial_i} = sprintf('%s-%02d', categories{categoryLabels(trial_i)}, ...
            exemplarLabels(trial_i));
    end

    % Add trials to fieldtrip structure
    for i = 1:ntrials
        % Store
        ft.trial{i} = permute(squeeze(X_3D(:, :, i)), [1, 2]);
        ft.time{i}  = t;
    end
    
    % Sort trials by condition for visualisation purposes
    [trial_conditions,sortorder] = sort(trial_conditions);
    ft.trial = ft.trial(sortorder);
    ft.time  = ft.time(sortorder);
    
    cd(out_dir);    

    % Convert to SPM format and pre-process
    % Output: CfdmD_**.mat
    % ---------------------------------------------------------------------           
    name = ['D_' subjects{s}];    
    D = apply_preprocessing(ft,name,trial_conditions);
    save(D);

    % Reduce dimensionality
    nmodes = 7;
    S = struct();
    S.D = D;              
    S.method = 'pca';
    S.settings.ncomp = nmodes;     
    S.settings.threshold = 0;
    D = spm_eeg_reduce(S);
    D = chantype(D, ':', 'LFP'); 
    save(D);    
    
    % Flip the sign of the ERPs if needed (PCA loses the sign)
    % pc = 1;
    % condition = 1;    
    % yy = squeeze(D(pc,:,condition));
    % max_pve = max(abs(yy(yy > 0)));
    % max_nve = max(abs(yy(yy < 0)));
    % if max_pve < max_nve
    %     flip = '*';
    %     D(:,:,:) = -D(:,:,:);
    %     save(D);
    % else
    %     flip = '';
    % end        
            
    % Plot timeseries
    subplot(3,4,s);
    yy = squeeze(D(:,:,1));
    plot(D.time,yy,'LineWidth',2);
    hold on
    yy = squeeze(D(:,:,2));
    plot(D.time,yy,'LineWidth',2);        
    set(gca,'FontSize',12);
    title(sprintf('Subject %s%s',subjects{s}),'FontSize',12);
    set(gca,'YTickLabel',[]);
    xline(0); yline(0); axis square;   
    xlabel('Time (s)');
    ylabel('ERP');
    cd(start_dir);
    
end

% Save figure
saveas(f,fullfile(out_dir,'ERPs.png'));

% Delete intermediate files
delete(fullfile(out_dir,'D_*'));
delete(fullfile(out_dir,'mD_*'));
delete(fullfile(out_dir,'dmD_*'));
delete(fullfile(out_dir,'fdmD_*'));
for s = 1:length(subjects)
    delete(fullfile(out_dir,['CfdmD_' subjects{s} '_2conditions*']));
end
    
% -------------------------------------------------------------------------
function D = apply_preprocessing(ft,name,condition_names)
    % Convert to SPM
    D = spm_eeg_ft2spm(ft, name);

    % Set channel types
    D = chantype(D, ':', 'LFP'); 

    % Set condition labels
    D = conditions(D, 1:(length(ft.trial)), condition_names);

    % Average over trials within each condition
    S = struct();
    S.D = D;               
    S.robust = true;
    D = spm_eeg_average(S); 