function vRSA_demo_preprocess
% Imports EEG data from:
%
% Miller, K.J., Schalk, G., Hermes, D., Ojemann, J.G. and Rao, R.P., 2016. 
% Spontaneous decoding of the timing and content of human object perception 
% from cortical surface recordings reveals complementary information in the 
% event-related potential and broadband spectral change. 
% PLoS computational biology, 12(1), p.e1004660.
% 
% Data downloaded from:
% https://purl.stanford.edu/xd109qh3109
%
% Outputs:
%
% subjects/CfdmD_**.mat - channel-wise, stimulus-wise ERP data for vRSA
%
% subjects/RCfdmD_**_2conditions_1channel.mat - single ERP summaries of
% faces and houses, used for DCM analysis to derive basis functions
% 
% by Peter Zeidman
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
subjects = {'ca','de','fp','ja','mv','wc','zt'};
nsubjects = length(subjects);

% Convert data to Fieldtrip and then to SPM format
% -------------------------------------------------------------------------
start_dir = pwd;

% Which brain area IDs to retain (see raw_data/fhpred_master.m for key)
areas = 1:19;

f = figure;
for s = 1:nsubjects 
    disp(subjects{s});
    
    % Prepare subject's data in Fieldtrip format
    % ---------------------------------------------------------------------
    
    % Load data (times x channels)
    load(fullfile(raw_dir,'data',subjects{s},[subjects{s} '_faceshouses.mat']));    
    
    % Identify channels corresponding to the selected brain areas
    load(fullfile(raw_dir,'locs',[subjects{s} '_xslocs.mat']));
    to_retain = find(ismember(elcode,areas));       
    
    % Limit to those channels
    data = data(:,to_retain);
    [ntimes,nchannels] = size(data);
    
    % Set channel labels    
    labels = cell(nchannels,1);
    for i = 1:nchannels
        labels{i} = num2str(i);
    end    
    ft = struct();
    ft.label = labels;
    
    % Identify start of each 400ms stimulus in samples
    % We record each trial as starting 50 samples (50ms) before the 
    % stimulus onset and lasting until 200 samples (200ms) after the 
    % stimulus ends, thus 650ms total.
    trial_start = [];
    trial_end   = [];
    trial_conditions = {};
    last_stim   = 0;
    pre_stim_isi_samples  = 50;
    post_stim_isi_samples = 200;    
    for i = 1:ntimes
        if stim(i) ~= last_stim
            is_picture = stim(i) > 0 && stim(i) < 101;
            is_ISI     = stim(i) == 101;
            
            if is_ISI && (last_stim > 0 && last_stim < 101) 
                % We've transitioned to a pre-stimulus ITI from a stimulus            
                trial_end(end+1) = i-1+post_stim_isi_samples;
            elseif is_picture
                % We've transitioned to a stimulus
                trial_start(end+1) = i-pre_stim_isi_samples;
                if stim(i) <= 50
                    trial_conditions{end+1} = sprintf('House%02d',stim(i));
                else
                    trial_conditions{end+1} = sprintf('Face%02d',stim(i));
                end
            end
        end
        last_stim = stim(i);        
    end
    ntrials = length(trial_start);
    
    % Acquisition times relative to stimulus
    stim_dur_samples  = 400;    
    trial_dur_samples = pre_stim_isi_samples+stim_dur_samples+post_stim_isi_samples;
    timeonset_samples = -pre_stim_isi_samples;
    t = ((0:(trial_dur_samples-1)) + timeonset_samples)/srate;

    for i = 1:ntrials
        % Get data [channels x time]
        y = data(trial_start(i):trial_end(i),:)';
        
        % Baseline correct             
        bl = mean(y(:,1:pre_stim_isi_samples), 2);
        y  = y - bl;
        
        % Store
        ft.trial{i} = y;
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
       
    % Repeat with trials collapsed into face and houses
    % Output: CfdmD_**_2conditions.mat
    % ---------------------------------------------------------------------       
    
    % Re-label trials
    is_house = startsWith(trial_conditions,'House');
    cond = cell(1,ntrials);
    for i = 1:length(cond)
        if is_house(i), cond{i} = 'House'; else, cond{i} = 'Face'; end
    end    
    
    % Preprocess
    name = ['D_' subjects{s} '_2conditions'];
    D = apply_preprocessing(ft,name,cond);           
    
    save(D);    
    
    % Repeat with trials collapsed into face and houses and just 1 channel
    % Output: CfdmD_**_2conditions_1channel.mat
    % ---------------------------------------------------------------------       
    
    % Re-label trials
    is_house = startsWith(trial_conditions,'House');
    cond = cell(1,ntrials);
    for i = 1:length(cond)
        if is_house(i), cond{i} = 'House'; else, cond{i} = 'Face'; end
    end    
    
    % Preprocess
    name = ['D_' subjects{s} '_2conditions_1channel'];
    D = apply_preprocessing(ft,name,cond);
    
    % Reduce dimensionality
    nmodes = 1;
    S = struct();
    S.D = D;              
    S.method = 'pca';
    S.settings.ncomp = nmodes;     
    S.settings.threshold = 0;
    D = spm_eeg_reduce(S);
    D = chantype(D, ':', 'LFP'); 
    save(D);    
    
    % Flip the sign of the ERPs if needed (PCA loses the sign)
    pc = 1;
    condition = 1;    
    yy = squeeze(D(pc,:,condition));
    max_pve = max(abs(yy(yy > 0)));
    max_nve = max(abs(yy(yy < 0)));
    if max_pve < max_nve
        flip = '*';
        D(:,:,:) = -D(:,:,:);
        save(D);
    else
        flip = '';
    end        
            
    % Plot timeseries
    subplot(2,4,s);
    yy = squeeze(D(:,:,1));
    plot(D.time,yy,'LineWidth',2);
    hold on
    yy = squeeze(D(:,:,2));
    plot(D.time,yy,'LineWidth',2);        
    set(gca,'FontSize',12);
    title(sprintf('Subject %s%s',subjects{s},flip),'FontSize',12);
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

    % Downsample
    S = struct();
    S.D = D;
    S.fsample_new = 200;
    D = spm_eeg_downsample(S);

    % Low-pass filter (removes high frequencies)
    S = struct();
    S.D = D;
    S.band = 'low';
    S.freq = 30;
    D = spm_eeg_filter(S);     
            
    % Apply a Hanning window (like that used in DCM)
    S = struct();
    S.D = D;
    S.detrend = 0;
    S.hanning = true;
    S.chtype  = 'LFP';
    D = spm_eeg_erp_correction(S);