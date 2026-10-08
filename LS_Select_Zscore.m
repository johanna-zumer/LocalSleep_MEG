
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% this script is adjustable (frequency / channel etc) - make sure all the
%%% variables are what you indend to look at before running!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path_setting_nandi;

subdir_adults = subdir_adults_ALL();
subdir_child = subdir_child_ALL();


selected_ppts = [3];
subdir = subdir_adults(selected_ppts);
%subdir = subdir_child(selected_ppts);

selected_waves_delta_all = cell(1,numel(subdir));

for ind = 1:numel(subdir)
    cd(subdir{ind});

    clear all_Waves_mags*
    files_detect = dir('*_all_detected_waves_task.mat');
    if isempty(files_detect)
        warning('No matching .mat file found in %s', subdir{ind});
        continue
    end
    load(files_detect(1).name);

    files_preproc = dir('*_preproc_noLP_task.mat');
    if isempty(files_preproc)
        warning('No matching .mat file found in %s', subdir{ind});
        continue
    end
    load(files_preproc(1).name);

    files_maxfilt_for_grads = dir('*_task_mc.fif');
    if isempty(files_maxfilt_for_grads)
        warning('No matching .fif file found in %s', subdir{ind});
        continue
    end
    maxfilt_file = fullfile(subdir{ind}, files_maxfilt_for_grads(1).name);

    %% params

    %%%%%%%%%%%%%%%%% CHANGE FREQ HERE FOR SELECTION %%%%%%%%%%%%%%%%
    paramSW.LimFrqW=[0.8 4]; % frequency range: delta [0.8 4] or theta [4.1 7]
    paramSW.max_Freq=4; % max frequency considered
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    paramSW.prticle_Thr=90; % percentile of slow waves selected
    paramSW.AmpCriterionIdx=4; % Amplitude criterion: 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=[]; % to select an absolute and not relative threshold
    paramSW.art_ampl=2e-12; % threshold to discard too high amplitude events - arbitrary
    paramSW.max_posampl=1e-12; % threshold for artifact on positive peak (blinks typically) - arbitrary

    all_Waves_mags=double(all_Waves_mags); % load in all_waves_mags and create a double
    all_freq_mags=1./(abs((all_Waves_mags(:,5)-all_Waves_mags(:,7)))./hdr.Fs); % select start and end of each wave and 'right array divide' by sample frequency to obtain frequency of each wave
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq_mags<paramSW.LimFrqW(1) | all_freq_mags>paramSW.LimFrqW(2) | all_freq_mags>paramSW.max_Freq)*100) % present the percentage of waves that are outside of the 0.8-4 freq window
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves_mags(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100) % present % of waves with P2P amplitude that exceeds the artefact threshold (paramSW.art_ampl)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves_mags(:,11)>paramSW.max_posampl | all_Waves_mags(:,14)>paramSW.art_ampl| abs(all_Waves_mags(:,15))>paramSW.art_ampl)*100) % present % of waves with positive peaks that exceed maxpos amplitude threshold
    if exist('all_Waves_mags_orig')
        all_Waves_mags=all_Waves_mags_orig; %prticle_Thr  % unsure what's happening here? Just stating prticle_Thr?
    end

    % now actually reject the above
    % keep 'orig' list aside;
    all_Waves_mags_orig=all_Waves_mags;
    % first reject those based on frequency
    all_Waves_mags(all_freq_mags<paramSW.LimFrqW(1) | all_freq_mags>paramSW.LimFrqW(2) | all_freq_mags>paramSW.max_Freq, :) =[];
    % then AmpCrit (often =4 = P2P amplitude)
    all_Waves_mags(all_Waves_mags(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl,:)=[];
    % then 11 = 'amplitude of positive peak', 14 = 'max amplitude in the window' and 15 = 'min amplitude in the window'
    all_Waves_mags(all_Waves_mags(:,11)>paramSW.max_posampl| all_Waves_mags(:,14)>paramSW.art_ampl| abs(all_Waves_mags(:,15))>paramSW.art_ampl,:)=[];

    thr_Wave_mags=[]; % create empty struct
    slow_Waves_mags=[]; % create empty struct
    ERP_SW_mags=[]; % create empty struct
    total_thresh=prctile(all_Waves_mags(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr); %%%%%%%%%%%%%%%%%%% unused... ?


    %% z score for frontal channel

    % do z-value of data just in the one channel that you are using below.
    % then find the value (fT or T) that corresponds to the 3*std, then set
    % that to be paramSW.fixThr

    %%%% Choose which channel you would like to use
    % FRONTAL = MEG0621
    % Central = MEG0731
    % Central (left leaning) = MEG0741 % exploratory
    % Posterior = MEG2121 (adults) or MEG2111 (child)  - lowest sensor may be
    % on child neck so used sensor above

    channel = match_str(data_clean_mags.label,'MEG0621'); % select individual channel
    thisE_Waves_mags_chan = all_Waves_mags(all_Waves_mags(:,3)==channel,:); % extract the data for this specific channel from the detected waves file

    % Extract the amplitude values (P2P) for this channel
    channel_amplitudes = thisE_Waves_mags_chan(:,paramSW.AmpCriterionIdx); % for the detected waves in this electrode, extract the data on P2P amplitudes

    % filter original timeseries data before using to clean up plots
    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 20; % just for the LS events
    orig_timeseries_LP = ft_preprocessing(cfg,data_clean_mags);

    % extract time series for this channel
    cfg = [];
    cfg.channel = channel;
    channel_data = ft_selectdata(cfg, orig_timeseries_LP); % extract the original time series from the preprocessed data file
    channel_timeseries = channel_data.trial{1}(1,:); % time series will just be one trial since it's resting state

    % Calculate z-scores for the time series data
    mean_amplitude = mean(channel_timeseries); % mean of preprocessed time series - should be 'zero', usually around e-18
    std_amplitude = std(channel_timeseries);

    % Calculate z-scores
    z_scores = (channel_timeseries - mean_amplitude) / std_amplitude;

    % Find the amplitude value that corresponds to 3 standard deviations (3SD
    % above mean)
    threshold_z = 3;
    amplitude_threshold = mean_amplitude + (threshold_z * std_amplitude);

    % Set this as the fixed threshold
    paramSW.fixThr = amplitude_threshold;

    total_waves_channel = length(channel_amplitudes); % measured by p2p amplitudes
    waves_above_threshold = sum(channel_amplitudes > amplitude_threshold); %number of waves with p2p amps above the 3sd


    %% sensor

    channel_counter=1;
    for nE_megmags=channel % select the 'electrode'
        thisE_Waves_mags=all_Waves_mags(all_Waves_mags(:,3)==nE_megmags,:); % for this electrode, select all waves
        temp_p2p_mags=thisE_Waves_mags(:,paramSW.AmpCriterionIdx); % p2p amplitudes of all detected slow waves

        % original way / now using total_thresh over all channels
        if ~isempty(paramSW.fixThr)   % with z-score, this will be set to 3*std, and not be empty
            thr_Wave_mags(nE_megmags)=paramSW.fixThr; %
        else
            % this is a channel-specific percentile threshold - not used
            % when z-scoring
            thr_Wave_mags(nE_megmags)=prctile(thisE_Waves_mags(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
        end
        slow_Waves_mags=[slow_Waves_mags ; thisE_Waves_mags(temp_p2p_mags>thr_Wave_mags(nE_megmags),:)]; % for this electrode, select the p2p amplitudes that are greater than the set threshold (thr_waves)

        temp_ERP_mags=[]; % empty struct
        temp_Waves_mags=thisE_Waves_mags(temp_p2p_mags>thr_Wave_mags(nE_megmags),:); % for all the waves from this electrode, select waves with amplitudes above the threshold
        for nW_mags=1:size(temp_Waves_mags,1) % for loop over selected waves, from every channel
            wave_onset_mags=temp_Waves_mags(nW_mags,5); % extract the start time of the wave
            if min(wave_onset_mags+(-1*hdr.Fs:1*hdr.Fs))<1 || max(wave_onset_mags+(-1*hdr.Fs:1*hdr.Fs))>size(trialdata_mags,2)
                disp(['Wave onset out of bounds for electrode' num2str(nE_megmags)]);
                continue;
            end % decide if wave onset is out of bounds for detection
            vec_EEG_mags=trialdata_mags(nE_megmags,wave_onset_mags+(-1*hdr.Fs:1*hdr.Fs)); % from the detected waves (trialdata_mags), for this channel, look at -1*sample frequency to 1*sample frequency ??
            vec_EEG_mags=vec_EEG_mags-mean(vec_EEG_mags(1:hdr.Fs/2));  % baseline subtraction
            if max(abs(vec_EEG_mags))>150 % if the max SOMETHING? is over 150, continue.
                continue;
            end
            temp_ERP_mags=[temp_ERP_mags ; vec_EEG_mags]; % temp_ERP_mags is whatever os going on with this threshold??
        end

        if size(temp_ERP_mags,1)>5
            ERP_SW_mags(channel_counter,:)=mean(temp_ERP_mags,1);
        else
            ERP_SW_mags(channel_counter,:)=nan(1,length((-1*hdr.Fs:1*hdr.Fs)));
        end
        channel_counter=channel_counter+1;
    end

    fprintf('Channel-specific z-score threshold calculation:\n');
    fprintf('Total waves detected in channel: %d\n', total_waves_channel);
    fprintf('Waves above threshold (waves_above_threshold): %d\n', waves_above_threshold);
    fprintf('Percentage of waves above threshold: %.1f%%\n', (waves_above_threshold / total_waves_channel) * 100);
    fprintf('Mean amplitude: %.2e\n', mean_amplitude);
    fprintf('Std amplitude: %.2e\n', std_amplitude);
    fprintf('Z-score threshold: %.1f\n', threshold_z);
    fprintf('Amplitude threshold (%.1f std above mean): %.2e\n', threshold_z, amplitude_threshold);

    fprintf('Waves above threshold (nW_mags): %d\n', nW_mags);

    % for plotting
    cfg=[];
    cfg.layout='neuromag306mag.lay';
    layout=ft_prepare_layout(cfg);
    numchan_mags=length(data_clean_mags.label);


    %%  JZ trigger marker adding

    artifact=nan(size(slow_Waves_mags,1),2);
    artifact=[slow_Waves_mags(:,5) slow_Waves_mags(:,7)];

    % set time window - ORIGINAL (MAGS)
    cfg=[];
    cfg.trl=[artifact(:,1)-1500 artifact(:,1)+2000 -1500*ones(size(artifact,1),1)];
    data_LS_longtrl=ft_redefinetrial(cfg,data_clean_mags); % can filter this again to 0.5 to 40Hz now to clean up plots
    
    if any(any(isnan(data_LS_longtrl.trial{1})))
        cfg=[];
        cfg.trials=2:length(data_LS_longtrl.trial);
        data_LS_longtrl=ft_selectdata(cfg,data_LS_longtrl);
    end
    if any(any(isnan(data_LS_longtrl.trial{length(data_LS_longtrl.trial)})))
        cfg=[];
        cfg.trials=1:(length(data_LS_longtrl.trial)-1);
        data_LS_longtrl=ft_selectdata(cfg,data_LS_longtrl);
    end
    
    
    %%%%%%%%%%%%%%%%  NEW  - grads  %%%%%%%%%%%%%%

    cfg = [];
    cfg.dataset = maxfilt_file;
    cfg.channel = 'MEG'; % all chans
    cfg.continuous = 'yes';
    cfg.demean = 'yes';
    cfg.hpfilter = 'yes';
    cfg.hpfiltord = 3;
    cfg.hpfreq = 0.5;
    hdr_grads = ft_read_header(cfg.dataset);
    % IF NEEDED: skip first second (adjust if needed) - if trouble loading in
    % cfg.trl = [hdr_grads.Fs*5, hdr.nSamples, 0]; % change here as needed
    data_initial = ft_preprocessing(cfg);

    cfg=[];
    cfg.channel = ft_channelselection({'*2', '*3'}, data_initial.label); % grads end in 2 or 3
    data_grads = ft_selectdata(cfg, data_initial);
   
    cfg=[];
    cfg.trl=[artifact(:,1)-1500 artifact(:,1)+2000 -1500*ones(size(artifact,1),1)];
    data_LS_longtrl_grads=ft_redefinetrial(cfg,data_grads);

    if any(any(isnan(data_LS_longtrl_grads.trial{1})))
        cfg=[];
        cfg.trials=2:length(data_LS_longtrl_grads.trial);
        data_LS_longtrl_grads=ft_selectdata(cfg,data_LS_longtrl_grads);
    end
    if any(any(isnan(data_LS_longtrl_grads.trial{length(data_LS_longtrl_grads.trial)})))
        cfg=[];
        cfg.trials=1:(length(data_LS_longtrl_grads.trial)-1);
        data_LS_longtrl_grads=ft_selectdata(cfg,data_LS_longtrl_grads);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subjectID = str2double(subjectID);
    slow_Waves_mags(:,1) = subjectID;

    % check variables first, then save if happy

    %save([num2str(subjectID) '_selected_waves_delta_frontal_task_withGrads.mat'], 'subjectID','data_LS_longtrl','data_LS_longtrl_grads','slow_Waves_mags','paramSW','total_waves_channel','waves_above_threshold','amplitude_threshold', 'mean_amplitude', 'std_amplitude', 'threshold_z', 'channel', '-v7.3');
    
end


%% BELOW IS JUST TO CHECK INDIVIDUAL PARTICIPANTS - exploratory / sanity checking

% % check individual wave events
 % cfg=[];
 % cfg.channel = channel;
 % cfg.lpfilter = 'yes'; % to clean up plots
 % cfg.lpfreq = 10;
 % ft_databrowser(cfg,data_LS_longtrl)
 % title('frontal sensor')
 
%%%%% plots - ORIG

 cfg=[];
 LS_avg_long=ft_timelockanalysis(cfg,data_LS_longtrl);
% figure;plot(LS_avg_long.time,LS_avg_long.avg)
% title('Avg LS events - from frontal sensor')

 channel_mag = match_str(data_clean_mags.label,'MEG0621');

% check avg of events
  cfg=[];
  cfg.channel = channel_mag;
  cfg.lpfilter = 'yes';
  cfg.lpfreq = 10;
  cfg.linecolor = [0 0 0];
  cfg.linewidth = 4;
  ft_databrowser(cfg,LS_avg_long)
  title('frontal sensor')

 % topo
 % cfg=[];
 % cfg.layout = 'neuromag306mag.lay';
 % cfg.parameter = 'avg';
 % cfg.zlim = 'maxabs';
 % ft_topoplotER(cfg,LS_avg_long)
 % title('LS events - from frontal sensor')

%%%%% NEW

% average grad LS events
%  cfg=[];
% %cfg.preproc.lpfilter = 'yes'; %ADDED
% %cfg.preproc.lpfreq = 10; %ADDED
% %cfg.preproc.demean = 'yes'; %ADDED
%  LS_avg_long_grads=ft_timelockanalysis(cfg,data_LS_longtrl_grads);
%  figure;plot(LS_avg_long_grads.time,LS_avg_long_grads.avg)
%  title('Avg LS events GRADS')

 % %%% ADDED - combine grads
 % cfg =[];
 % LS_avg_long_combined = ft_combineplanar(cfg,LS_avg_long_grads);

% channel_grad = match_str(data_grads.label,'MEG2042');

% avg time series - SEP
  % cfg=[];
  % cfg.channel = channel_grad;
  % cfg.lpfilter = 'yes';
  % cfg.lpfreq = 10;
  % cfg.linecolor = [0 0 0];
  % cfg.xlim = [-0.5 1];
  % ft_databrowser(cfg,LS_avg_long_grads)
  % title('frontal sensor')

 % topo - sep
 % cfg=[];
 % cfg.layout = 'neuromag306planar.lay';
 % cfg.parameter = 'avg';
 % cfg.zlim = 'maxabs';
 % ft_topoplotER(cfg,LS_avg_long_grads)
 % title('LS events - from frontal sensor')
     
 % time series - cmb
 % cfg=[];
 % cfg.channel = channel_grad;
 % cfg.lpfilter = 'yes';
 % cfg.lpfreq = 10;
 % cfg.linecolor = [0 0 0];
 % cfg.xlim = [-0.5 1];
 % ft_databrowser(cfg,LS_avg_long_combined)
 % title('frontal sensor')

 % topo - cmb
 % cfg = [];
 % cfg.layout = 'neuromag306cmb.lay';
 % cfg.parameter = 'avg';
 % cfg.zlim = 'maxabs';
 % ft_topoplotER(cfg,LS_avg_long_combined)
 % title('LS events - from frontal sensor')
 %%%%

% freq analysis
% wavelet
 % cfg = [];
 % cfg.method     = 'wavelet';  %mtmconvol - do this one for higher freqs
 % cfg.width      = 2; % dont go lower than 2
 % cfg.output     = 'pow';
 % cfg.foi        = 1:1:20; % -- doesn't include gamma!
 % cfg.toi        = -1.5:0.05:2;
 % TFRwave = ft_freqanalysis(cfg, data_LS_longtrl);

 % cfg = [];
 % cfg.baseline     = [-1.0 -0.5];
 % cfg.baselinetype = 'relative';
 % cfg.masknans     = 'yes';
 % cfg.maskstyle    = 'saturation';
 % cfg.layout       = 'neuromag306mag.lay';
 % cfg.zlim         = 'maxabs';
 % cfg.colorbar     = 'yes';
 
% Plot
 % figure;
 % ft_multiplotTFR(cfg, TFRwave);
 % title('Avg local sleep event (I think avg?)');

% % mtmconvol
% % freq dependent window
 % cfg              = [];
 % cfg.output       = 'pow';
 % cfg.method       = 'mtmconvol';
 % cfg.taper        = 'hanning';
 % cfg.foi          = 1:1:20; % doesn't include gamma!
 % cfg.t_ftimwin    = 2./cfg.foi;  % 2 cycles per time window
 % cfg.toi          = -1.5:0.05:2;
 % TFRmtm = ft_freqanalysis(cfg, data_LS_longtrl);
 
% % freq dependent - visuals
 % cfg              = [];
 % cfg.baseline     = [-1.0 -0.5];
 % cfg.baselinetype = 'absolute';
 % cfg.maskstyle    = 'saturation';
 % cfg.zlim         = 'maxabs';
 % cfg.interactive  = 'yes';
 % cfg.layout       = 'neuromag306mag.lay';
 % figure
 % ft_singleplotTFR(cfg, TFRmtm);
 
