

% LW and JZ notes for local sleep detection in MEG

%%%%%
% preprocessing overview
%%%%%%

% 1 - maxfilter - used on Elekta MEGIN data
% usually SSS is enough, used tSSS where needed (normally for children with excessive movement) - LSW detection works fine on both.

% 2 - basic preprocessing
% demeaned
% hpfilter order = 3
% hp freq = 0.5

% then split into magnetometers and non-combined gradiometers

% run fastICA on mags only, rejecting cardiac and occular artefacts -
% normally around 4 components, tried to avoid over cleaning

% saved out using --- save([subjectID '_preproc_noLP_task.mat'], 'subjectID', 'data_clean_mags', 'hdr');
% where 'data_clean_mags' is the output variable of the component
% rejection and 'hdr' is simply extracted from the original maxfiltered dataset using ft_read_header.

% We created the gradiometer versions of the waves later, but could try preprocessing the grads and running the detection on them directly.
% Initial attempts of running the LSW detection on gradiometers gave confusing results so we did not pursue it.

%%%%%
% LS_Detect.m
% detection of LS waves (across all magnetometer sensors)
%%%%%

% LW split the original algorithm (from Andrillon) into 'detection' and 'selection' scripts to help manage keeping track of each aspect.


% 'path_setting_nandi' loads in...(set this specific to your own computing environment) 
% ---- fieldtrip version 'fieldtrip-20240916' - could be fine with earlier or later versions, but definetely works with this.
% ----'LSCPtools' found here - 'https://github.com/andrillon/LSCPtools' -
% needed for running detection, we also used v3 of
% the twalldetectnew_TA script for our adjustments (found in the sleeptools subfolder)

% the main detection script is largely unchanged, we made some slight adjustments to
% the twalldetect function (see 'MEG_LW_twalldetectnew_TA_v3.m')

%%%%%
% LS_Select_Zscore.m
% Selection of LS waves based on threhold params and from specific sensors etc
%%%%%

% see comments throughout script

% adjusted 90th percentile threshold to use individualised Z-scoring
% instead - Wienke et al (2021) have done something similar.

% also created wave-based epochs / pseudo trials around each individual wave

% used the wave-based epoch sample points to create parallel gradiometer
% data epochs -- see previous note about alternatively running detection on
% gradiometer data directly


% from here, saved out selected waves can then be used (alongside the preprocessed data file) for further task-based or
% rest-state analyses.



