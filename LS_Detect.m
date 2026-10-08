%% DETECT slow waves

path_setting_nandi;
adult_all = subdir_adults_ALL;
child_all = subdir_child_ALL;

selected_ppts = [10];
subdir = child_all(selected_ppts);

detectwaves_all = cell(1,numel(subdir));

for ind = 1:numel(subdir)
    cd(subdir{ind});

    files = dir('*_preproc_noLP_task.mat');
    if isempty(files)
        error('No matching .mat file found in %s', subdir{ind});
    end
    load(files(1).name);

% if including LP filt, make sure this like is correctly named first
mags_data_cleaned = data_clean_mags;

trialdata_mags = mags_data_cleaned.trial{1};
[twa_results_mags]= MEG_LW_twalldetectnew_TA_v3(trialdata_mags,hdr.Fs,0);  % this '0' 3rd parameter tweak later to find best threshold
all_Waves_mags=[];
for nE_mags=1:size(trialdata_mags,1)
    all_Waves_mags=[all_Waves_mags ; [repmat([1 0 nE_mags],length(abs(cell2mat(twa_results_mags.channels(nE_mags).maxnegpkamp))),1) abs(cell2mat(twa_results_mags.channels(nE_mags).maxnegpkamp))'+abs(cell2mat(twa_results_mags.channels(nE_mags).maxpospkamp))' ...
        cell2mat(twa_results_mags.channels(nE_mags).negzx)' ...  % column 5
        cell2mat(twa_results_mags.channels(nE_mags).poszx)' ...
        cell2mat(twa_results_mags.channels(nE_mags).wvend)' ...
        cell2mat(twa_results_mags.channels(nE_mags).maxnegpk)' ...
        cell2mat(twa_results_mags.channels(nE_mags).maxnegpkamp)' ...
        cell2mat(twa_results_mags.channels(nE_mags).maxpospk)' ...
        cell2mat(twa_results_mags.channels(nE_mags).maxpospkamp)' ...
        cell2mat(twa_results_mags.channels(nE_mags).mxdnslp)' ...
        cell2mat(twa_results_mags.channels(nE_mags).mxupslp)' ...
        cell2mat(twa_results_mags.channels(nE_mags).maxampwn)' ...
        cell2mat(twa_results_mags.channels(nE_mags).minampwn)' ...
        ]];
    % Columns of the all_Waves matrix
    %  1: subject number (as in the for loop), you could replace with a subject ID if numerical
    %  2: 0 here but could be a block number
    %  3: electrode number (same order as the data)
    %  4: peak-to-peak amplitude
    %  5: start (in sample)
    %  6: half wave position
    %  7: end
    %  8: position of negative peak
    %  9: amplitude of negative peak
    % 10: position of positive peak
    % 11: amplitude of positive peak
    % 12: maximum downward slope
    % 13: maximum upward slope
    % 14: max amplitude on the slow wave window
    % 15: min amplitude on the slow wave window
end

save([subjectID '_all_detected_waves_task.mat'], 'subjectID', 'all_Waves_mags','trialdata_mags', 'hdr');

end

