%% description
% This script reads the data, filters and transfomrs them into 
% the source space by the model provided by Anton Tokariev 
% (https://github.com/babyEEG/AED-exposed-infants).
%
% NOTE: Some lines must be changed based on your data path.

%% initialization
close all; clear; clc;

subjectList   = setdiff(1:68, [10 19 24:25 45 59:60 66]); % excludes 8 neonates
dataPath      = '.\Data\Healthy\';
savePath      = '.\Data\';
headModelPath = '.\head model\';

% loading source space stuff
load([headModelPath 'Atlas.mat'])
load([headModelPath 'CollapseOperator.mat'])
load([headModelPath 'InverseOperator.mat'])

% leakage correction and dipole resolving initalization
Fs    = 100;           % sampling rate
[p,q] = rat(Fs / 250); % for resampling

%% designing band pass filter (0.3-25 Hz)
lpButter = designfilt('lowpassiir', 'FilterOrder', 7, ...
'HalfPowerFrequency', 22, 'SampleRate', 250, 'DesignMethod', 'butter');     
hpButter = designfilt('highpassiir', 'FilterOrder', 7, ...
'HalfPowerFrequency', 0.4, 'SampleRate', 250, 'DesignMethod', 'butter');

%% transforming the sensor-level signals into the source-level signals
mat_files = cell(length(subjectList),1);
T_all     = cell(length(subjectList),1);

for dataID = 1:length(subjectList)
    tic
    % loads .SET files for each subject  
    EEG_H_AS = pop_loadset([dataPath 'AS_H_' num2str(subjectList(dataID)) '.set']);
    EEG_H_QS = pop_loadset([dataPath 'QS_H_' num2str(subjectList(dataID)) '.set']);
    
    % changes the montage to the common average montage
    % in line with the preprocessing step in
    % (Tokariev et. al, 2019 https://www.nature.com/articles/s41467-019-10467-8)
    EEG_H_AS.data = infans_convert_to_averge_montage(EEG_H_AS.data);
    EEG_H_QS.data = infans_convert_to_averge_montage(EEG_H_QS.data);
    
    % Reorders EEG channels based on the "Inverse Operator" proposed by Anton
    % The order of channels for that model should be as follows:
    % Fp1, Fp2, F7, F3, Fz, F4, F8, T3, C3, Cz, C4, T4, T5, P3, Pz, P4, T6, O1, O2
    eeg_H_AS = infans_reorder_channels(EEG_H_AS.data); % channels X samples
    eeg_H_QS = infans_reorder_channels(EEG_H_QS.data); % channels X samples
    
    clear EEG_H_AS EEG_H_QS

    % removes the DC part from all EEG signals (based on Anton's suggestion)
    eeg_H_AS = eeg_H_AS - mean(eeg_H_AS, 2); % channels X samples
    eeg_H_QS = eeg_H_QS - mean(eeg_H_QS, 2); % channels X samples
    
    % filters, resmaples, and selects 3 minutes
    filtered_H_AS = filtfilt(hpButter,filtfilt(lpButter,eeg_H_AS')); % samples X channels
    filtered_H_QS = filtfilt(hpButter,filtfilt(lpButter,eeg_H_QS')); % samples X channels
    
    filtered_H_AS = resample(filtered_H_AS,p,q); % chnages Fs from 250 to 100 Hz
    filtered_H_QS = resample(filtered_H_QS,p,q); % chnages Fs from 250 to 100 Hz
    
    filtered_H_AS = filtered_H_AS(1: Fs * 180, :); % selects just 3 minutes
    filtered_H_QS = filtered_H_QS(1: Fs * 180, :); % selects just 3 minutes
     
    clear eeg_H_AS eeg_H_QS
    
    % reconstructs source activities
    Source_H_AS = infans_get_source_signals(filtered_H_AS, InverseOperator, CollapseOperator, Atlas); % samples X parcels
    Source_H_QS = infans_get_source_signals(filtered_H_QS, InverseOperator, CollapseOperator, Atlas); % samples X parcels
    
    % removes the DC part from all souce activities 
    Source_H_AS = Source_H_AS - mean(Source_H_AS); % samples X parcels
    Source_H_QS = Source_H_QS - mean(Source_H_QS); % samples X parcels

    % concatenates AS and QS for each subject
    Source = vertcat(Source_H_AS, Source_H_QS);

    clear filtered_H_AS filtered_H_QS
    
	% makes input structures for tde-hmm algorithm
	% read the documentaion on the folloing link:
	% https://github.com/OHBA-analysis/HMM-MAR/wiki/User-Guide
	
    mat_files{dataID, 1} = Source;                                    % samples X channels
    T_all{dataID, 1}     = [size(Source_H_AS,1) size(Source_H_QS,1)]; % no. of samples for each epoch
    
    clear Source Source_H_AS Source_H_QS

    elapsed_time = toc;   
    fprintf(['Subject ' num2str(dataID) ' out of ' num2str(length(subjectList)) ', elapsed time: ' num2str(elapsed_time) '\n'])
end
save([savePath 'Data_Source.mat'], 'mat_files', 'T_all');
fprintf('Source reconstruction is done!\n')
fprintf('------------------------------\n')