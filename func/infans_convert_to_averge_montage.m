function eegout = infans_convert_to_averge_montage(eegin)
% This function converts the original montage to common average
%
% INPUT: 
%   eegin  : a matrix contains EEG data (channels X samples)
%
% OUTPUT:
%   eegout : the changed version of eegin in term of EEG montage (channels X samples)

    eegout = double(eegin - mean(eegin));
end