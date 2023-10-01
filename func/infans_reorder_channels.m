function eegout = infans_reorder_channels(eegin)
% This function reorders the channels to make them in line with the
% inverse operator proposed by Anton (see the tutorial file in the following link):
% https://github.com/babyEEG/AED-exposed-infants
% Channels for the inverse operator:
% Fp1, Fp2, F7, F3, Fz, F4, F8, T3, C3, Cz, C4, T4, T5, P3, Pz, P4, T6, O1, O2
%
% The order of channels in the sleep dataset is as follows:
% Fp1, Fp2, F3, F4, C3, C4, P3, P4, O1, O2, F7, F8, T3, T4, T5, T6, Fz, Cz, Pz 
%
%
% INPUT:
%   eegin  : a matrix contains EEG data (channels X samples)
%
% OUTPUT:
%   eegout : the reordered version of eegin (channels X samples)

    inverseOrder = ["Fp1"; "Fp2"; "F7"; "F3"; "Fz"; ...
                    "F4"; "F8"; "T3"; "C3"; "Cz";   ...
                    "C4"; "T4"; "T5"; "P3"; "Pz";   ...
                    "P4"; "T6"; "O1"; "O2"];

    mainOrder    = ["Fp1"; "Fp2"; "F3"; "F4"; "C3"; ...
                    "C4"; "P3"; "P4"; "O1"; "O2"; ...
                    "F7"; "F8"; "T3"; "T4"; "T5"; ...
                    "T6"; "Fz"; "Cz"; "Pz"];

    eegout = zeros(size(eegin));           
    for i = 1:length(inverseOrder)
        index = find(inverseOrder == mainOrder(i));
        eegout(index,:) = eegin(i,:);
    end
    eegout = double(eegout);
end