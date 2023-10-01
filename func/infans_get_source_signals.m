function Source = infans_get_source_signals(Sensor, InverseOperator, CollapseOperator, MyAtlas)
% This function reconstructs source activities.
%
% see Tokariev et al., 2019, Cerebral Cortex
%
% INPUTS:
%  Sensor           : filtered EEG data (samples X channels)
%  InverseOperator  : Inverse solution for infant head with N = 19 electrodes
%  CollapseOperator : weights for source signals within 'host' parcels
%  MyAtlas          : assignment of src/verticies to parcels (in MyAtlas.Parcels)
%
% OUTPUTS:
%  Source           : filtered parcel signals (samples X parcels (58))
 
    Np = length(MyAtlas.Parcels);  % number of parcels
    L  = size(Sensor, 1);          % EEG length 
    CollapseOperator = repmat(CollapseOperator, 1, L);

    % source signals  
    src_buf = (InverseOperator * Sensor') .* CollapseOperator;      % (sources X samples)

    % parcel/cortical signals  
    parcel_buf = zeros(Np, L);

    for j = 1:Np
        parcel_buf(j, :) = mean(src_buf(MyAtlas.Parcels{j, 1}, :)); % (parcels X samples) 
    end

    Source = parcel_buf'; % (samples X parcels)   
end

function [data,M] = infans_leakcorr(data,T,order)
% Correction for volume conduction, when working with source space M/EEG
% If order > 0, it performs Pasqual-Marqui's method (order refers to the MAR order)
%   Pascual-Marqui et al. (2017)
% If order < 1, it performs symmetric orthogonalisation as specified in 
%   Colclough et al. (2015)
% Output M is only defined for Pasqual-Marqui's method

if order > 0 
    if isstruct(data)
        [XX,Y] = formautoregr(data.X,T,1:order,order,1);
    else
        [XX,Y] = formautoregr(data,T,1:order,order,1);
    end
    B = XX \ Y;
    eta = Y - XX * B;
    epsilon = closest_orthogonal_matrix(eta);
    M = pinv(epsilon) * eta;
    if isstruct(data)
        data.X = data.X * inv(M);
    else
        data = data * inv(M);
    end
else
    if isstruct(data)
        data.X = closest_orthogonal_matrix(data.X);
    else
        data = closest_orthogonal_matrix(data);
    end
    M = [];
end

end