%######################################################

function numTcells = getModelPrediction(...
    x, innate_params, scheme, data_index)
%---------------------------------------
% Inputs: Optimized kinetic parameters, Fixed kinetic parameters,
% struct of dosing schemes, and indices of the schemes we are interested
% Outputs: Number of T cells, Total DCs
%---------------------------------------

%innate_params = [innate_params, x']; % parameters
% Get model predictions

numTcells = zeros(length(data_index),1);   % for storing model prediction

for j=1:length(data_index)
    i = data_index(j);
    %innate_params_temp =  [innate_params(plotIndex,:), x'];
    [time,conc] = getInnateDynamics( ...
        scheme.T{i}, scheme.numshot{i}, scheme.k{i}, innate_params);
    numTcells(j) = conc(2100, end); % Take day 21 T cell number
end

end