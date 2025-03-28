function numTcells = getNumTcells(param)
%--------------------------------
% Calculates the number of T cells based on the given parameters.
%--------------------------------

% Pre-defined kinetic parameters
x_optim = dlmread('optimized_parameters.txt', '\t');
kinetic_params = [3, 3, 0.1, 1.5, 0.22, 1.2, 10, x_optim(1), x_optim(2)];

% Define T_slow based on the pSER parameter value.
%   If pSER is 0, there is no extended antigen release.
%   If pSER is 1, only the last dose is extended antigen release.
%   If pSER is 2, final two doses are extended antigen release. 
if param.pSER == 0
    T_slow = [];
else
    T_slow = repmat([param.pSER_ag_release_time, param.pSER_adj_release_time], ...
        param.pSER,1);
end

% Retrieve time and concentration values
[time,conc] = getInnateDynamics(param.T*[1,1], ...
                                param.numshot*[1,1], ...
                                param.k*[1,1], ...
                                kinetic_params, ...
                                T_slow);

% Calculate the number of T cells
numTcells = conc(:,7:8)*10; % Scaling factor of 10 is applied
end