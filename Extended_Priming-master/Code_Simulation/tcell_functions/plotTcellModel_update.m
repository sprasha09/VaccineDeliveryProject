% This code is used to quickly obtain the plots of innate immune cells and
% Tfh cells using pre-defined parameters. 

% Pre-defined kinetic parameters
x_optim = dlmread('optimized_parameters.txt', '\t');
% Schemes to compare
scheme = defineSchemes();

% Load experimental data
Tfh_data = readmatrix(fullfile(fileparts(mfilename('fullpath')),...
  '../../Data/Tfh_Number_Data.xlsx'), 'Sheet', 'Sheet1', 'Range', 'B1:F7');            

% Plotting
% fit_parameter = 0;

% innate_inputs = [3, 3, 0.1, 1.5, 0.22, 1.2, 10];

innate_inputs = [3, 3, 0.1, 1.5, 0.22, 1.2, 0.1;
                  0.1, 3, 0.1, 1.5, 0.22, 1.2, 0.1;
                  10, 3, 0.1, 1.5, 0.22, 1.2, 0.1];

% [x_optim, chisq, params_guess, res, exitflags] =...
%     tCellModel(fit_parameter, x_optim, innate_inputs);

% Plotting
f1 = plotInnateDyn_update(x_optim, innate_inputs, scheme, 1); % Plot the innate immune cell dynamics
%f2 = plotTcellNum_update(x_optim, innate_inputs, scheme, Tfh_data, 1); % Compare the model with data

% innate_inputs = [3, 3, 0.1, 1.5, 0.22, 1.2, 10];
%d_Ag, d_Adj, S_Adj, alpha, delta, mu, k
% 
% [x_optim, chisq, params_guess, res, exitflags] =...
%     tCellModel(fit_parameter, x_optim, innate_inputs);
