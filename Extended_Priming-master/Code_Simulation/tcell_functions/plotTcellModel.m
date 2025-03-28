% This code is used to quickly obtain the plots of innate immune cells and
% Tfh cells using pre-defined parameters. 

% Pre-defined kinetic parameters
x_optim = dlmread('optimized_parameters.txt', '\t');

% Plotting
fit_parameter = 0;
[x_optim, chisq, params_guess, res, exitflags] =...
    tCellModel(fit_parameter, x_optim);