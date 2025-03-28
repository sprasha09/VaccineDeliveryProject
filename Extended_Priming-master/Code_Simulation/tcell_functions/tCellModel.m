% Author: Leerang Yang
% Last modified 2023/07/03
% Here, I use the ODE model of innate immune cells and T cells to study the
% dependence on the dosing kinetics. Unknown parameters are fitted to the
% data observation if desired. 

function [x_optim, chisq, params_guess, res, exitflags] = tCellModel(varargin)
%  Analyzes the ODE model of innate immune cells and T cells to 
%  study the dependence on the dosing kinetics. Fitting unknown parameters 
%  to data observation is possible.
%
% Inputs:
%   - varargin: Variable number of input arguments.
%               1st argument: fit_parameters (1 if parameters are to be fitted, 
%                                           0 otherwise)
%               2nd argument (optional): x_optim value
%               3rd argument (optional): innate_inputs 
%
% Outputs:
%   - x_optim: Optimal parameters.
%   - chisq: Chi-squared statistic.
%   - params_guess: Initial guess parameters.
%   - res: Result of optimization.
%   - exitflags: Exit flags of optimization.

% Parse Input Arguments
fit_parameters = varargin{1};
if nargin == 3
    if fit_parameters == 1
        error('fit_parameters must be 0 when x_optim value is given');
    end
    x_optim = varargin{2};
    INNATE_PARAMS = varargin{3};
elseif nargin > 3
    error('Incorrect number of input arguments.');
end

% Constants
% INNATE_PARAMS = [3, 3, 0.1, 1.5, 0.22, 1.2, 10]; % Kinetic parameters
                %d_Ag, d_Adj, S_Adj, alpha, delta, mu, k

% Schemes to compare
scheme = defineSchemes();
    %Schemes represented:
    % (1)7-ED, (2)6-ED, (3)4-ED, (4)3-ED, (5)2-ED(20_80, d0_d7), 
    % (6)Bolus, (7)7-ED+Adjuvant bolus

% Load experimental data
Tfh_data = readmatrix(fullfile(fileparts(mfilename('fullpath')),...
  '../../Data/Tfh_Number_Data.xlsx'), 'Sheet', 'Sheet1', 'Range', 'B1:F7');            

% Fitting parameters to data
if fit_parameters==1
    [x_optim, chisq, params_guess, res, exitflags] = fitData(INNATE_PARAMS, scheme, Tfh_data);
else
    chisq = [];
    params_guess = [];
    res = [];
    exitflags = [];
end

% Plotting
f1 = plotInnateDyn(x_optim, INNATE_PARAMS, scheme, [1,2,3,4,5,6]); % Plot the innate immune cell dynamics
f2 = plotTcellNum(x_optim, INNATE_PARAMS, scheme, Tfh_data, [1,2,3,4,5,6]); % Compare the model with data

% Plot 7-ED vs. 7-ED(Adjuvant Bolus)
% f1 = plotInnateDyn(x_optim, INNATE_PARAMS, scheme, [1,7]);
% f2 = plotTcellNum(x_optim, INNATE_PARAMS, scheme, Tfh_data, [1,7]);

%% Save results
% Define the output directory
outputDir = fullfile('..', '..', 'Outputs');

% Save figures in .fig format
savefig(f1, fullfile(outputDir, 'innate_dynamics.fig'));
savefig(f2, fullfile(outputDir, 'Tcell_number.fig'));

% Save figures in .png format with a resolution of 300 dpi
print(f1, fullfile(outputDir, 'innate_dynamics.png'), '-dpng', '-r300'); 
print(f2, fullfile(outputDir, 'Tcell_number.png'), '-dpng', '-r300'); 

% Save optimized parameters to a text file
dlmwrite('optimized_parameters.txt', x_optim, 'delimiter', '\t');
clear f1 f2
end


%% Subfunctions -- Parameter Fitting

function scheme = defineSchemes()
% Define the dosing schemes, colors, and names.
scheme.T = {[12,12],[11,11],[12,12],[12,12],[7,7],[0,0],[12,0]}; 
scheme.numshot = {[7,7],[6,6],[4,4],[3,3],[2,2],[1,1],[7,1]}; 
scheme.k = {[1,1],[1,1],[1,1],[1,1],[log(4),log(4)],[0,0],[1,0]};

scheme.colors = {[0,0,0], [253, 128, 8]/256, [204, 102, 255]/256, [64, 0, 128]/256, ...
    [251, 2, 128]/256, [128, 64, 3]/256, [1,0,0]};

scheme.names = {'7-ED', '6-ED', '4-ED', '3-ED', '2-ED', 'Bolus', 'Adjuvant Bolus'};
end

%#######################################################

function [x_optim, chisq, params_guess, res, exitflags] = fitData(INNATE_PARAMS, scheme, Tfh_data)
%---------------------------------------
% Fit the data using Monte Carlo technique with fmincon.
%---------------------------------------
num_MC = 10;
data_index = [1,2,3,4,5,6];
[x_optim, chisq, params_guess, res, exitflags] = monte_carlo_fmincon(...
    @(x) objFunMLE(x, INNATE_PARAMS, scheme, Tfh_data, data_index), ...
    [], [], [], [], ...
    [1e5, 1e0]', [1e8, 1e3]', num_MC);

confidence = 0.95;
DOF = 7 - 2;
threshold = chi2inv(confidence, DOF);
if chisq < threshold
    disp('Model consistent with data');
else
    disp('Model not consistent with data');
end
end

%#######################################################

function [params, chisq, params_guess, res, exitflags] = monte_carlo_fmincon(...
        objfun, A, b, Aeq, beq, lb, ub, num)
%---------------------------------------
% Summary: Given the object function and linear constraints, repeat fmincon
%     optimization with random initialization for "num" times.
% Inputs: Object function to minimize, linear constraints, and number of
%     random initializations
% Output: Best optimized parameters, Chi-squared value, Array of all local
%     minimums found from MC approach, Array of residuals, Array of exit flags
%---------------------------------------

% Run optimizations
n_params = length(lb);
params_guess = zeros(n_params, num);
res = zeros(1,num);
exitflags = zeros(1,num);
for i=1:num
    x0 = exp(log(lb) + (log(ub)-log(lb)).*rand(n_params,1));
    [x, r, f] = fmincon(@(x)objfun(x), x0, A, b, Aeq, beq, lb, ub);
    params_guess(:,i) = x;
    res(i) = r;
    exitflags(i) = f;
end

% Obtain the best guess for the global minimum
[chisq,k] = min(res);
params = params_guess(:,k);
end

%#######################################################

function residual = objFunMLE(x, innate_params, scheme, Tfh_data, data_index)
%---------------------------------------
% Summary: Define the cost function used for the optimization of model 
%    parameters based on the data. We use MLE with gaussian error assumption. 
% Inputs: Optimized kinetic parameters, Fixed kinetic parameters, 
%    struct of dosing schemes, experimentally observed T cell number data, 
%    and indices of the schemes we are interested
% Outputs: residual
%---------------------------------------

% Get model predictions
numTcells = getModelPrediction(...
    x, innate_params, scheme, data_index);

% Get the chi-squared residual
data = Tfh_data(data_index,:);
data_mean = mean(data,2);
V = diag(var(data')); % We ignore off-diagonal correlations
V_inv = inv(V);
residual = (data_mean-numTcells)'*V_inv*(data_mean-numTcells);
end


% %######################################################
% 
% function numTcells = getModelPrediction(...
%     x, innate_params, scheme, data_index)
% %---------------------------------------
% % Inputs: Optimized kinetic parameters, Fixed kinetic parameters, 
% % struct of dosing schemes, and indices of the schemes we are interested
% % Outputs: Number of T cells, Total DCs 
% %---------------------------------------
% 
% innate_params = [innate_params, x']; % parameters
% % Get model predictions
% numTcells = zeros(length(data_index),1);   % for storing model prediction
% for j=1:length(data_index)
%     i = data_index(j);
%     [time,conc] = getInnateDynamics( ...
%         scheme.T{i}, scheme.numshot{i}, scheme.k{i}, innate_params);
%     numTcells(j) = conc(2100, end); % Take day 21 T cell number
% end
% end



% %% Subfunctions -- Plotting
% %#######################################################
% 
% function f = plotTcellNum(x, innate_params, scheme, Tfh_data, data_index, varargin)
% %---------------------------------------
% % Inputs: Optimized kinetic parameters, Fixed kinetic parameters, 
% %    struct of dosing schemes, experimentally observed T cell number data, 
% %    and indices of the schemes we are interested
% % Outputs: None
% % Summary: Plots a bar graph showing model-predicted and experimentally
% %     observed T cell numbers 
% %---------------------------------------
% 
% % Determines whether experimental data is plotted
% if length(varargin)==0
%     plot_exp_data = 1;
% else
%     plot_exp_data = varargin{1};
% end
% 
% % Plot the predicted T cell number along with experimentally measured data
% numTcells = getModelPrediction(...
%     x, innate_params, scheme, data_index);
% 
% f = figure;
% set(gcf, 'Units', 'centimeters','Position',[1,1,6,6])
% % Plot simulation data
% for i=1:length(data_index)
%     j = data_index(i);
%     b = bar(i, numTcells(i));
%     b.EdgeColor = scheme.colors{j};
%     b.FaceColor = [1,1,1];
%     hold on
% end
% 
% % Plot exprimental data
% if plot_exp_data
%     for i=1:length(data_index(data_index<=7))
%         j = data_index(i);
%         scatter(repmat(i,1,size(Tfh_data,2)), Tfh_data(j,:), ...
%         30, 'filled', 'square', 'MarkerFaceColor', scheme.colors{j})
%     end
% end
% 
% % Decorate
% ylabel('#Tfh cells','fontweight','bold')
% ylim([10^3, 10^6])
% xticks(1:length(data_index))
% xticklabels(scheme.names(data_index))
% a = get(gca, 'XTickLabel');
% set(gca, 'XTickLabel', a, 'FontWeight', 'bold')
% set(gca,'YScale','log')
% box off 
% end

% %#######################################################
% 
% function f_TcellDyn = plotInnateDyn(x, innate_params, scheme, scheme_idx)
% %---------------------------------------
% % Inputs: Optimized kinetic parameters, Fixed kinetic parameters, 
% %    struct of dosing schemes, experimentally observed T cell number data, 
% %    and indices of the schemes we are interested
% % Outputs: None
% % Summary: Plots a bar graph showing model-predicted and experimentally
% %     observed T cell numbers 
% %---------------------------------------
% set_log_scale = 1;
% innate_params = [innate_params, x'];
% 
% f_TcellDyn = figure;
% set(gcf, 'Units', 'centimeters','Position',[1,1,15.5,4])
% for i=flip(scheme_idx)
%     [time,conc] = getInnateDynamics( ...
%         scheme.T{i}, scheme.numshot{i}, scheme.k{i}, innate_params);
%     figure(f_TcellDyn)
%     subplot(1,4,4) %Tfh cells
%     plot(time, conc(:,8), 'LineWidth', 1, 'Color', scheme.colors{i})
%     if set_log_scale
%         set(gca, 'YScale', 'log')
%         ylim([10^2, 10^5])
%         yticks([10^2, 10^3, 10^4, 10^5])
%     end    
%     hold on
%     xlim([0,21])
%     xticks([0,7,14,21])
%     box off 
% 
%     subplot(1,4,3) %Proliferating T cells
%     plot(time, conc(:,7), 'LineWidth', 1, 'Color', scheme.colors{i})
%     if set_log_scale
%         set(gca, 'YScale', 'log')
%         yticks([10^2, 10^3, 10^4, 10^5])
%         ylim([10^2, 10^5])
%     end
%     hold on
%     xlim([0,21])
%     xticks([0,7,14,21])
%     box off 
% 
%     set_log_scale = 0;
% 
% 
%     subplot(1,4,2) %Ag+ DCs
%     plot(time, conc(:,6), 'LineWidth', 1, 'Color', scheme.colors{i})
%     hold on
%     xlim([0,14])
%     xticks([0,7,14])
%     box off 
% 
%     subplot(1,4,1) %Total DCs
%     plot(time, sum(conc(:,4:6),2), 'LineWidth', 1, 'Color', scheme.colors{i})
%     hold on
%     xlim([0,14])
%     xticks([0,7,14])
%     box off
% end
% subplot(1,4,4)
% axis_setting(gca, 'Number of Days', '# Tfh cells')
% subplot(1,4,3)
% axis_setting(gca, 'Number of Days', '# T cells')
% subplot(1,4,2)
% axis_setting(gca, 'Number of Days', '# aDC^{Ag+}s')
% leg = legend(scheme.names(flip(scheme_idx)), 'location', 'best');
% leg.ItemTokenSize = [10,5];
% subplot(1,4,1)
% axis_setting(gca, 'Number of Days', '# DCs')
% end
% 
% %#########################################################
% function axis_setting(ax, xname, yname)
%     xlabel(xname)
%     ylabel(yname)
%     hXAxis = get(ax, 'XAxis');
%     set(hXAxis.Label, 'FontSize', 9, 'FontWeight', 'bold')
%     hYAxis = get(ax, 'YAxis');
%     set(hYAxis.Label, 'FontSize', 9, 'FontWeight', 'bold')
% end