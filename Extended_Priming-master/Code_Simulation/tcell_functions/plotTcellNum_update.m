%% Subfunctions -- Plotting
%#######################################################

function f = plotTcellNum_update(x, innate_params, scheme, Tfh_data, data_index, varargin)
%---------------------------------------
% Inputs: Optimized kinetic parameters, Fixed kinetic parameters,
%    struct of dosing schemes, experimentally observed T cell number data,
%    and indices of the schemes we are interested
% Outputs: None
% Summary: Plots a bar graph showing model-predicted and experimentally
%     observed T cell numbers
%---------------------------------------

% Determines whether experimental data is plotted
if length(varargin)==0
    plot_exp_data = 1;
else
    plot_exp_data = varargin{1};
end

x_dim = size(innate_params,1);
f = figure;

for plotIndex = 1:x_dim

    innate_params_temp =  [innate_params(plotIndex,:), x'];
    % Plot the predicted T cell number along with experimentally measured data
    numTcells = getModelPrediction(...
        x, innate_params_temp, scheme, data_index);


    set(gcf, 'Units', 'centimeters','Position',[1,1,6,6])
    % Plot simulation data
    for i=1:length(data_index)
        j = data_index(i);
        b = bar(i, numTcells(i));
        b.EdgeColor = scheme.colors{j};
        b.FaceColor = [1,1,1];
        hold on
    end

    % Plot exprimental data
    if plot_exp_data
        for i=1:length(data_index(data_index<=7))
            j = data_index(i);
            scatter(repmat(i,1,size(Tfh_data,2)), Tfh_data(j,:), ...
                30, 'filled', 'square', 'MarkerFaceColor', scheme.colors{j})
        end
    end
end

% Decorate
ylabel('#Tfh cells','fontweight','bold')
ylim([10^3, 10^6])
xticks(1:length(data_index))
xticklabels(scheme.names(data_index))
a = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', a, 'FontWeight', 'bold')
set(gca,'YScale','log')
box off
end
