%#######################################################

function f_TcellDyn = plotInnateDyn_update(x, innate_params, scheme, scheme_idx)
%---------------------------------------
% Inputs: Optimized kinetic parameters, Fixed kinetic parameters,
%    struct of dosing schemes, experimentally observed T cell number data,
%    and indices of the schemes we are interested
% Outputs: None
% Summary: Plots a bar graph showing model-predicted and experimentally
%     observed T cell numbers
%---------------------------------------
set_log_scale = 1;
%innate_params = [innate_params, x'];

x_dim = size(innate_params,1);
f_TcellDyn = figure;

for plotIndex = 1:x_dim


    set(gcf, 'Units', 'centimeters','Position',[1,1,15.5,4])

    for i=flip(scheme_idx)
        innate_params_temp =  [innate_params(plotIndex,:), x'];
        [time,conc] = getInnateDynamics_update( ...
            scheme.T{i}, scheme.numshot{i}, scheme.k{i}, innate_params_temp);

        if plotIndex == 1
            color_ED = [0, 0, 0];% Black for Base

        elseif plotIndex == 2
            color_ED = [1, 0.4, 0.6];  % Pink for 7-ED low
        else
            color_ED = [0, 1, 0];  % green for 7-ED high
        end



        figure(f_TcellDyn)
        subplot(1,6,6) %Tfh cells
        plot(time, conc(:,8), 'LineWidth', 1, 'Color', color_ED)
        % if set_log_scale
        %     set(gca, 'YScale', 'log')
        %     %ylim([10^2, 10^5])
        %     %yticks([10^2, 10^3, 10^4, 10^5])
        % end
        hold on
        xlim([0,21])
        xticks([0,7,14,21])
        box off

        subplot(1,6,5) %Proliferating T cells
        plot(time, conc(:,7), 'LineWidth', 1, 'Color', color_ED)
        % if set_log_scale
        %     set(gca, 'YScale', 'log')
        %     yticks([10^2, 10^3, 10^4, 10^5])
        %     ylim([10^2, 10^5])
        % end
        hold on
        xlim([0,21])
        xticks([0,7,14,21])
        box off

        set_log_scale = 0;


        subplot(1,6,4) %Ag+ DCs
        plot(time, conc(:,6), 'LineWidth', 1, 'Color', color_ED)
        hold on
        xlim([0,14])
        xticks([0,7,14])
        box off

        subplot(1,6,3) %Total DCs
        plot(time, sum(conc(:,4:6),2), 'LineWidth', 1, 'Color', color_ED)
        hold on
        xlim([0,14])
        xticks([0,7,14])
        box off

        subplot(1,6,3) %Total Ags
        plot(time, conc(:,1), 'LineWidth', 1, 'Color', color_ED)
        hold on
        xlim([0,14])
        xticks([0,7,14])
        box off

        subplot(1,6,3) %Total Ags
        plot(time, conc(:,2), 'LineWidth', 1, 'Color', color_ED)
        hold on
        xlim([0,14])
        xticks([0,7,14])
        box off



    end
end
subplot(1,6,6)
axis_setting(gca, 'Number of Days', '# Tfh cells')
legend_entries = [scheme.names(flip(scheme_idx)), scheme.names_low(flip(scheme_idx)),scheme.names_high(flip(scheme_idx))];
legend(legend_entries, 'Location', 'Best');
subplot(1,6,5)
axis_setting(gca, 'Number of Days', '# T cells')
legend_entries = [scheme.names(flip(scheme_idx)), scheme.names_low(flip(scheme_idx)),scheme.names_high(flip(scheme_idx))];
legend(legend_entries, 'Location', 'Best');
subplot(1,6,4)
axis_setting(gca, 'Number of Days', '# aDC^{Ag+}s')
legend_entries = [scheme.names(flip(scheme_idx)), scheme.names_low(flip(scheme_idx)),scheme.names_high(flip(scheme_idx))];
legend(legend_entries, 'Location', 'Best');
subplot(1,6,3)
axis_setting(gca, 'Number of Days', '# DCs')
legend_entries = [scheme.names(flip(scheme_idx)), scheme.names_low(flip(scheme_idx)),scheme.names_high(flip(scheme_idx))];
legend(legend_entries, 'Location', 'Best');
subplot(1,6,2)
axis_setting(gca, 'Number of Days', '# Ags')
%legend_entries = [scheme.names(flip(scheme_idx)), scheme.names_low(flip(scheme_idx)),scheme.names_high(flip(scheme_idx))];
%legend(legend_entries, 'Location', 'Best');
subplot(1,6,1)
axis_setting(gca, 'Number of Days', '# Adjs')
%legend_entries = [scheme.names(flip(scheme_idx)), scheme.names_low(flip(scheme_idx)),scheme.names_high(flip(scheme_idx))];
%legend(legend_entries, 'Location', 'Best');
end

%#########################################################
function axis_setting(ax, xname, yname)
xlabel(xname)
ylabel(yname)
hXAxis = get(ax, 'XAxis');
set(hXAxis.Label, 'FontSize', 9, 'FontWeight', 'bold')
hYAxis = get(ax, 'YAxis');
set(hYAxis.Label, 'FontSize', 9, 'FontWeight', 'bold')
end