%#######################################################

function f_TcellDyn = plotInnateDyn(x, innate_params, scheme, scheme_idx)
%---------------------------------------
% Inputs: Optimized kinetic parameters, Fixed kinetic parameters, 
%    struct of dosing schemes, experimentally observed T cell number data, 
%    and indices of the schemes we are interested
% Outputs: None
% Summary: Plots a bar graph showing model-predicted and experimentally
%     observed T cell numbers 
%---------------------------------------
set_log_scale = 1;
innate_params = [innate_params, x'];

f_TcellDyn = figure;
set(gcf, 'Units', 'centimeters','Position',[1,1,15.5,4])
for i=flip(scheme_idx)
    [time,conc] = getInnateDynamics( ...
        scheme.T{i}, scheme.numshot{i}, scheme.k{i}, innate_params);
    figure(f_TcellDyn)
    subplot(1,4,4) %Tfh cells
    plot(time, conc(:,8), 'LineWidth', 1, 'Color', scheme.colors{i})
    if set_log_scale
        set(gca, 'YScale', 'log')
        ylim([10^2, 10^5])
        yticks([10^2, 10^3, 10^4, 10^5])
    end    
    hold on
    xlim([0,21])
    xticks([0,7,14,21])
    box off 

    subplot(1,4,3) %Proliferating T cells
    plot(time, conc(:,7), 'LineWidth', 1, 'Color', scheme.colors{i})
    if set_log_scale
        set(gca, 'YScale', 'log')
        yticks([10^2, 10^3, 10^4, 10^5])
        ylim([10^2, 10^5])
    end
    hold on
    xlim([0,21])
    xticks([0,7,14,21])
    box off 

    set_log_scale = 0;


    subplot(1,4,2) %Ag+ DCs
    plot(time, conc(:,6), 'LineWidth', 1, 'Color', scheme.colors{i})
    hold on
    xlim([0,14])
    xticks([0,7,14])
    box off 

    subplot(1,4,1) %Total DCs
    plot(time, sum(conc(:,4:6),2), 'LineWidth', 1, 'Color', scheme.colors{i})
    hold on
    xlim([0,14])
    xticks([0,7,14])
    box off
end
subplot(1,4,4)
axis_setting(gca, 'Number of Days', '# Tfh cells')
subplot(1,4,3)
axis_setting(gca, 'Number of Days', '# T cells')
subplot(1,4,2)
axis_setting(gca, 'Number of Days', '# aDC^{Ag+}s')
leg = legend(scheme.names(flip(scheme_idx)), 'location', 'best');
leg.ItemTokenSize = [10,5];
subplot(1,4,1)
axis_setting(gca, 'Number of Days', '# DCs')
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