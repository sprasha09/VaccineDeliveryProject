% Author: Leerang Yang
% Last modified 2023/08/12
% This code visualizes the impact of antigen release duration during
% extended release on the (Figs 6G-I)

%% Loading Data
% Load pSER simulation summary data and simply redundant data
load(fullfile('..','pSER_summary.mat'))
xlength = length(tspan{1,1});
tspan = tspan{1,1};
pSER_ag_release_times = [1, 3, 5, 10, 14];
n = length(pSER_ag_release_times);

% Summary statistics: Ag+ and Ag- B cells
Bcells = repmat({zeros(1, length(pSER_ag_release_times)+1)},1,2);
agconc = repmat({zeros(1, length(pSER_ag_release_times)+1)},1,2);

% We also want to plot the case when ag_release_time is zero. 
main_result = load(fullfile('..','summary.mat'));
%% Summarizing data and plotting
set(groot, 'DefaultAxesFontName', 'Arial')
set(groot, 'DefaultTextFontName', 'Arial')
f = []; % figures
colors = {[17, 128, 2]/256, [216,179,101]/256};

labels = cell(1,length(pSER_ag_release_times)+1);
labels{1} = '2-ED';
for j=1:length(pSER_ag_release_times)
    labels{j+1} = sprintf('Ag in %d days', pSER_ag_release_times(j));
end
%% ---------------- B cell number ----------------------
% First plot: Temporal trajectories of B cells for all cases
% Second plot: GC B cell fractions as a function of release time at a
%              specific time
% Third plot: GC B cell numbers as a function of release time at a
%              specific time

f(1) = figure; % First plot: Time Trajectories
set(gcf, 'Units', 'centimeters','Position',[1,1,14,7])
for ep=1:2
    d14_post_dose2_idx = xlength*(ep-1)+4*(7+14);
    subplot(1,2,ep)
    plot(tspan, main_result.totalnum{2}(xlength*(ep-1)+1:xlength*ep,1),...
        '-', 'LineWidth', 1);
    hold on
    Bcells{ep}(1,1) = main_result.totalnum{2}(d14_post_dose2_idx,1);
    for j=1:length(pSER_ag_release_times)
        plot(tspan, totalnum{1,j}(xlength*(ep-1)+1:xlength*ep,1),...
            '-', 'LineWidth', 1);
        xlabel('Day', 'fontsize', 9, 'fontweight', 'bold')
        % Get the final B cell numbers
        Bcells{ep}(1,j+1) = totalnum{1,j}(d14_post_dose2_idx,1);
    end
end
subplot(1,2,1)
legend(labels, 'location', 'best')
ylabel('Number of Ag+ GC B cells', 'fontsize', 9, 'fontweight', 'bold')
box off
subplot(1,2,2)
ylabel('Number of Ag- GC B cells', 'fontsize', 9, 'fontweight', 'bold')
savefig(f(1), 'Bcells')
box off

f(2) = figure; % Second plot: GC B cell fractions as a function of release time
set(gcf, 'Units', 'centimeters','Position',[6,6,3.5,5])
linestyles = {'-o', '-o'};
lines = gobjects(1,2);
for ep=1:2
    lines(1,ep) = plot([0,pSER_ag_release_times], 100*Bcells{ep}(1,:)./ ...
        (Bcells{1}(1,:)+Bcells{2}(1,:)), linestyles{ep}, 'Color', colors{ep},...
        'MarkerSize', 4, 'MarkerFaceColor', colors{ep}, 'LineWidth', 1.5);
    hold on
end
xlim([0,10])
ylim([0, 100])
xticks([0,5,10,15])
xlabel('Duration (Days)', 'fontsize', 9, 'fontweight', 'bold')
ylabel('% Ag^+ GC B cells', 'fontsize', 9, 'fontweight', 'bold')
leg = legend([lines(1,1), lines(1,2)],...
    {'Native', 'Non-Native'}, 'location', 'best');
leg.ItemTokenSize = [15,5];
leg.BoxFace.ColorType='truecoloralpha';
leg.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
box off

f(3) = figure; % Third plot: GC B cell numbers as a function of release time
set(gcf, 'Units', 'centimeters','Position',[1,1,3.5,5])
linestyles = {'-o', '-o'};
lines = gobjects(1,2);
for ep=1:2
    lines(1,ep) = plot([0,pSER_ag_release_times], Bcells{ep}(1,:), linestyles{ep}, 'Color', colors{ep},...
        'MarkerSize', 4, 'MarkerFaceColor', colors{ep}, 'LineWidth', 1.5);
    hold on
end
xlim([0,10])
ylim([0, 15*10^5])
xticks([0,5,10,15])
xlabel('Duration (Days)', 'fontsize', 9, 'fontweight', 'bold')
ylabel('#GC B cells', 'fontsize', 9, 'fontweight', 'bold')  
leg = legend([lines(1,1), lines(1,2)],...
    {'Native', 'Non-Native'}, 'location', 'best');
leg.ItemTokenSize = [15,5];
leg.BoxFace.ColorType='truecoloralpha';
leg.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
box off

savefig(f(1), 'Bcell_trajectories')
savefig(f(2), 'Bcell_fractions')
savefig(f(3), 'Bcell_numbers')

%% -------------------- Antigen concentration ----------------------
agconc{1}(1) = main_result.agconc_mean{2}(3,4*(7+14));
agconc{2}(1) = main_result.agconc_mean{2}(4,4*(7+14));
for j=1:length(pSER_ag_release_times)
    agconc{1}(j+1) = agconc_mean{j}(3,4*(7+14));
    agconc{2}(j+1) = agconc_mean{j}(4,4*(7+14));
end

f(4) = figure;
set(gcf, 'Units', 'centimeters','Position',[1,1,3.5,5])

linestyles = {'-o', '-o'};
lines = gobjects(2,2);
for ep=1:2
    lines(1,ep) = plot([0,pSER_ag_release_times], 100*agconc{ep}(1,:)./(agconc{1}(1,:)+agconc{2}(1,:)), ...
        linestyles{ep}, 'Color', colors{ep},...
        'MarkerSize', 4, 'MarkerFaceColor', colors{ep}, 'LineWidth', 1.5);
    hold on
end
xlim([0,10])

leg = legend([lines(1,1), lines(1,2)],...
    {'Native', 'Non-Native'}, 'location', 'best')
leg.ItemTokenSize = [15,5];
leg.BoxFace.ColorType='truecoloralpha';
leg.BoxFace.ColorData=uint8(255*[1 1 1 0.75]');
xlabel('Duration (Days)', 'fontsize', 9, 'fontweight', 'bold')
ylabel('% of Ag on FDC', 'fontsize', 9, 'fontweight', 'bold')
box off
savefig(f(3), 'Antigen')


