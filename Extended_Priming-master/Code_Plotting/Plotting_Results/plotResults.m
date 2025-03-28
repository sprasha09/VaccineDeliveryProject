% Author: Leerang Yang
% Last modified 2023/09/14
% Description: This code is used to generate various data visualization plots
% summarizing simulation results. It imports data from the "summary.mat" file
% and creates plots to visualize the simulation results.

% SUMMARY OF THE FIGURES PLOTTED:
% 1. Antigen Amount over time: Soluble and IC concentrations for bolus,
% 2-ED, and 7-ED. (Fig. 4B)
% 2. Antigen Amount Bar Plot: Maximum antigen amount (native, non-native)
% on FDC for bolus, 2-ED, ED. (Fig. 4D)
% 3. Antibody Concentration Over Time: Antibody titer for bolus, 2-ED, and
% ED. (Fig. 4C)
% 4. GC B cell: Displays total GC B cells, Antigen+ GC B cells, and
% Fraction of Antigen+ GC B cells. (Fig. 4E and F)
% 5. Ag+ B cell fraction bar graph: Shows the fraction of Ag+ GC B cells
% for different conditions. (Fig. 4G)
% 6. Antibody and Antigen Amount Over Time For Extended 2nd Dose: Displays
% titer for Native and Non-Native for the extended second dose. (Fig. 6B-C)
% 7. Ag amount on FDC bar graphs for extended 2nd dose: Shows the antigen
% amount on FDC for the regular and extended second doses. (Fig. 6D)
% 8. Total GC B cells and Ag+ B cell fraction bar graph for Extended 2nd
% Dose (Fig. 6E-F)
% 9. Number of GC B cells when Both doses are extended for 2-dose ED (Fig.
% S5D-E)



% -------------------------------------------
%% Figure Settings
% Set the default font to Arial for all figures
set(groot, 'DefaultAxesFontName', 'Arial');
set(groot, 'DefaultTextFontName', 'Arial');

% Define custom colors for plot lines
colors = {[128, 64, 3]/256;    % Bolus
          [251, 2, 128]/256;   % 0-7
          [102, 204, 255]/256; % 0-12
          [0, 0, 0];           % 7-ED
          [17, 128, 2]/256;    % pSER 2nd dose
          [1, 0, 0]};          % pSER both doses


%% Import Data
% Import summary results from GC/EGC simulations
% The following variables are imported:
%   - abtiter_mean: Mean antibody titer
%   - agconc_mean: Mean antigen
%   - gcnum: Number of B cells in individual GCs
%   - totalnum: Total number of GC B cells
%   - tspan: Time span for simulation
%
% Each variable is a 1x6 cell array, with indices corresponding to:
%   1: Bolus
%   2: 0-7
%   3: 0-12
%   4: 2 weeks
%   5: pSER 2nd dose only
%   6: pSER both doses
load(fullfile('..', 'summary.mat'));

%% Antigen Amount over time (Fig. 4B)
% Plot shows soluble and IC concentrations for bolus, 2-ED, and 7-ED over
% time in three different subplots.
figure
tiledlayout(1,4);
set(gcf, 'Units', 'centimeters','Position',[1,1,20,5])
titles = {'Bolus', '2-ED', '', '7-ED'};
markers = {'-.',':','-','--'};
for i=[1,2,4]
    nexttile;
    for j=1:4
    plot(tspan{i}, agconc_mean{i}(j,:), markers{j}, 'LineWidth', 1, color=colors{i})
    hold on
    end
    set(gca,'Yscale', 'log')
    xlabel('Number of Days', 'fontweight', 'bold')
    ylabel('Antigen Amount', 'fontweight', 'bold')
    ylim([10^-6, 1])
    xlim([0,14])
    title(titles{i})
    yticks(10.^[-6, -4, -2, 0])
    xticks([0,7,14])
    box off
    grid on
end
leg = legend({'Soluble Native', 'Soluble Non-Native', 'IC Native', 'IC Non-Native'});
leg.ItemTokenSize = [15,5];
leg.Position(1:2) = [0.7, 0.3];


%% Antigen Amount Bar Plot (Fig. 4D)
%  Bar plots of the maximum antigen amount (native, non-native) on FDC
%  for bolus, 2-ED, ED
f = figure;
set(f, 'Units', 'centimeters','Position',[1,1,5,5])
idcs = [1,2,4];
h = {};
agMeanMax = zeros(length(idcs), 2);
for i=1:length(idcs)
    agMeanMax(i,:) = [max(agconc_mean{idcs(i)}(3,:)), max(agconc_mean{idcs(i)}(4,:))];
    h{i} = bar([i], agMeanMax(i,:));
    hold on
end
[CI, pvals] = agStatisticalAnalysis(agconc, agconc_mean, idcs);
numGroups = 3;
numBars = 2;
% Calculate the width of all the bars
groupWidth = min(0.8, numBars/(numBars + 1.5));
figure(f)
hold on;
for j = 1:numBars
    % Calculate center of each bar
    x = (1:numGroups) - groupWidth/2 + (2*j-1) * groupWidth / (2*numBars);
    errorbar(x, agMeanMax(:,j), CI(:,j), 'k', 'linestyle', 'none');
end

for i=1:3
    set(h{i}(1), 'FaceColor', colors{idcs(i)})
    set(h{i}(2), 'EdgeColor', colors{idcs(i)})
    set(h{i}(2), 'FaceColor', [1,1,1])
    set(h{i}(2), 'LineWidth', 0.5)
end
hLegend = legend([h{3}(1), h{3}(2)],...
    {'Native', 'Non-Native'}, 'orientation', 'horizontal', 'NumColumns', 1,...
    'location', 'northoutside');
hLegend.ItemTokenSize = [10, 5];
xticks([1,2,3])
xticklabels({'Bolus', '2-ED', '7-ED'})
xtl = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', xtl, 'FontSize', 9, 'FontWeight', 'bold')
set(gca, 'Yscale', 'log')
ylim([10^-6, 1])
ylabel({'Antigen Amount on FDC'})
hXAxis = get(gca, 'XAxis');
set(hXAxis.Label, 'FontSize', 9, 'FontWeight', 'bold')
hYAxis = get(gca, 'YAxis');
set(hYAxis.Label, 'FontSize', 9, 'FontWeight', 'bold')
box off



%% Antibody Concentration Over Time (Fig. 4C)
% Antibody titer for bolus, 2-ED, and ED over time, shown in three separate
% panels
figure
tiledlayout(1,4);
set(gcf, 'Units', 'centimeters','Position',[1,1,20,5])
titles = {'Bolus', '2-ED', '', '7-ED'};
for i=[1,2,4]
    nexttile;
    ln(1) = semilogy(tspan{i}, abtiter_mean{i}(1,:), '-','LineWidth', 1.5, color=colors{i});
    hold on
    ln(2) = semilogy(tspan{i}, abtiter_mean{i}(2,:), '--', 'LineWidth', 1.5, color=colors{i});
    xlim([0,14])
    ylim([10^-3, 10^3])
    xticks([0,7,14])
    xlabel('Number of Days', 'fontweight', 'bold')
    ylabel('Titer (Arbitrary Unit)', 'fontweight', 'bold')
    title(titles{i})
    box off
    grid on
    yticks([10^-3, 10^-1, 10^1, 10^3])
end
leg = legend({'Native Antigen Binding', 'Non-Native Antigen Binding'});
leg.Position(1:2) = [0.7, 0.3];

%% GC B cell (Fig. 4E, F)
% (1) Total GC B cells
% (2) Antigen+ GC B cells

g = gobjects(1,2);
for i=1:2
    g(i) = figure;
    set(gcf, 'Units', 'centimeters','Position',[1,1,4,5])
    box off
end
styles = {'-', '--'};
for i=[1,2,4]
    num = {};
    for ep=1:2
        figure(g(ep));
        xlength = length(tspan{i});
        num{ep} = totalnum{i}(xlength*(ep-1)+1:xlength*ep,1);
    end
    figure(g(1))
    semilogy(tspan{i}, num{1}, 'LineWidth', 1.5, 'Color', colors{i});
    hold on

    figure(g(2))
    semilogy(tspan{i}, num{1}+num{2}, 'LineWidth', 1.5, 'Color', colors{i});
    hold on
end
figure(g(1))
xlabel('Number of Days', 'fontweight', 'bold', 'fontsize', 9)
ylabel('#Ag^+ GC B cells', 'fontweight', 'bold', 'fontsize', 9)
xlim([0,21])
ylim([1e4, 1e7])
xticks([0,7,14,21])
% legend({'Bolus', '2-ED', 'ED'}, 'location', 'best')
box off
grid on

figure(g(2))
xlabel('Number of Days', 'fontweight', 'bold', 'fontsize', 9)
ylabel('# GC B cells', 'fontweight', 'bold', 'fontsize', 9)
xlim([0,21])
ylim([1e4, 1e7])
xticks([0,7,14,21])
% legend({'Bolus', '2-ED', 'ED'}, 'location', 'best')
box off
grid on



%% Ag+ B cell fraction bar graph (Fig. 4G)
nums = repmat({[0,0]},1,4);
xlength = length(tspan{1});
for ep=1:2
    for j=1:4
        nums{j}(ep) = totalnum{j}(xlength*(ep-1)+4*21+1);
    end
end

[CI,pvals] = gcStatisticalAnalysis(gcnum, xlength, idcs);

f = figure
set(gcf, 'Units', 'centimeters','Position',[1,1,4.3,5])
x = [1,2,3];
y = cell2mat(nums([1,2,4])')./sum(cell2mat(nums([1,2,4])'),2);
colors1 = colors([1,2,4]);


for i=1:3
b = bar(x(i),y(i,1)*100, 0.4, "basevalue", 0.05);
b.FaceColor = colors1{i};
hold on
end

errorbar(x, y(:,1)*100, CI(:,3), 'k', 'linestyle', 'none');

xlim([0.5,3.5])
ylim([1, 100])
set(gca,'YScale', 'log')
xticks(x);
xticklabels({'Bolus', '2-ED', '7-ED'})
ylabel({'% Ag^+ GC B Cells'}, 'fontsize', 9)
xtl = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', xtl, 'FontSize', 9, 'FontWeight', 'bold')
box off


%% Antibody and Antigen Amount Over Time For Extended 2nd Dose (Fig. 6B, 6C)
figure
set(gcf, 'Units', 'centimeters','Position',[1,1,4.5,5])
colors2 = {colors{5}, [216,179,101]/256}; % Native and Non-native colors
for ep=1:2
    semilogy(tspan{5}, abtiter_mean{5}(ep,:), '-','LineWidth', 1.5, color=colors2{ep})
    hold on
end
ylabel('Titer (Arbitrary Unit)')
xlabel('Number of Days')
leg = legend({'Native', 'Non-Native'}, 'location', 'best');
leg.ItemTokenSize = [15,5];
box off
hXAxis = get(gca, 'XAxis');
set(hXAxis.Label, 'FontSize', 9, 'FontWeight', 'bold')
hYAxis = get(gca, 'YAxis');
set(hYAxis.Label, 'FontSize', 9, 'FontWeight', 'bold')
xlim([0,21])
xticks([0,7,14,21])


figure
set(gcf, 'Units', 'centimeters','Position',[1,1,7.2,5])
markers = {'-.',':','-','--'};
for i=1:4
    semilogy(tspan{5}, agconc_mean{5}(i,:), markers{i}, 'LineWidth', 1.5, color=colors2{2-mod(i,2)})
    hold on
end
xlim([0,21])
xticks([0,7,14,21])
ylim([10^-6, 1])
xlabel('Number of Days')
ylabel('Antigen Amount')
box off
hLegend=legend({sprintf('Soluble\nNative'), sprintf('Soluble\nNon-Native'),...
    sprintf('IC\nNative'), sprintf('IC\nNon-Native')}, 'location', 'eastoutside');
set(hLegend.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.7]));
hLegend.ItemTokenSize = [15,5];
hXAxis = get(gca, 'XAxis');
set(hXAxis.Label, 'FontSize', 9, 'FontWeight', 'bold')
hYAxis = get(gca, 'YAxis');
set(hYAxis.Label, 'FontSize', 9, 'FontWeight', 'bold')
%% Ag amount on FDC bar graphs for extended 2nd dose (Fig. 6D)
f = figure
set(gcf, 'Units', 'centimeters','Position',[1,1,5,5])
h = {};
idcs = [2,5];
agMeanMax = zeros(length(idcs), 2);
for i=1:length(idcs)
    agMeanMax(i,:) = [max(agconc_mean{idcs(i)}(3,:)), max(agconc_mean{idcs(i)}(4,:))];
    h{i} = bar([i], agMeanMax(i,:));
    hold on
end
[CI, pvals] = agStatisticalAnalysis(agconc, agconc_mean, idcs);
numGroups = size(agMeanMax,1);
numBars = size(agMeanMax,2);
% Calculate the width of all the bars
groupWidth = min(0.8, numBars/(numBars + 1.5));
figure(f)
hold on;
for j = 1:numBars
    % Calculate center of each bar
    x = (1:numGroups) - groupWidth/2 + (2*j-1) * groupWidth / (2*numBars);
    errorbar(x, agMeanMax(:,j), CI(:,j), 'k', 'linestyle', 'none');
end

% h{1} = bar([1], [max(agconc_mean{2}(3,:)), max(agconc_mean{2}(4,:))]);
% hold on
% h{2} = bar([2], [max(agconc_mean{5}(3,:)), max(agconc_mean{5}(4,:))]);
% idcs = [2,5];
for i=1:2
    set(h{i}(1), 'FaceColor', colors{idcs(i)})
    set(h{i}(2), 'EdgeColor', colors{idcs(i)})
    set(h{i}(2), 'FaceColor', [1,1,1])
    set(h{i}(2), 'LineWidth', 0.5)
end
hLegend = legend([h{2}(1), h{2}(2)],...
    {'Native', 'Non-Native'}, 'orientation', 'horizontal', 'NumColumns', 1,...
    'location', 'northwest');
hLegend.ItemTokenSize = [10, 5];
xticks([1,2,3])
xticklabels({'2-ED', 'Dose 2 Extended'})
xtl = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', xtl, 'FontSize', 9, 'FontWeight', 'bold')
set(gca, 'Yscale', 'linear')
ylim([0, 1])
ylabel({'Antigen Amount on FDC'})
hXAxis = get(gca, 'XAxis');
set(hXAxis.Label, 'FontSize', 9, 'FontWeight', 'bold')
hYAxis = get(gca, 'YAxis');
set(hYAxis.Label, 'FontSize', 9, 'FontWeight', 'bold')
box off

disp(sprintf('ratio increases from %.2f to %.2f', ...
    max(agconc_mean{2}(3,:))/(max(agconc_mean{2}(3,:))+max(agconc_mean{2}(4,:))), ...
    max(agconc_mean{5}(3,:))/(max(agconc_mean{5}(3,:))+max(agconc_mean{5}(4,:))) ))


%% Total GC B cells and Ag+ B cell fraction bar graph for Extended 2nd Dose (Fig. 6E, 6F)

idcs = [1,2,5];
nums = zeros(length(idcs),2);
for ep=1:2
    for j=1:length(idcs)
        nums(j,ep) = totalnum{idcs(j)}(xlength*(ep-1)+4*21+1);
    end
end
[CI,pvals] = gcStatisticalAnalysis(gcnum, xlength, idcs);

numGroups = 3;
numBars = 2;
% Calculate the width of all the bars
groupWidth = min(0.8, numBars/(numBars + 1.5));


f = figure
set(gcf, 'Units', 'centimeters','Position',[1,1,5,5])
h = {};
for i=1:3
    h{i} = bar(i, nums(i,:));
    hold on
end
figure(f)
hold on;
for j = 1:numBars
    % Calculate center of each bar
    x = (1:numGroups) - groupWidth/2 + (2*j-1) * groupWidth / (2*numBars);
    errorbar(x, nums(:,j), CI(:,j), 'k', 'linestyle', 'none');
end

for i=1:3
    set(h{i}(1), 'FaceColor', colors{idcs(i)})
    set(h{i}(2), 'EdgeColor', colors{idcs(i)})
    set(h{i}(2), 'FaceColor', [1,1,1])
    set(h{i}(2), 'LineWidth', 0.5)
end
hLegend = legend([h{3}(1), h{3}(2)],...
    {'Native', 'Non-Native'}, 'orientation', 'horizontal', 'NumColumns', 1,...
    'location', 'northwest');
hLegend.ItemTokenSize = [10, 5];
xlim([0,4])
ylim([0,2*10^6])
xticks(x);
xticklabels({'Bolus', '2-ED', 'Dose 2 Extended'})
ylabel({'# GC B cells'})
hXAxis = get(gca, 'XAxis');
set(hXAxis.Label, 'FontSize', 9, 'FontWeight', 'bold')
hYAxis = get(gca, 'YAxis');
set(hYAxis.Label, 'FontSize', 9, 'FontWeight', 'bold')
xtl = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', xtl, 'FontSize', 9, 'FontWeight', 'bold')
set(gca,'YScale','Linear')
box off



f = figure
set(gcf, 'Units', 'centimeters','Position',[1,1,5,5])
xlim([0,4])
x = [1,2,3];
y = nums(:,1)./sum(nums,2);
colors1 = colors([1,2,5]);
for i=1:3
    b = bar(x(i),y(i)*100, 0.5, 'BaseValue', 10^-1);
    b.FaceColor = colors1{i};
    hold on
end
errorbar(x, y*100, CI(:,3), 'k', 'linestyle', 'none');
xticks(x);
xticklabels({'Bolus', '2-ED', 'Dose 2 Extended'})
ylabel({'% Ag^+ GC B cells'})
hXAxis = get(gca, 'XAxis');
set(hXAxis.Label, 'FontSize', 9, 'FontWeight', 'bold')
hYAxis = get(gca, 'YAxis');
set(hYAxis.Label, 'FontSize', 9, 'FontWeight', 'bold')
xtl = get(gca, 'XTickLabel');
set(gca, 'XTickLabel', xtl, 'FontSize', 9, 'FontWeight', 'bold')
set(gca,'YScale','Linear')
box off

%% Number of GC B cells when Both doses are extended for 2-dose ED
% (1) Total GC B cells
% (2) Antigen+ GC B cells

g = gobjects(1,2);
for i=1:2
    g(i) = figure;
    set(gcf, 'Units', 'centimeters','Position',[1,1,6,5])
end

xlength = length(tspan{2});
for i=[6,5]
    num = {};
    for ep=1:2
        figure(g(ep));
        num{ep} = totalnum{i}(xlength*(ep-1)+1:xlength*ep,1);
    end
    figure(g(1))
    semilogy(tspan{i}, num{1}, 'LineWidth', 1.5, 'Color', colors{i});
    hold on

    figure(g(2))
    semilogy(tspan{i}, num{1}+num{2}, 'LineWidth', 1.5, 'Color', colors{i});
    hold on
end

figure(g(1))
xlabel('Number of Days', 'fontweight', 'bold')
ylabel('# Ag+ GC B cells', 'fontweight', 'bold')
ylim([0,5*10^6])
xlim([0,21])
xticks([0,7,14,21])
leg = legend({'Alum-pSer (both doses)', 'Alum-pSer (dose 2)'}, 'location', 'best');
leg.ItemTokenSize = [15,5];

figure(g(2))
xlabel('Number of Days', 'fontweight', 'bold')
ylabel('# GC B cells', 'fontweight', 'bold')
ylim([0,5*10^6])
xlim([0,21])
xticks([0,7,14,21])
leg = legend({'Alum-pSer (both doses)', 'Alum-pSer (dose 2)'}, 'location', 'best');
leg.ItemTokenSize = [15,5];






%% Subfunction

function [CI,pvals] = gcStatisticalAnalysis(gcnum, xlength, idcs)
    n_group = length(idcs);
    n_repeat = size(gcnum{1},1);
    native = zeros(n_repeat,n_group);
    nonNative = zeros(n_repeat,n_group);
    for i=1:n_group
        native(:,i) = squeeze(gcnum{idcs(i)}(:,xlength*(1-1)+4*21+1,1));
        nonNative(:,i) = squeeze(gcnum{idcs(i)}(:,xlength*(2-1)+4*21+1,1));
    end
    blockSum = @(data) squeeze(sum(reshape(data, [200, 10, n_group]),1));
    native = blockSum(native);
    nonNative = blockSum(nonNative);
    getSEM = @(data_matrix) std(data_matrix) ./ sqrt(size(data_matrix, 1));
    
    fraction = native./(native+nonNative)*100;
    pvals.native = anova_posthoc(native);
    pvals.nonNative = anova_posthoc(nonNative);
    pvals.fraction = anova_posthoc(fraction);

    CI = [getSEM(native)', getSEM(nonNative)', getSEM(fraction)']*tinv(1-0.05/2, 10-1);
end

function [CI, pvals] = agStatisticalAnalysis(agconc, agconc_mean, idcs)
% Description: Analyzes the provided data to compare native and non-native groups using one-way ANOVA and Tukey's post hoc test.
% Input: 
%   - agconc: Cell array containing the concentration data for all groups.
%   - agconc_mean: Cell array containing the mean concentration data.
%   - idcs: Indices of the groups to be considered for statistical analysis.
% Output:
%   - pvals: Structure containing the results of Tukey's post hoc test for both native and non-native groups.
    n_group = length(idcs);
    n_repeat = length(agconc{1});
    native = zeros(n_repeat,n_group);
    nonNative = zeros(n_repeat,n_group);
    for i=1:n_group
        [~,j1] = max(agconc_mean{idcs(i)}(3,:));
        [~,j2] = max(agconc_mean{idcs(i)}(4,:));
        native(:,i) = cellfun(@(x) x(3,j1), agconc{idcs(i)});
        nonNative(:,i) = cellfun(@(x) x(4,j2), agconc{idcs(i)});
    end

    getSEM = @(data_matrix) std(data_matrix) ./ sqrt(size(data_matrix, 1));
    CI = [getSEM(native)', getSEM(nonNative)']*tinv(1-0.05/2, 10-1);

    if n_group>2
        pvals.native = anova_posthoc(native);
        pvals.nonNative = anova_posthoc(nonNative);
    elseif n_group==2
        [~, pvals.native, ~, ~] = ttest2(native(:,1), native(:,2));
        [~, pvals.nonNative, ~, ~] = ttest2(nonNative(:,1), nonNative(:,2));
    end
end

function c = anova_posthoc(data)
% Function: anova_posthoc
% Description: Performs one-way ANOVA on the given data and then runs Tukey's 
% post hoc analysis.
% Input: 
%   - data: 2D array where each column contains data for one group.
% Output:
%   - c: Results from Tukey's post hoc test.
    hiddenFig = figure('Visible', 'off');
    [p,tbl,stats] = anova1(data, [], 'off');
    [c,m,h,nms] = multcompare(stats, 'Display', 'off');
    close(hiddenFig);
end