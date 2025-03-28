% ------------------------------------------------------
% SUMMARY:
% This script analyzes and summarizes the main simulation results across 
% different immunization scenarios, including bolus, 2 doses (0-7 days apart),
% 2 doses (0-12 days apart), and every day. The results include metrics like
% antigen and antibody concentrations, and numbers of immune cells over time. 
% Once the results for all scenarios have been summarized, they are saved in a
% summary.mat file.
% ------------------------------------------------------

% Initialization
initializeDependencies();
p = baseCaseParameters();
scenarios = getScenarios();

% Analyze each scenario
for i = 1:numel(scenarios)
    p = updateParameters(p, scenarios{i});
    [result, ~] = loadResult(p);
    [tspan{i}, totalnum{i}, agconc{i}, abtiter{i}, agconc_mean{i}, abtiter_mean{i}, IgM_mean{i}, IgG_mean{i}, gcnum{i}] = summarize(result);
end

% Save results
save(fullfile('..','summary.mat'), 'tspan', 'totalnum',...
    'IgM_mean', 'IgG_mean', 'agconc', 'abtiter', 'agconc_mean', 'abtiter_mean', 'gcnum', '-v7.3');


%% Subfunctions

function initializeDependencies()
    % Add dependencies to path
    addpath(fullfile('..','..','Code_Parameter_Generation'));
    addpath(fullfile('..','..','Code_Simulation','sim_initialization_functions'));
    addpath(fullfile('..','..','Code_Simulation','helper_functions'));
end


function scenarios = getScenarios()
    scenarios = {
        struct('T', 0, 'k', 0, 'numshot', 1); % Bolus
        struct('T', 7, 'k', log(4), 'numshot', 2); % 2-ED d0_d7
        struct('T', 12, 'k', log(4), 'numshot', 2); % 2-ED d0_d12
        struct('T', 12, 'k', 1, 'numshot', 7); % 7-ED
        struct('T', 7, 'k', log(4), 'numshot', 2, 'pSER', 1); % 2-ED d0_d7 alum-pSER 2nd dose only
        struct('T', 7, 'k', log(4), 'numshot', 2, 'pSER', 2); % 2-ED d0_d7 alum-pSER both doses
    };
end


function p = updateParameters(p, scenario)
    fields = fieldnames(scenario);
    for j = 1:numel(fields)
        p.(fields{j}) = scenario.(fields{j});
    end
end
