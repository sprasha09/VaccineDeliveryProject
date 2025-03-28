% ------------------------------------------------------
% Script to summarize simulations specifically regarding varying pSER
% release conditions. 
% This script analyzes the scanning results of the pSER parameter and generates relevant summary statistics. 
% ------------------------------------------------------

% Add necessary paths for dependencies
initializeDependencies();

% Set fixed parameters
p = baseCaseParameters();
p.T = 7;
p.k = 1.3863;
p.numshot = 2;
p.pSER_adj_release_time = 0;
p.pSER = 1;

% Load results and obtain summary statistics
summaryFile = fullfile('..', 'pSER_summary.mat');
if ~exist(summaryFile, 'file')
    % Parameters that are varied
    pSER_ag_release_times = [1, 2, 3, 5, 7, 10, 14];

    % Initialize result storage
    n = length(pSER_ag_release_times);
    tspan = cell(1,n);
    totalnum = cell(1,n);
    agconc = cell(1,n);
    abtiter = cell(1,n);
    agconc_mean = cell(1,n);
    abtiter_mean = cell(1,n);


    % Summarize results for each pSER_ag_release_time
    for j = 1:n
        [tspan{j}, totalnum{j}, agconc{j}, abtiter{j}, agconc_mean{j}, abtiter_mean{j}] = ...
            summarizeResultsForParameter(j, pSER_ag_release_times, p);
    end
    
    % Save Results
    save(summaryFile, 'tspan', 'agconc', 'abtiter', 'abtiter_mean', 'agconc_mean', 'totalnum', '-v7.3');
end

%% Subfunctions

function initializeDependencies()
    % Add dependencies to path
    addpath(fullfile('..','..','Code_Parameter_Generation'));
    addpath(fullfile('..','..','Code_Simulation','sim_initialization_functions'));
    addpath(fullfile('..','..','Code_Simulation','helper_functions'));
end


function [tspan, totalnum, agconc, abtiter, agconc_mean, abtiter_mean] =...
    summarizeResultsForParameter(j, pSER_ag_release_times, p)
    % Summarizes the results for a given pSER_ag_release_time.
    
    p.pSER_ag_release_time = pSER_ag_release_times(j);

    % Define file path
    dirnm = sprintf('pSER_Ag_%d_Adj_%d/', p.pSER_ag_release_time, p.pSER_adj_release_time);
    fnm = getFileLocation(p);
    fnm = fnm{1};
    fnm = insertAfter(fnm, 'Data/', dirnm);
    fnm = fullfile('..','..','Code_Simulation',fnm);
    dirnm = fileparts(fnm);

    % Load and combine results
    result = compileGCRunResults(p.first, p.last, dirnm);
    result = combineResult(result);

    % Summarize result
    [tspan, totalnum, agconc, abtiter, agconc_mean, abtiter_mean, ~, ~, ~] = summarize(result);
end