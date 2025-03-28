function [result,dirname] = loadResult(p)
% Summary: Load simulation result data based on given parameters.
%
% INPUT:
% p : structure containing parameters for determining the data location.
%
% OUTPUT:
% result : combined simulation results
% dirname : directory name where the data was loaded from
%
% DESCRIPTION:
% The function determines the file path for the simulation results based on 
% the provided parameters. It then loads the result data and combines it 
% based on the number of vaccinations.

% Initialize parameters 
param = initializeParameters(p.vaxnum, p.T, p.k, p.numshot, p.E1h, p.dE12, p.p, ...
 p.masking, p.C0, p.w1, p.w2, p.steric, p.memToGCFrac,...
 p.outputprob, p.outputpcfrac, p.rho, p.pSER, p.tmax, p.first, p.last);

% Construct Data Directory Path and Validate the existence of the directory
fnm = getFileLocation(param);
dirname = fileparts(fullfile('..', '..', 'Code_Simulation', fnm{1}));
fprintf('%s\n', dirname)
if not(isfolder(dirname))
    error('Data Path Does Not Exist');
end

% Load Result
result = compileGCRunResults(param.first, param.last, dirname);
result = combineResult(result);
end