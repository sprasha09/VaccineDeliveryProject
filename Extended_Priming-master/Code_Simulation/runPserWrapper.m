function runPserWrapper(T, k, numshot, pSER_ag_release_time, pSER_adj_release_time, first, last)
%--------------
% Wrapper function used to study pSER cases with variations in parameters.
%--------------

% Define parameters
addpath(fullfile('..','Code_Parameter_Generation'));
addpath('sim_initialization_functions')
addpath('helper_functions')
p = baseCaseParameters();
p.T = T; 
p.k = k;
p.numshot = numshot;
% ######Check this!!########
p.pSER = 1;
% ##########################
param = initializeParameters(p.vaxnum, p.T, p.k, p.numshot, p.E1h, p.dE12, p.p, ...
          p.masking, p.C0, p.w1, p.w2, p.steric, p.memToGCFrac,p.outputprob, ...
          p.outputpcfrac, p.rho, p.pSER, p.tmax, first, last);
param.pSER_ag_release_time = pSER_ag_release_time;
param.pSER_adj_release_time = pSER_adj_release_time;


% Run simulation
result = runGCsMain(param);

% Define the file path
fnm = getFileLocation(param);
fnm = fnm{param.vaxnum};
dirnm = sprintf('pSER_Ag_%d_Adj_%d/', pSER_ag_release_time, pSER_adj_release_time);
fnm = insertAfter(fnm, 'Data/', dirnm);
mkdir(fileparts(fnm));

% Save result
save(fnm, 'result')
end