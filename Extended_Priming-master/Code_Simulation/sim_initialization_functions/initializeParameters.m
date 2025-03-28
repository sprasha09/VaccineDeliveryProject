function param = initializeParameters(varargin)
%% Documentation
% Summary:
%   Initializes the variables based on inputs.
%   Also initialize parameters whose values are fixed. 

%% Initializes Parameters for Simulation
%
% This function initializes variables and parameters based on the provided inputs.
%
% ## Usage
% param = initializeParameters(vaxnum, T, k, numshot, E1h, dE12, p, ...
%                              masking, C0, w1, w2, steric, memToGCFrac,...
%                              outputprob, outputpcfrac, rho, pSER, tmax, first, last);

% Inputs
%
% vaxnum - Integer (1-4): Indicates dose number. ----- ONLY 1 FOR THIS STUDY
% T, k, numshot - Scalars: Dosing scheme parameters
% E1h - Scalar (6-8): Defines naive B cell germline affinities.
% dE12, p - Scalars (0-1): Naive B cell germline affinities.
% masking - Binary: 1 for epitope masking, 0 otherwise. ----- ONLY 0 FOR THIS STUDY
% C0 - Scalar: Reference antigen concentration.
% w1 - Scalar: Defines saturation of antigen capture (typically 0).
% w2 - Scalar (0.3-1): Selection stringency (typically 0.5).
% steric - Scalar (0-1): 0 if no masking; varies if masking is 1. ----- ONLY 0 FOR THIS STUDY
% memToGCFrac - Scalar: Fraction of pre-existing memory cells entering GC. ----- ONLY 0 FOR THIS STUDY
% outputprob, outputpcfrac - Scalars: B cell exit and conversion to plasma cell.
% rho - Scalar (typically 0.95): Epitope conservation level for subdominant. ----- NOT RELEVANT FOR THIS STUDY
% pSER - Scalar: Number of doses given as slow-release (from the final dose).
% tmax - Scalar: Time duration of simulation.
% first, last - Integers: Indexes of GCs (e.g., 201-400).

% Output:
% param - Struct: Contains simulation parameters and constants. Refer to the function body for detailed fields.

%% Initialization
[vaxnum, T, k, numshot, E1h, dE12, p, ...
 masking, C0, w1, w2, steric, memToGCFrac,...
 outputprob, outputpcfrac, rho, pSER, tmax, first, last] ...
 = deal(varargin{:});
param = struct();

%% Simulation numbers 
param.first = first; 
param.last = last;
param.M_GC = last-first+1; % Number of GCs;

%% Data storage
param.N_GC_MAX = 3000; % Max B cell number per GC
param.N_PC_MAX = 10000; % Max PC number per GC 
param.N_PB_MAX = 1000; % Max PB number per GC 
% (However PB production is not considered in this simulation) 
param.naivefieldnum = 7; % Number of properties that naive B cells have 
param.gcfieldnum = 5; % Number of properties that GC B cells have
param.pcfieldnum = 6; % Number of properties that PCs have
param.memfieldnum = 9; % Number of properties that memory cells have

%% Dosing scheme
param.vaxnum = vaxnum;
param.T = T;
param.k = k;
param.numshot = numshot;
param.pSER = pSER;

%% Naive B cell affinity Distribution
param.E1h = E1h;
param.dE12 = dE12;
param.p = p;
param.n_ep = 2; % Number of epitopes. Dominant and subdominant epitope
param.NaiveMax = 2010; % Max number of Naive B cells per GC 

%% Re-entry of memory B cells into GCs
param.memToGCFrac = memToGCFrac;
param.MemoryReentryMax = 200; % Max number of pre-existing memory cells 
                              % that can re-enter GC, per GC

%% Antibody Production and Effect
param.production = 10; % Parameter that is used to change the rate 
param.delay = 2;
param.masking = masking;
param.steric = steric;
param.outputprob = outputprob;
param.outputpcfrac = outputpcfrac;
param.outputpbfrac = outputpcfrac;
param.mutprob = 1; % Probability of mutation of a daughter cell. Fixed 
param.IgM0 = 0.01; % nM; Initial IgM amount 
param.r_IgG = param.IgM0*param.production*1; %nM per PC per day;
param.r_IgG_EGC = param.r_IgG*1;
param.r_IgM = param.IgM0*param.production*1; %nM per PB per day;


%% Antigen Dose and Effect
param.Ag0 = 10; % nM; Initial soluble Ag amount
param.Ageff = 0.01; % Effectiveness of soluble Ag at B cell activation
param.C0 = C0;
param.F0 = nan; % Not used in this simulation; for slow delivery
param.dose = []; % Not used in this simulation; for slow delivery
param.dose_t = []; % Not used in this simulation; for slow delivery


%% Time
param.tmax = tmax;
param.dt = 0.01; % day; time interval of simulation
param.tspan_summary = 0:(25*param.dt):param.tmax; % time points at which
                               % statistics from the simulation are saved
param.pSER_ag_release_time = 10;
param.pSER_adj_release_time = 0;

%% Activation and Competition Parameters
param.f0 = 6; % Reference binding affinity 
param.activation_threshold = 1; % Amount of antigen captured for activation
param.w1 = w1;
param.w2 = w2;

%% Mutation Parameters
param.MutationPDF = [3.1, 1.2, 3.08];
param.rho = rho;
param.n_res = 80; % Number of residues that each B cell has

%% Reaction Rates
param.k_deposit = 24; % per day; Rate at which IC is transported to FDC
param.d_Ag = 3; % per day; Rate of antigen decay 
param.d_IC = 0.15; % per day; Rate of decay for IC on FDC
param.d_IgM = log(2)/4; % per day; Decay rate of IgM; Half-life 4 days 
param.d_IgG = log(2)/4; %per day; Decay rate of IgG; Half-life 28 days
param.d_pc = log(2)/4; %per day; Decay rate of PCs; Half-life 4 days

%% Recruitment, birth, death rates
param.lambdamax = 1; % per day; Max rate of GC entry for naive B cells
param.betamax = 4; % per day; Max rate of GC B cell positive selection.
                      % Division of up to ~4 times a day is possible
param.betamaxEGC = 4; % per day
param.mu = 0.5; % per day; Death rate of GC B cells 
% param.Nmax = 10; % per day; Defines number of naive B cells that enter GC
param.numTmax = 1200; % Max non-dimensionalized helper T cell availability
param.naiveprolnum = 4; % Number of copies made when Naive B cell enters GC
end