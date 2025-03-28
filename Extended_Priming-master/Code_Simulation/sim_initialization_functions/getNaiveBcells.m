function [naivebcells, mutations] = getNaiveBcells(param)
%% Documentation
% Summary
%  Initialize the naive B cells with their lineages, targets, affinities,
%  and mutational fitness landscapes

% Outputs
%  naiveBcells: 200 x 2210 x 7 array
%   Dim 1 - GCs, Dim 2 - B cells, Dim 3 - Properties
%  mutations - 1x2 cell array of 200 x 2210 x 80 array 
%              Dim1 - GCs; Dim2 - B cells; Dim3 - Residues
%              Effects of mutations against WT and variant
%       
% Inputs
%  param: parameter struct 

%% Parameters and Constants
M_GC = param.M_GC;
E1h = param.E1h;
if param.E1h > 8 
    E1h = 7.2;
end
dE12 = param.dE12;
p = param.p;
f0 = param.f0;
NaiveMax = param.NaiveMax; 
MemoryMax = param.MemoryReentryMax;

n_res = param.n_res; % number of residues 
dE = 0.2;            % class size
classnum = 11;       % total number of class bins
N = 2000;            % total number of naive precursors

%% Find out the number of B cells in each class
maxclasses = round(linspace(E1h-f0, E1h-f0-dE12, param.n_ep)/dE)+1;
% Bin at which the value of frequency is 1 for dominant/subdominant cells
fitnessarray = linspace(f0, f0+2, classnum); %6 to 8 with interval of 0.2
options = optimoptions('fsolve','Display','off');
r = zeros(param.n_ep,1); % Slopes of the geometric distribution
r(1) = fsolve(@(x) (N)-(x^maxclasses(1)-1)/(x-1), 1.1, options);
for i=2:param.n_ep
    if maxclasses(i)>1
        r(i) = fsolve(@(x) (N)-(x^maxclasses(i)-1)/(x-1), 1.1, options);
    else
        r(i) = N;
    end
end %NOTE: the code is written to be able to generalize for n_ep>2
naivebcellsarr = zeros(param.n_ep,classnum); %2 x 11 array, number of 
                                    % naive B cells in each fitness class
naivebcellsarr(1,1:classnum) = ...
    (1-p*(param.n_ep-1))*r(1).^(maxclasses(1)-(1:(classnum)));
for i=2:param.n_ep
    if maxclasses(i) > 1
        naivebcellsarr(i,1:classnum) = p*r(i).^(maxclasses(i)-(1:(classnum)));
    elseif maxclasses(i)==1
        naivebcellsarr(i,1) = p*N;
    end
end

naivebcells = zeros(M_GC, NaiveMax+MemoryMax, param.naivefieldnum);
            % Array of naive B cells 
            %Dimension 1: GC number
            %Dimension 2: Lineage number
            %Dimension 3: properties % (1)lineage, (2)target, 
            %(3)WT affinity, (4)Variant affinity, (5)activated time;

MutationPDF = param.MutationPDF;
Sigma{1} = [1, 0.4;
        0.4, 1];
Sigma{2} = [1, param.rho;
        param.rho, 1];
mutations = cell(1,2);
mutations{1} = zeros(M_GC, NaiveMax+MemoryMax, n_res); %WT mutations
mutations{2} = zeros(M_GC, NaiveMax+MemoryMax, n_res); %Variant mutations

for i=1:M_GC %For each GC
    naivebcellsint = floor(naivebcellsarr+rand(size(naivebcellsarr)));
    if param.E1h > 8 
        naivebcellsint(1,end) = 1;
    end
     %Stochastically round up or down the frequency to get integer numbers
    idx = 1;
    for type=1:param.n_ep %For dominant and subdominant epitopes
        for j=1:length(fitnessarray)
            idx_new = idx+naivebcellsint(type,j);
            naivebcells(i, idx:idx_new-1, 1) = idx:idx_new-1; %Lineage
            naivebcells(i, idx:idx_new-1, 2) = type; %Target
            naivebcells(i, idx:idx_new-1, 3) = fitnessarray(j); %WT aff
            naivebcells(i, idx:idx_new-1, 4) = f0; %Variant aff

            X = MutationPDF(1) + mvnrnd([0,0],MutationPDF(2)^2*...
                Sigma{1+(type==param.n_ep)}, (idx_new-idx)*n_res);
            dE = log10(exp(1))*(exp(X) - MutationPDF(3));
            for strain=1:2
                mutations{strain}(i, idx:idx_new-1, :) = ...
                    reshape(dE(:,strain), idx_new - idx, n_res);
            end
            idx = idx_new;
        end
    end
end
end