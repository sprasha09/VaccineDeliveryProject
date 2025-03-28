function result = updateHistory(result, gcBcells, plasmaCells, memoryCells,...
        plasmaCellsEGC, memoryCellsEGC, agconc, abconc, Ka, Ka_var,...
        agconc_Epmask, tspan_summary, storeidx, param)
%% Documentation
% Summary:
%   Store the summary statistics of the current status of the simulation.
%   This function is called repeatedly at a given time interval
% Output:
%   result: Struct array that summarizes the outcome of the simulation. 
%   ***** See the documentation for "runGCs.m" for details. *****
%       
% Inputs:
%   result, and various arrays and values from the simulation

%% Initialization
n = length(tspan_summary);
M_GC = param.M_GC;
if isempty(result)
    result.param = param;
    result.gc.numbytime = zeros(M_GC, param.n_ep*n, 4);
    result.gc.affbytime = zeros(M_GC, param.n_ep*n, 4);
    result.gc.numbylineage = zeros(M_GC, n, 2010);
    result.conc.concarray = zeros(4, param.n_ep+2, n);
    result.conc.concarray_Epmask = zeros(2, param.n_ep, n);
    result.conc.Kaarray = zeros(3, param.n_ep, n);
    result.conc.Kaarray_var = zeros(3, param.n_ep, n);
    result.output.pcnumbytime = zeros(2, param.n_ep*n, 4); %dim1: GC-target1, GC-target2, EGC-target1, EGC-target2
                                                           %dim3: affinities                                                           %Affinities
    result.output.memnumbytime = zeros(2, param.n_ep*n, 4);
    result.output.pcaffbytime = zeros(2, param.n_ep*n, 4);
    result.output.memaffbytime = zeros(2, param.n_ep*n, 4);
end

%% Recording results
% Concentration
result.conc.concarray(:,:,storeidx) = [agconc; [zeros(3,2),abconc]];
result.conc.Kaarray(:,:,storeidx) = Ka;
result.conc.Kaarray_var(:,:,storeidx) = Ka_var;
result.conc.concarray_Epmask(:,:,storeidx) = agconc_Epmask;

% Number and affinity of GC B cells
if any(any(gcBcells))
    [numbyaff, affprct] = cellsNumAff(gcBcells, M_GC, param);
    result.gc.numbytime(:, storeidx+n*(0:(param.n_ep-1)), :) = numbyaff;
    result.gc.affbytime(:, storeidx+n*(0:(param.n_ep-1)), :) = affprct;
    
    lineage = gcBcells(1:M_GC, :);
    for k=1:M_GC
        result.gc.numbylineage(k,storeidx,:) = ...
            histcounts(lineage(k,:), 1:2011);
    end
end
buildup3D = @(A,d3) permute(reshape(A',size(A',1),[],d3),[2,1,3]);
flatten1D = @(A) reshape(permute(A,[3,2,1]), size(A,3),[]);

% Number and affinity of memory and plasma cells
PCs = {flatten1D(buildup3D(plasmaCells, param.pcfieldnum)); %GC-derived
       plasmaCellsEGC}; %EGC-derived
MEMs = {flatten1D(buildup3D(memoryCells, param.memfieldnum)); %GC-derived
        memoryCellsEGC}; %EGC-derived
for i=1:2 %GC-derived, EGC-derived
    [numbyaff, affprct] = cellsNumAff(PCs{i}, 1, param);
    result.output.pcnumbytime(i, storeidx+n*(0:(param.n_ep-1)), :) =...
        numbyaff;
    result.output.pcaffbytime(i, storeidx+n*(0:(param.n_ep-1)), :) =...
        affprct;
    [numbyaff, affprct] = cellsNumAff(MEMs{i}, 1, param);
    result.output.memnumbytime(i, storeidx+n*(0:(param.n_ep-1)), :) =...
        numbyaff;
    result.output.memaffbytime(i, storeidx+n*(0:(param.n_ep-1)), :) =...
        affprct;
end

end

function [numbyaff, affprct] = cellsNumAff(cellsarr, M, param)
% Obtain the summary of number and affinities of B cells
% Outpus:
%   numbyaff: 1x2x4 array; Dim1,2,3 - GC, Epitope, 
%             # of B cells with affinities greater than 6, 7, 8, 9
%   affprct: 1x2x4 aray; Dim1,2,3 - GC, Epitope,
%             100, 90, 75, 50 percentiles affinities
% Inputs:
%   cellsarr: 2D array of B cells, each column representing a B cell.
%             Can be GC, memory, or plasma cells
%   M: Number of GC/EGC
%   param: parameter struct

    thresholds = [6,7,8,9];
    percentile = [100, 90, 75, 50];
    aff = cellsarr(M*2+1:M*3,:);
    target = cellsarr(M*1+1:M*2,:);
    
    numbyaff = zeros(M, param.n_ep, 4);
    affprct = zeros(M, param.n_ep, 4);
    
    for i=1:4
        for ep = 1:param.n_ep
            numbyaff(:,ep,i) = sum((aff.*(target==ep))>thresholds(i),2);
            for k=1:M
                affprct(k, ep, i) = prctile((aff(k, target(k,:)==ep))',percentile(i))';
            end
        end
    end
end