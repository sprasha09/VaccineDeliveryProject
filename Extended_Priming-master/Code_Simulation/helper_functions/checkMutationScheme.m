function flag = checkMutationScheme(gcBcells, k, i, naiveBcells, ...
    mutations, gcBcellsMutation, M_GC, param)
%% Documentation
% Summary:
%  For a specific GC B cell, check if the affinities are correct, by
%  comparing it against affinities calculated from the mutation state and
%  corresponding fitness landscape.

% Outputs:
%  flag: 1 if the affinities are not calculated correctly.

% Inputs:
%  gcBcells, naiveBcells: Arrays of GC and naive B cells
%  k, i: Index of GC, Index of B cell
%  mutations: 1x2 cell array containing mutation sizes
%  gcBcellsMutation: 2D array containing the mutation states of GC B cells
%  M_GC: Number of GCs 
%  param: parameter struct

flag = 0;
BcellMutations = gcBcellsMutation(k,convertMutNum(i,1:param.n_res, param));
lineage = gcBcells(k,i);
E_from_mut(1) = naiveBcells(k+M_GC*2,lineage) - sum(mutations{1}(k,lineage,BcellMutations==1));
E_from_mut(2) = param.f0 - sum(mutations{2}(k,lineage,BcellMutations==1));
if round(E_from_mut(1)*1000) ~= round(gcBcells(k+2*M_GC, i)*1000)
    flag = 1;
    fprintf(fid, sprintf('WT affinity not correct; k=%d, i=%d\n', k, i));
    fprintf(fid, ['mutations: ', num2str(find(BcellMutations)), '\n']);
end
if round(E_from_mut(2)*1000) ~= round(gcBcells(k+3*M_GC, i)*1000)
    flag = 1;
    fprintf(fid, sprintf('Variant affinity not correct; k=%d, i=%d\n', k, i));
    fprintf(fid, ['mutations: ', num2str(find(BcellMutations)), '\n']);
end
if flag==1
    error('affinity is wrong')
end
end