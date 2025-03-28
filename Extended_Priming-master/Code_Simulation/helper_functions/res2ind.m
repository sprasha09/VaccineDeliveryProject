function k = res2ind(gcBcellNum,residues,param)
%% Documentation
% Summary:
%   Helper function that allows easy indexing of gcBcellMutations array.
%   It has two uses: 
%    (1) If either all, first half, or second half of the residue
%        indices for a group of B cell is wanted
%    (2) If the index of specific residues of a B cell is wanted

% Inputs:
%  For use (1): 
%   gcBcellNum: row of indices of one or more B cells from same GC
%   residues: only 1:80, 1:40, 41:80 are allowed
%  For use (2):
%   gcBcellNum: single index of a B cell
%   residues: array of any numbers between 1 and 80
%  Common:
%   param: parameter struct

% Outputs: 
%   k: a row array containing desired indices of the residues

%First use: 
if isequal(residues, 1:param.n_res) || isequal(residues, 1:param.n_res/2) ||...
        isequal(residues, param.n_res/2+1:param.n_res) 
    k = repmat((gcBcellNum'-1)*param.n_res, 1, length(residues)) + residues;
    k = reshape(k',1,[]);
%Second use: 
else
    k = (gcBcellNum-1)*param.n_res+residues;
end
end

