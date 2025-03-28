function result = combineResult(resultarr)
%   Combines multiple result structures into one.
% 
%   INPUTS:
%       resultarr - Cell array of result structures to be combined.
%       numDosing - Integer indicating the number of doses.
%
%   OUTPUTS:
%       result - Combined result structure.


% Extract common parameters
result.param = resultarr{1}.param;

% Combine concentration data
result.conc = cellfun(@(obj) obj.conc, resultarr, 'UniformOutput', false);

% Combine germinal center (GC) data
result.gc.affbytime = combineFields(resultarr, 'gc.affbytime');
result.gc.numbytime = combineFields(resultarr, 'gc.numbytime');
result.gc.numbylineage = combineFields(resultarr, 'gc.numbylineage');

% Combine naive B cell data
result.naive = combineFields(resultarr, 'naive');

% Combine output data
result.output.finalpc = combineFields(resultarr, 'output.finalpc');
result.output.finalmem = combineFields(resultarr, 'output.finalmem');
result.output.finalpb = combineFields(resultarr, 'output.finalpb');
result.output.pcnumbytime = combineFields(resultarr, 'output.pcnumbytime');
result.output.pcaffbytime = combineFields(resultarr, 'output.pcaffbytime');
result.output.memnumbytime = combineFields(resultarr, 'output.memnumbytime');
result.output.memaffbytime = combineFields(resultarr, 'output.memaffbytime');

% Combine dead plasma cell data
result.dead.plasmaCells = alignLeft(combineFields(resultarr, 'dead.plasmaCells', false), 1);

result.finalmemEGC = alignLeft(combineFields(resultarr, 'memoryCellsEGC', false), 1);
result.plasmaCellsEGC = alignLeft(combineFields(resultarr, 'plasmaCellsEGC', false), 1);
result.dead.plasmaCellsEGC = alignLeft(combineFields(resultarr, 'dead.plasmaCellsEGC', false), 1);
end


%% Subfunctions

function combined = combineFields(resultarr, field, isVertical)
%COMBINEFIELDS Combines fields of the result array.
% 
%   INPUTS:
%       resultarr - Cell array of result structures.
%       field - Name of the field to be combined.
%       isVertical (optional) - Boolean indicating if the combination should be vertical (default is true).
%
%   OUTPUTS:
%       combined - Combined data.

% Check for default isVertical value
if nargin < 3
    isVertical = true;
end

data = cellfun(@(obj) getNestedField(obj, field), resultarr, 'UniformOutput', false);
if isVertical
    combined = padvertcat(data);
else
    combined = horzcat(data{:});
end

end


function value = getNestedField(structure, field)
%   Gets the value of a nested field from a structure.
% 
%   INPUTS:
%       structure - Structure from which the field value should be extracted.
%       field - Name of the nested field.
%
%   OUTPUTS:
%       value - Value of the nested field.

fields = strsplit(field, '.');
for i = 1:length(fields)
    structure = structure.(fields{i});
end
value = structure;

end


function A = padvertcat(A)
%   Summary: Vertically concatenates numerical arrays in a cell array,
%   ensuring they have the same size along the second dimension.
%
%   A = PADVERTCAT(A) takes a cell array A, where each cell contains a
%   numerical array (matrix) with potentially different sizes along the
%   second dimension. It vertically concatenates these arrays while
%   ensuring they have the same size along the second dimension. If the
%   arrays have different sizes, the function pads them with zeros to match
%   the size of the largest array along the second dimension before
%   concatenation.
%
%   Inputs:
%       A - Cell array containing numerical arrays (matrices).
%
%   Outputs:
%       A - Vertically concatenated numerical array with uniform size along
%           the second dimension.

lengths = cellfun(@(x) size(x, 2), A);
maxlength = max(lengths);

% If all arrays have the same length, combine them vertically
if all(lengths == maxlength)
    A = vertcat(A{:});
else
    % If not, pad them to the same size with zeros and combine
    A = cellfun(@(x) [x, zeros(size(x, 1), maxlength - size(x, 2), size(x, 3))], A, 'UniformOutput', false);
    A = vertcat(A{:});
end

end