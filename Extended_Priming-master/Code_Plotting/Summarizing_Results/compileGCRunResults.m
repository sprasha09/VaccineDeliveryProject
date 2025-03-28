function compiledResults = compileGCRunResults(firstRun, lastRun, directory)
% Combines smaller result files from a folder to create a compiled dataset.
%
% INPUTS:
%   firstRun  - Starting number of the GC run range
%   lastRun   - Ending number of the GC run range
%   directory - Name of the directory containing the .mat files
%
% OUTPUT:
%   compiledResults (1xN cell array) - Each cell contains the result of a single GC run

% Initialize variables
expectedLength = lastRun - firstRun + 1;
existingRuns = zeros(1, expectedLength);
files = dir(fullfile(directory, '*.mat'));
fileNames = {files.name};
compiledResults = cell(1, length(fileNames));

for i = 1:length(fileNames)
   % Regular expression to parse the file name
   filePattern = '(?<startRun>\d*)_to_(?<endRun>\d*).mat';
   token = regexp(fileNames{i}, filePattern, 'names');
   startOfFile = str2double(token.startRun);
   endOfFile = str2double(token.endRun);
   
   if startOfFile >= firstRun && endOfFile <= lastRun
       existingRuns(startOfFile:endOfFile) = 1;
       loadedData = load(fullfile(directory, fileNames{i}));
       compiledResults{i} = loadedData.result;
       compiledResults{i}.startRun = startOfFile;
   end
end

% Check if any runs are missing in the range
missingRuns = find(existingRuns == 0);
if ~isempty(missingRuns)
   fprintf('Missing files for runs: %s in directory: %s\n', num2str(missingRuns), directory)
   error('Incomplete dataset due to missing files.') 
end

% Remove empty cells
nonEmptyIndices = cellfun(@(x) ~isempty(x), compiledResults);
compiledResults = compiledResults(nonEmptyIndices);
end
