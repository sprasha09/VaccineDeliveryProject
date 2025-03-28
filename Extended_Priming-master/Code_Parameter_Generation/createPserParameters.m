% ---------------------------------------------------
% Create a text file "pSER_params.txt" containing simulation parameters for 
% simulating extended dosing cases with control over antigen or adjuvant
% release time. 
% ---------------------------------------------------

% Parameters that are fixed
T = 7; k = log(4); numshot = 2; pSER_adj_release_time = 0;

% Define text file path
filePath = fullfile('..', 'parameters', 'pSER_params.txt');
if exist(filePath, "file")
    delete(filePath);
end
fid = fopen(filePath, 'a');

% Define the repeat ranges
first = 1; last = 2000; numfrag = 10;
first_arr = linspace(first, last+1, numfrag+1);
first_arr(end) = [];
numrepeat = (last-first+1)/numfrag;

% Write the text file
for pSER_ag_release_time = [1, 2, 3, 5, 7, 10, 14]
    for first = first_arr
        last = first + numrepeat - 1;
        parameters = num2cell([T, k, numshot, pSER_ag_release_time, pSER_adj_release_time, first, last]);
        parameters = cellfun(@(x) num2str(x), parameters, 'UniformOutput', false);
        fprintf(fid, [strjoin(parameters(1:end),'\t'),'\n']);
    end
end

% Close the text file
fclose(fid);