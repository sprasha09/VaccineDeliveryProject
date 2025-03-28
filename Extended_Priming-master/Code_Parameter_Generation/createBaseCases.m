% ---------------------------------------------------
% Creates two parameter files "Base_cases.txt" and
% "Base_cases_summarize.txt" in the folder "../parameters/"
% They are used to run simulations and summarize the results, respectively.
% ---------------------------------------------------

p = baseCaseParameters();
configs = {struct('numfrag', 10, 'name', 'Base_cases');
    struct('numfrag', 1, 'name', 'Base_cases_summarize')};

scenarios = {
    struct('T', 0, 'k', 0, 'numshot', 1); % Bolus
    struct('T', 7, 'k', log(4), 'numshot', 2); % 2-ED d0_d7 20_80
    struct('T', 12, 'k', log(4), 'numshot', 2); % 2-ED d0_d12 20_80
    struct('T', 12, 'k', 1, 'numshot', 7); % 7-ED
    struct('T', 7, 'k', log(4), 'numshot', 2, 'pSER', 1); % 2-ED + only the last dose slow release
    struct('T', 7, 'k', log(4), 'numshot', 2, 'pSER', 2); % 2-ED + both doses slow release
    };

for i=1:2
    config = configs{i};
    txt_name = getFileName(config.name);
    p.numfrag = config.numfrag;
    for j=1:length(scenarios)
        p = updateParameters(p, scenarios{j});
        createParameters(txt_name, 'a', p);
    end
end


%% Subfunctions


function p = updateParameters(p, scenario)
% Summary: Overwrites or adds fields in a parameter struct based on a scenario struct.
%
% Inputs:
%   p        : Struct containing initial parameters.
%   scenario : Struct containing values that need to be overwritten or added to 'p'.
%
% Outputs:
%   p        : Updated parameter struct with values from the scenario struct.

    fields = fieldnames(scenario);
    for i=1:numel(fields)
        p.(fields{i}) = scenario.(fields{i});
    end
end



function [txt_name] = getFileName(baseName)
% Summary: Generates the full path for a parameter file in the 'parameters' directory.
%
% Inputs:
%   baseName : Base name of the file (without the '.txt' extension).
%
% Outputs:
%   txt_name : Full path of the file, assuming it's in a 'parameters' directory in the parent folder.
%
% Description:
%   The function constructs the file name by appending '.txt' to the provided base name, then 
%   determines the full path in the 'parameters' directory. If the file already exists, it will be deleted.

    fileName = sprintf('%s.txt', baseName);
    txt_name = fullfile('..', 'parameters', fileName);
    if exist(txt_name, 'file')
        delete(txt_name);
    end
end



function createParameters(txt_name, writeoption, p)
% Summary: Writes parameter values into a text file for simulations.
%
% Inputs:
%   txt_name    : Name of the parameter text file. It can be with or without the '.txt' extension.
%   writeoption : File write mode. Use 'a' for appending to existing content or 'w' to overwrite.
%   p           : Struct containing simulation parameters. The struct should have the following fields:
%                 - first: Starting range
%                 - last: Ending range
%                 - numfrag: Number of fragments or divisions between 'first' and 'last'
%                 - vaxnum, T, k, numshot, E1h, dE12, p, masking, C0, w1, w2, steric,
%                   memToGCFrac, outputprob, outputpcfrac, rho, pSER, tmax
%
% Description:
%   The function will generate a list of parameters based on the struct 'p'. These parameters will be 
%   written into the specified text file (txt_name) for use in simulations. The range specified by 
%   'first' and 'last' is fragmented into smaller parts defined by 'numfrag', and parameters are written 
%   for each fragment.
%
%   If the 'parameters' directory does not exist in the parent folder, it will be created.
%
% Note:
%   The function does not return any value. It writes the parameters directly into the text file.
%


% Ensure the parameters directory exists
if ~isfolder(fullfile('..', 'parameters'))
    mkdir(fullfile('..', 'parameters'));
end

% Ensure the filename ends with '.txt'
txt_name = erase(txt_name, '.txt');
fid = fopen([txt_name, '.txt'], writeoption);

% Determine the range
first_arr = linspace(p.first, p.last + 1, p.numfrag + 1);
first_arr(end) = [];
numrepeat = (p.last - p.first + 1) / p.numfrag;

% Write the parameters
for first = first_arr
    last = first + numrepeat - 1;
    inputs = num2cell([p.vaxnum, p.T, p.k, p.numshot, p.E1h, p.dE12, p.p, ...
          p.masking, p.C0, p.w1, p.w2, p.steric, p.memToGCFrac,p.outputprob, ...
          p.outputpcfrac, p.rho, p.pSER, p.tmax, first, last]);
    inputs = cellfun(@(x) num2str(x), inputs, 'UniformOutput', false);
    fprintf(fid, [strjoin(inputs(1:end),'\t'),'\n']);
end

fclose(fid);

end