function fnm = getFileLocation(param)
%% Documentation
% Summary: Determines the location/name of data files based on input parameters.
%
% Inputs:
%  - param: A structure containing various parameters
%
% Outputs:
%  - fnm: 1 x vaxnum cell array. Contains the file location/name of data
%        files for all previous vaccinations and the current vaccination.
%
% Notes:
% The function constructs the filename based on various conditions and 
% parameters, and checks if the file already exists.

% Default vaccination schedule (in days)
vaxtiming = [0, 28, 180, 180]; 

% Initialize directories based on vaccination number
data_locs = {'Data_Prime', 'Data_Secondary', 'Data_Tertiary', 'Data_Vax4'};

for i = 1:param.vaxnum
    % Check for memory cell re-entry
    if param.memToGCFrac > 0 && i > 1
        data_locs{i} = fullfile('memToGC', sprintf('memToGCFrac_%.3f', param.memToGCFrac), data_locs{i});
    end
    
    % Check for pSER condition
    if param.pSER
        data_locs{i} = fullfile(sprintf('pSER%d', param.pSER), data_locs{i});
    end
    
    % Check for epitope overlap
    if param.steric ~= 0
        data_locs{i} = fullfile('steric', sprintf('steric_%.2f', param.steric), data_locs{i});
    end
    
    % Check for alternative model of Ag capture
    if param.w1 > 0
        data_locs{i} = fullfile('agCaptureSaturation', data_locs{i});
    end
    
    data_locs{i} = fullfile('Data', data_locs{i});
end

% Convert input parameters to string for use in filename
inputs = num2cell([param.vaxnum, param.T, param.k, param.numshot, param.E1h, param.dE12, param.p, ...
    param.masking, param.C0, param.w1, param.w2, param.steric, param.outputprob, ...
    param.outputpcfrac, param.rho, param.tmax]);
inputs = cellfun(@num2str, inputs, 'UniformOutput', false);

dirnm = cell(1, param.vaxnum);
fnm = cell(1, param.vaxnum);

dirnm{param.vaxnum} = fullfile(data_locs{param.vaxnum}, strjoin(inputs, '_'));

if param.vaxnum > 1
    dirnm{param.vaxnum-1} = fullfile('..', data_locs{param.vaxnum-1}, ...
        strjoin([num2str(param.vaxnum-1), inputs(2:end-1), num2str(vaxtiming(param.vaxnum))], '_'));
end

for i = 1:length(dirnm)
    fnm{i} = fullfile(dirnm{i}, [num2str(param.first), '_to_', num2str(param.last), '.mat']);
end

% Display constructed input string
fprintf('%s\n', strjoin(inputs, '_'));
end