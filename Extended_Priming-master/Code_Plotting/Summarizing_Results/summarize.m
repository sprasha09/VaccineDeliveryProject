function [tspan, totalnum, agconc, abtiter, agconc_mean, abtiter_mean, IgM_mean, IgG_mean, gcnum] = summarize(result)
% This function summarizes the results.
%
% Input:
%   result : The results structure containing param, gc, and conc data.
%            Obtained from concatenating repeat simulations.
%
% Output:
%   tspan       : Time span summary from result.
%   totalnum    : Total number of GC B cells over time.
%   agconc_mean : Mean Antigen concentration.
%   abtiter_mean : Mean Antibody titer.
%   IgM_mean    : Mean IgM titer.
%   IgG_mean    : Mean IgG titer.
%   gcnum       : GC B cell numbers in individual GCs.

tspan = result.param.tspan_summary;
totalnum = squeeze(sum(result.gc.numbytime, 1));
gcnum = result.gc.numbytime;
[agconc, agconc_mean] = getAgConc(result.param, result);
[abtiter, abtiter_mean, IgM_mean, IgG_mean] = getAbTiter(result);
end

%% Subfunctions

function [agconc, agconc_mean] = getAgConc(param, result)
% Calculates the mean aggregate concentration from a set of simulation results.
% Output:
%   - agconc_mean: Mean aggregate concentration as a 4xN matrix, where N is the number of time points.
Ag0 = 10;
agconc = cell(1,length(result.conc));
agconc_mean = zeros(4,length(param.tspan_summary));
for i=1:length(result.conc)
    agconc{i} = squeeze(result.conc{i}.concarray(1,:,:))/Ag0;
    
    agconc_mean = agconc_mean + agconc{i};
    agconc_mean(isnan(agconc_mean)) = 0;
end
agconc_mean = agconc_mean/length(result.conc);


end



function [abtiter, abtiter_mean, IgM_mean, IgG_mean] = getAbTiter(result)
% Calculates the mean antibody titers and antibody class concentrations from simulation results.
%
% Output:
%   - abtiter_mean: Mean antibody titers as a 2xN matrix, where N is the number of time points.
%   - IgM_mean: Mean IgM class antibody concentrations as a 2xN matrix.
%   - IgG_mean: Mean IgG class antibody concentrations as a 2xN matrix.

abtiter = cell(1,length(result.conc));
IgM = cell(1,length(result.conc));
IgG = cell(1,length(result.conc));
    for i=1:length(result.conc)
    abtiter{i} = ...
    squeeze(sum(result.conc{i}.concarray(3:4,3:4,:).*result.conc{i}.Kaarray(2:3,:,:)));
    IgM{i} = squeeze((result.conc{i}.concarray(3,3:4,:).*result.conc{i}.Kaarray(2,:,:)));
    IgG{i} = squeeze((result.conc{i}.concarray(4,3:4,:).*result.conc{i}.Kaarray(3,:,:)));
    end
abtiter_mean = cell2mat(abtiter');
abtiter_mean = [mean(abtiter_mean(1:2:end-1,:)); mean(abtiter_mean(2:2:end,:))];
IgM_mean = cell2mat(IgM');
IgM_mean = [mean(IgM_mean(1:2:end-1,:)); mean(IgM_mean(2:2:end,:))];
IgG_mean = cell2mat(IgG');
IgG_mean = [mean(IgG_mean(1:2:end-1,:)); mean(IgG_mean(2:2:end,:))];
end