function [time,conc] = getInnateDynamics(T, numshot, k, kinetic_params, T_slow)
%---------------------------------------
% Summary:
% This function calls getInnateDynamicsFromDose.
% From the dosing scheme parameters T, numshot, and k, calculate the
% exact doses and timing of antigen and adjuvant and then pass them
% on to getInnateDynamicsFromDose. 

% Inputs: T, numshot, k are scalars that define the dosing scheme.
%   See also the documentation for getInnateDynamicsFromDose. 
% Outputs: time and concentration arrays
%---------------------------------------
rng(1)

if nargin < 5
    T_slow = []; % Not extended dosing
end

% Calculate the dose and timing
dose_t = cell(1,2); dose = cell(1,2); %Antigen and Adjuvant
for i=1:2
    dose_t{i} = [round(linspace(0,T(i),numshot(i)))];
    dose{i} = exp((0:(numshot(i)-1))*k(i))...
        /sum(exp((0:(numshot(i)-1))*k(i)));
end

% Need to handle the case when slow doses overlap - the next dose is given
% before the release is over.



% For the case when Ag and Adj injection timing are different, we want to
% combine them to get one overall dosing scheme. 
if ~isequal(dose_t{1}, dose_t{2}) || ~isequal(dose{1}, dose{2})
    t = union(dose_t{1}, dose_t{2});
    t_intersect = intersect(dose_t{1},dose_t{2});
    dose_1 = cell(1,2);
    for i=1:length(t)
        if ismember(t(i), t_intersect)
            dose_1{1}(i) = dose{1}(t(i)==dose_t{1});
            dose_1{2}(i) = dose{2}(t(i)==dose_t{2});
        elseif ismember(t(i), dose_t{1})
            dose_1{1}(i) = dose{1}(t(i)==dose_t{1});
            dose_1{2}(i) = 0;
        else
            dose_1{1}(i) = 0;
            dose_1{2}(i) = dose{2}(t(i)==dose_t{2});        
        end
    end
    dose = dose_1;
    dose_t{1} = t;
    dose_t{2} = t;
end

% Call getInnateDynamicsFromDose using the dosing scheme
[time, conc] = getInnateDynamicsFromDose(dose, dose_t, kinetic_params, T_slow);
end

%#######################################################


%% Subfunctions
function [time,conc] = getInnateDynamicsFromDose( ...
                        dose, dose_t, kinetic_params, T_slow)
%------------------------------------%
% Inputs: Antigen and adjuvant dose and timing, kinetic parameters, and
%   T_slow which is non-empty if any of the doses are slowly released.
%   In this case, T_slow is an m x 2 array, each row representing the
%   duration over which antigen and adjuvant are released. E.g. [5,1] means
%   that antigen is released over 5 days and adjuvant is released over 1 day. 
%   
% Outputs: Time and concentration arrays. 
% Summary: First, define multiple time intervals based on the number of
%   doses. Then sequentially solve the ODEs for the intervals. Each 
%   subsequent interval takes in as the initial condition the output of
%   the previous interval. 
% NOTE:
% Current version of the code assumes that dose_t is identical for antigen
% and adjuvant. 2023/3/7 Leerang Yang %
%------------------------------------%

% Deal with optional parameter T_slow
if nargin < 4
    T_slow = [];
elseif ~any(T_slow)
    T_slow = [];
end
    

% Define time intervals for which ODEs will be solved
tmp = dose_t{1};
dose_t{1} = [tmp, tmp(end)+30]; % Antigen. Add final time point
dose_t{2} = [tmp, tmp(end)+30]; % Adjuvant. Add final time point
clear tmp

% Initial condition for the concentrations
y0 = zeros(1,8);
%      1-Ag, 2-Adj, 3-Sentinel, 4-DC, 5-DC-Adj, 6-DC-Adj-Ag, 
%      7-T cell (T-B border), 8-T cell (Tfh)
y0(7) = kinetic_params(end); % T0 - controls the scale of the model

% Result storage arrays. Not assigning pre-defined size because the number 
% of points from solving ODE is unknown. There are only a few intervals,
% so this is not computationally too problematic. 
time = 0;
conc = y0;

% Solve ODE from each interval between doses and concatenate the results.
num_slow_dose = size(T_slow,1); % The last n shots are extended release
num_total_dose = length(dose_t{1})-1;

for i=1:num_total_dose % Iterate over the intervals
    if i > (num_total_dose-num_slow_dose) % Current dose is slowly released
        k = num_slow_dose-(num_total_dose-i); %index of slow-relase dose
        slow_final = dose_t{1}(i) + T_slow(k,:);
        slow_rates = dose{1}(i)./T_slow(k,:);
        if T_slow(k,2)==0                 % Bolus adjuvant
            y0(2) = y0(2) + dose{2}(i);
        end
    else                                  % Current dose is not slowly released
        y0(1) = y0(1) + dose{1}(i);
        y0(2) = y0(2) + dose{2}(i);
        slow_final = [];
        slow_rates = [];
    end
    tspan = dose_t{1}(i):0.01:dose_t{1}(i+1); % Use same interval as the B cell simulation
    [t,y] = ode45(@(t,y) odeInnate( ...       % Solve the ODE
           t, y, kinetic_params, slow_final, slow_rates), tspan, y0);
    y0 = y(end,:);

    % Concatenate the time and conc calculated from each interval 
    time = [time;t(2:end)];
    conc = [conc;y(2:end,:)];
end
end

%#######################################################

function dydt = odeInnate(t, y, kinetic_params, slow_final, slow_rates)
%---------------------------------------
% Inputs: t: time, y: vector of concentration of the species
%      1-Ag, 2-Adj, 3-Sentinel, 4-DC, 5-DC-Adj, 6-DC-Adj-Ag, 
%      7-T cell (T-B border), 8-T cell (Tfh)
%      k: cell array containing all kinetic parameters
%      slow_final: When non-empty, the odse is given as slow 
%         release. slow_final is 1x2 array, defining the durations
%         over which Ag and Adj are released. The values correspond to the
%         date when the slow-release ends (e.g. [12,8] means the final
%         times of slow-release for antigen and adjuvant are d12 and d8.)
%      slow_rates: defines the amounts of Ag/Adj given slowly
% Outputs: First-order time derivative of y
% Summary: Define the system of ODEs that govern the dynamics of 
%      antigen, adjuvant, innate immune cells, and T cell.
%---------------------------------------

% Initialization
dydt = zeros(8,1);
kinetic_params = num2cell(kinetic_params);
[d_Ag, d_Adj, K_Adj, alpha, delta, mu, k1, C0, T0]=deal(kinetic_params{:});

% Ag and Adj dynamics
d = [d_Ag, d_Adj];
for j=1:2
    if ~isempty(slow_rates) && t<slow_final(j) 
        dydt(j) = -d(j)*y(j) + slow_rates(j);
    else
        dydt(j) = -d(j)*y(j);
    end
end

% Sentinel cells dynamics
dydt(3) = y(2)/(K_Adj+y(2)) - mu*y(3);

% DC dynamics
dydt(4) = C0*y(3) - (mu + (1+k1*y(2)/(y(2)+K_Adj))*y(1))*y(4); % DCs
dydt(5) = 0; % DC-Adj
dydt(6) = (1+k1*y(2)/(y(2)+K_Adj))*y(1)*y(4) - mu*y(6); % DC-ag

% T cell dynamics
if y(7)>=T0
    dydt(7) = alpha*(y(6)/(y(6)+y(7)))*y(7) - delta*(y(7)-T0);
else
    dydt(7) = T0;
end
dydt(8) = delta*(y(7)-T0); 
end