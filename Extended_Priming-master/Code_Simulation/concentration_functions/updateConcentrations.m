function [agconc, abconc, Ka, Ka_var] = updateConcentrations(...
    agconc, abconc, Ka, Ka_var, plasmaBlasts, plasmaCells, plasmaCellsEGC, t, param)
%% Documentation
% Summary:
% Updates the antigen and antibody concentrations, antibody affinities

% Outputs: 
%   agconc: 1x3 vector; Concentrations of Antigens; 
%                       soluble, IC1, IC2
%   abconc: 3x2 vector; Concentration of antibodies
%           Dim1 - IgM(natural),IgM(immune),IgG; Dim2-Epitopes
%   Ka: 3x2 vector; Antibody binding affinities (Ka)
%   param: parameter struct
%
%     Note:
%     IC-FDC1: Deposited to FDC by Ig that targets dominant epitope
%     IC-FDC2: Deposited to FDC by Ig that targets subdominant epitope
%     However in this study there is no difference between the two

% Inputs: 
%   agconc, abconc, Ka, Ka_var are same as the outputs
%   plasmaBlasts, plasmaCells, plasmaCellsEGC are arrays of plasma cells
%   t: current time
%   param: parameter struct


%% Function definition
% Check if antigen should be given at current time 
if any(t==param.dose_t)
    if param.pSER==0 || (param.pSER==1 && t~=param.dose_t(end))
        agconc(1) = agconc(1) + param.dose(t==param.dose_t);
    end
end

% Get average affinity of the serum, and soluble Ag and Ab concentrations
R = agconc(1:2); L = sum(abconc);
Ka_avg = sum(Ka.*abconc);
Ka_avg(L>0) = Ka_avg(L>0)./L(L>0);

% Get rates for deposition of IC onto FDC
r = zeros(size(Ka));
IC_soluble = zeros(1,2);
for i=1:2
    if L(i)>0 && R(i)>0
        IC_soluble(i) = ((R(i)+L(i)+1/Ka_avg(i))  ... 
            - sqrt((R(i)+L(i)+1/Ka_avg(i))^2-4*R(i)*L(i))); %Equil.
        r(:,i) = param.k_deposit*IC_soluble(i)*(Ka(:,i) ... 
            .*abconc(:,i))/sum(sum(Ka(:,i).*abconc(:,i)));  
            % 3x2 array; rates of IC-FDC deposition by the Ab species
        if ~isreal(IC_soluble(i)) || isnan(IC_soluble(i)) % Check if im or nan
            error('imaginery number or nan for soluble IC concentration') 
        end
    end
end

abconc_old = abconc; %Temporarily save current Ab concentrations

% Check if any of the antibody species concentration will go to zero.
% If this is going to happen, then rescale the reaction rates.
decay_Ab = -r - [0;param.d_IgM;param.d_IgG].*abconc;
rescale_idx = (abconc<-decay_Ab*param.dt);
rescale_factor = abconc./(-decay_Ab*param.dt);
if any(rescale_idx(:))
    fprintf('Reaction rates rescaled. Time t = %.2f', t)
    r(rescale_idx) = r(rescale_idx).*rescale_factor(rescale_idx);
    decay_Ab(rescale_idx) = decay_Ab(rescale_idx).*rescale_factor(rescale_idx);
    if any(isnan(decay_Ab))
        error('Rescaled Ab decay rates contain NaN')
    end
end

% Check if the soluble Ag concentration will go to 0. Rescale if needed.
decay_Ag = zeros(1,2);
breakdown_Ag = zeros(1,2);
for i=1:2
decay_Ag(i) =  - param.d_Ag*agconc(i)- sum(r(:,i)); % Net consumption of soluble Ag
breakdown_Ag(i) = - param.d_Ag*agconc(i);
if agconc(i)<(-decay_Ag(i)*param.dt)
    rescale_factor = agconc(i)/(-decay_Ag(i)*param.dt);
    decay_Ag(i) = decay_Ag(i)*rescale_factor;
    breakdown_Ag(i) = breakdown_Ag(i)*rescale_factor;
    rnew = r(:,i)*rescale_factor;
    decay_Ab(:,i) = decay_Ab(:,i) + r(:,i) - rnew;
    r(:,i) = rnew;
end
end

%Update the antigen & antibody concentrations for injection / consumption
    agconc(1) = agconc(1) + decay_Ag(1)*param.dt; %+ param.F0*exp(param.k*t)*(t<param.T)*param.dt;
    if param.pSER==1 && t>=param.dose_t(end) && t<param.dose_t(end)+param.pSER_ag_release_time
        agconc(1) = agconc(1) + param.dose(end)/param.pSER_ag_release_time*param.dt;
    elseif param.pSER==2
        ag_release = find(t>=param.dose_t & t<param.dose_t+param.pSER_ag_release_time);
        if ~isempty(ag_release)
            agconc(1) = agconc(1) + sum(param.dose(ag_release))/param.pSER_ag_release_time*param.dt;
        end
    end
    agconc(2) = agconc(2) + decay_Ag(2)*param.dt - breakdown_Ag(1)*param.dt;
    agconc(3:end) = agconc(3:end) + sum(r,1)*param.dt - agconc(3:end)*param.d_IC*param.dt;
    
    abconc = abconc + decay_Ab*param.dt;
    abconc(abs(abconc)<10^-10)=0;
    if any(agconc(:)<0) || any(abconc(:)<0)
       error('Negative concentration')
    end
        
    
%Update the antibody concentration & affinity for production
M_GC = param.M_GC;
Ig_new = zeros(3,param.n_ep);
Ka_new = repmat({Ig_new}, 1, 2); 
Affinity = cell(2,3);
Target = cell(1,3);

for i=1:2 %WT and Variant
    Affinity{i,1} = plasmaBlasts((i+1)*M_GC+1:(i+2)*M_GC,:);
    Affinity{i,2} = plasmaCells((i+1)*M_GC+1:(i+2)*M_GC,:);
    Affinity{i,3} = plasmaCellsEGC(i+2,:);
    Target{1} = plasmaBlasts(M_GC+1:2*M_GC,:) .* (plasmaBlasts(4*M_GC+1:5*M_GC,:)<(t-param.delay)) .* (plasmaBlasts(4*M_GC+1:5*M_GC,:)>(t-3));
    Target{2} = plasmaCells(M_GC+1:2*M_GC,:) .* (plasmaCells(5*M_GC+1:6*M_GC,:)<(t-param.delay));
    Target{3} = plasmaCellsEGC(2,:); % .* (plasmaCellsEGC(6,:)<(t-param.delay));
end

for Ig_type = 1:3 %IgM, IgG-GCPC, IgG-EGCPC
    for target=1:param.n_ep
        if any(Target{Ig_type}(:)==target)
            % Ig production
            Ig_new(Ig_type,target) = sum(sum(Target{Ig_type}==target)) * ...
                                (param.r_IgM*(Ig_type==1)+param.r_IgG*(Ig_type==2)+param.r_IgG_EGC*(Ig_type==3));
            % Ka of new Ig
            for variant=1:2
               Ka_new{variant}(Ig_type,target) = ...
                   mean(10.^(Affinity{variant,Ig_type}(Target{Ig_type}==target)-9));
            end
        end
    end
end

%Update amounts
abconc(2:3, :) = abconc(2:3, :) + [Ig_new(1,:); sum(Ig_new(2:3, :))]*param.dt;

%update Ka
Ka = {Ka, Ka_var};
for variant=1:2
    current_sum = (abconc_old + param.dt*decay_Ab).*Ka{variant};
    new_sum = current_sum + [0,0,0;1,0,0;0,1,1]*(Ig_new.*Ka_new{variant}*param.dt);
    new_Ka = new_sum./abconc;
    new_Ka(abconc==0)=0;
    if any(new_Ka(:)<0) || any(abs(new_Ka(:))>10^11)
        warning('Error in Ka value: negative or too large')
    end
    new_Ka(isnan(new_Ka)) = 0;
    Ka{variant} = new_Ka;
end
Ka_var = Ka{2};
Ka = Ka{1};
epsilon = 10^-10;
abconc(abconc<epsilon) = 0;
agconc(agconc<epsilon) = 0;
end