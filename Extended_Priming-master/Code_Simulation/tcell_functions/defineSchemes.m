
%% Subfunctions -- Parameter Fitting

function scheme = defineSchemes()
% Define the dosing schemes, colors, and names.
scheme.T = {[12,12],[11,11],[12,12],[12,12],[7,7],[0,0],[12,0]}; 
scheme.numshot = {[7,7],[6,6],[4,4],[3,3],[2,2],[1,1],[7,1]}; 
scheme.k = {[1,1],[1,1],[1,1],[1,1],[log(4),log(4)],[0,0],[1,0]};

scheme.colors = {[0,0,0], [253, 128, 8]/256, [204, 102, 255]/256, [64, 0, 128]/256, ...
    [251, 2, 128]/256, [128, 64, 3]/256, [1,0,0]};

scheme.names = {'7-ED', '6-ED', '4-ED', '3-ED', '2-ED', 'Bolus', 'Adjuvant Bolus'};
scheme.names_low = {'7-ED Low', '6-ED Low ', '4-ED Low', '3-ED Low', '2-ED Low', 'Bolus Low', 'Adjuvant Bolus Low'};
scheme.names_high = {'7-ED High', '6-ED High', '4-ED High', '3-ED High', '2-ED High', 'Bolus High', 'Adjuvant Bolus High'};
end

%#######################################################
