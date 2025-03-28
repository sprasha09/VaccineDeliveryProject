function [numActiveGCs, avgNumBcells] = getGCStat(gcnum, idx)
    numB = squeeze(gcnum{idx}(:,113,1))+squeeze(gcnum{idx}(:,226,1));
    % The line above is giving an error (gcnum does not have enough input
    % arguments)
    avgNumBcells = mean(numB(numB>0));
    numActiveGCs = sum(sum(squeeze(gcnum{idx}(:,:,1)),2)>0);
end
