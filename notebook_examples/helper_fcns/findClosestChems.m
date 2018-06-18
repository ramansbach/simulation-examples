function [chems,sasas,pmfs,bothInds] = findClosestChems(chemDict,sasaTable,pmfTable,sasaTarget,pmfTarget)
% function that returns all chemistries (and their corresponding sasas and
% pmfs) that fall within a target sigma LJ range and with ranges
% encompassing the target pmf
%
% ----------
% Parameters
% ----------
% chemDict: 1D cell array of possible chemistries, length C
% sasaTable: C x 2 array of projected LJ sigma range from mean - err to
% mean + err
% pmfTable: C x 2 array of projected PMFs for NDI and PDI cores
% sasaTargetRange: 1 x 2 vector that defines a target range for the LJ
% sigma parameter
% pmfTarget: target value for LJ well depth
%
% -------
% Returns
% -------
% chems: cell array of candidate chemistries
% sasas: vector of corresponding sigma_SC's
% pmfs: vector of corresponding eps_SC's

sasaInds = ((sasaTarget - sasaTable(:,1)) > 0) & ((sasaTable(:,2) - sasaTarget) > 0);
pmfInds = ((pmfTarget - pmfTable(:,1)) > 0 ) & ((pmfTable(:,2) - pmfTarget) > 0);
bothInds = sasaInds.*pmfInds;
chems = chemDict(find(bothInds));
sasas = sasaTable(find(bothInds),:);
pmfs = pmfTable(find(bothInds),:);