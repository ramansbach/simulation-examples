function [sigma,ssigma] = computeSasaSig(fnames)
% Function that takes a (cell) list of file names, opens the files, 
% reads in the time course of the solvent accessible surface area for each
% one and returns the LJ sigma of the equivalent side chain made up of all
% the different residues in the file names
Vtot = 0;
sVtot = 0;
for i = 1:length(fnames)
   fname = fnames{i};
   sasa = importdata(fname,' ',23);
   sasa = sasa.data(:,2);
   msasa = mean(sasa);
   smsasa = std(sasa);
   rsasa = sqrt(msasa/(4*pi));
   srsasa = (smsasa/msasa)*rsasa;
   Vsasa =  (4/3)*pi*rsasa^3;
   sVsasa = (srsasa/rsasa)*Vsasa;
   Vtot = Vtot + Vsasa;
   sVtot = sVtot + sVsasa^2;
end
sVtot = sqrt(sVtot);
sigma = 2*(Vtot*(3/(4*pi)))^(1/3);
ssigma = 2*(sVtot/Vtot)*sigma;