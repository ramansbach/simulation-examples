function [hshift,rmsdiff,rmsfit,inds] = interp_rms(prof1,prof2)
% Takes an original (PMF) profile and another profile
% First linear-interpolate prof2 to have the same x-values as prof1
% then find the 
interp_prof = zeros(size(prof1));
interp_prof(:,1) = prof1(:,1);
useinds = ones(size(prof1,1),1);
for i=1:size(interp_prof,1)
    
    low = find(prof2(:,1)<prof1(i,1));
    high = find(prof2(:,1)>prof1(i,1));
    if (~isempty(low))&&(~isempty(high))
        low = low(end);
        high = high(1);
        interp_prof(i,2) = ((prof2(high,2)-prof2(low,2))/(prof2(high,1)-prof2(low,1)))*(prof2(high,1)-prof1(i,1))+prof2(low,2);
    else
        useinds(i) = 0;
    end
end
inds = find(useinds);
hshift = sum(prof1(inds,2)-interp_prof(inds,2))/size(prof1,1);

rmsfit = interp_prof;
rmsfit(:,2) = rmsfit(:,2)+hshift;
rmsdiff = (1/size(inds,1))*sqrt(sum((prof1(inds,2)-rmsfit(inds,2)).^2)); %rms difference per point
%rmsfit = rmsfit(inds,:);