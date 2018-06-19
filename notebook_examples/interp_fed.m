function [hshift,feddiff,rmsfit,inds] = interp_fed(prof1,prof2,lg,minmaxx)
% Takes an original (PMF) profile and another profile
% First linear-interpolate prof2 to have the same x-values as prof1
% then align them using rms difference of a portion of their ends
% lg toggles whether to adjust for entropic contribution
if lg

   prof1(:,2) = prof1(:,2) + 2*log(prof1(:,1));
   
   prof2(:,2) = prof2(:,2) + 2*log(prof2(:,1));
end
if ~exist('minmaxx','var')
   minmaxx = []; 
end
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
if isempty(minmaxx)
    figure();
    hold on;
    plot(prof1(inds,1),prof1(inds,2));
    plot(interp_prof(inds,1),interp_prof(inds,2));

    mn = input('What is the approximate minimum x-value to use for RMS-fitting? ');
    mx = input('What is the approximate maximum x-value to use for RMS-fitting? ');
    close;
else
    mn = minmaxx(1);
    mx = minmaxx(2);
end
minind = find(prof1(inds,1)>mn, 1 );
maxind = find(prof1(inds,1)<mx, 1, 'last' );
hshift = sum(prof1(minind:maxind,2)-interp_prof(minind:maxind,2))/(maxind-minind);

rmsfit = interp_prof;
rmsfit(:,2) = rmsfit(:,2)+hshift;
feddiff = (1/size(prof1,1))*sqrt(sum((prof1(inds,2)-rmsfit(inds,2)).^2));
%rmsdiff = (1/size(prof1,1))*sqrt(sum((prof1(inds,2)-rmsfit(:,2)).^2)); %rms difference per point
%rmsfit = rmsfit(inds,:);