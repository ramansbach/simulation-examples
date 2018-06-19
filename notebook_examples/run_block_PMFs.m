function [prof,errs] = run_block_PMFs(wdir,tstart,tinterval,tend,remove_entropic,minmaxx)
%this is a helper function for prepping and doing dimer PMF blocking for
%rotationally-restrained PMF simulations for patchy particle doing
prof = importdata([wdir '/prof.xvg']); %full profile
prof = prof.data;
intprofs = zeros(length(prof),2*(tend-tstart)/tinterval);
ind = 1;
for i = tstart:tinterval:tend-1
   currf = [wdir '/prof' int2str(i) '-' int2str(i+tinterval) '.xvg']; 
   data = importdata(currf);
   intprofs(:,ind:ind+1) = data.data;
   ind = ind+2;
end
errs = fed_error(prof,intprofs,remove_entropic,minmaxx); 
figure()
hold on
set(gca,'fontsize',20)
if remove_entropic
   prof(:,2) = prof(:,2)+2*log(prof(:,1)); 
end
errorbar(prof(:,1),prof(:,2)-min(prof(:,2)),errs)
plot(prof(:,1),prof(:,2)-min(prof(:,2)),'--','linewidth',2)
