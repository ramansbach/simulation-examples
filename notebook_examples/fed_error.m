function errs = fed_error(prof,bprofs,lg,minmaxx)
%Takes the original profile and a matrix of profiles for smaller distances,
%calls interp_rms to get back the rms versions of the profiles, and then
%finds the error in the original profile as the std dev of the values of
%the corresponding interpolated blocks
%lg controls whether to take out the entropy contribution of the PMFs
rmsprofs = zeros(size(prof,1),size(bprofs,2)/2);
figure();
if lg
    plot(prof(:,1),prof(:,2)+2*log(prof(:,1)));
else
    plot(prof(:,1),prof(:,2));
end
hold on;
for i=1:2:size(bprofs,2)
    if exist('minmaxx','var')
        [~,~,rmsprof,inds] = interp_fed(prof,bprofs(:,i:i+1),lg,minmaxx);
    else
        [~,~,rmsprof,inds] = interp_fed(prof,bprofs(:,i:i+1),lg);
    end
    rmsprofs(inds,(i-1)/2+1) = rmsprof(inds,2);
    plot(rmsprof(inds,1),rmsprof(inds,2));
end
errs = std(rmsprofs,0,2);