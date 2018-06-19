function errs = rms_error(prof,bprofs)
%Takes the original profile and a matrix of profiles for smaller distances,
%calls interp_rms to get back the rms versions of the profiles, and then
%finds the error in the original profile as the std dev of the values of
%the corresponding interpolated blocks
rmsprofs = zeros(size(prof,1),size(bprofs,2)/2);
figure(2);
plot(prof(:,1),prof(:,2));
hold on;
for i=1:2:size(bprofs,2)
    [hshift,rmsdiff,rmsprof,inds] = interp_rms(prof,bprofs(:,i:i+1));
    rmsprofs(inds,(i-1)/2+1) = rmsprof(inds,2);
    plot(rmsprof(inds,1),rmsprof(inds,2));
end
errs = std(rmsprofs,0,2);