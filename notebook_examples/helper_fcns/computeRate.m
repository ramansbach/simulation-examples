function rate = computeRate(a,i,cSizes)
%for a cluster of size ai+ai, compute the rate at which clusters of size
%ai+ai form from clusters of size ai and ai
%find molecules in clusters of size ai and see how many times they go to
%molecules in clusters of size 2ai, then divide by the size of the larger
%cluster because you've overcounted by that much

rate = 0;
sz = size(cSizes);
cSizes0 = [zeros(1,sz(2));cSizes(1:(end-1),:)];
cDiffs = cSizes - cSizes0;
lcDiffs = cDiffs == (a*i);
lcSizes = cSizes == 2*a*i;
lInds = lcDiffs.*lcSizes;
trate = sum(sum(lInds))/(2*a*i);
rate = trate/sz(1);