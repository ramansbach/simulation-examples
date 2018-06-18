function cldistrib = clDistrib(ldistrib,cids,maxlen)
%given a row of the lengths of the cluster each molecule is involved in and
%a row of the cluster ID that each molecule is involved in, return a
%set of unique cluster lengths

ucids = unique(cids);
if ~exist('maxlen','var')
    maxlen = max(ldistrib)+1;
end

nclusts = length(ucids);
cldistrib = zeros(1,nclusts);

for i = 1:nclusts
   cid = ucids(i);
   ls = ldistrib(cids==cid);
   lcurr = ls(1);
   if lcurr < maxlen %remove clusters with a length over the maximum
       cldistrib(i) = ls(1); 
   else
       cldistrib(i) = NaN;
   end
end
cldistrib = cldistrib(~isnan(cldistrib)); %remove clusters with a length over the maximum