%performs chi-squared test on a matrix whose rows are different COUNTS
%(non-normalized histograms from histcounts)
function [chi2s,ps] = chiMat(mh)

chi2s = zeros(size(mh,1),size(mh,1));
ps = zeros(size(mh,1),size(mh,1));
for i = 1:size(mh,1)
   for j = 1:size(mh,1)
       [chi2,p,nu] = chi22(mh(i,:),mh(j,:));
       chi2s(i,j) = chi2;
       ps(i,j) = p;
   end
end