function [chi2,p,nu] = chi22(Rs,Ss)
%function that takes two histograms of counts with the same length,
%and returns the chi-squared value and the corresponding p-value
%citation: Numerical Recipes in C: The Art of Scientific Computing
%null hypothesis: R & S are drawn from the same distribution
%large chi2 value means this is unlikely
%if p < significance, reject the null hypothesis (significance is typically
%~ 0.05)
R = sum(Rs);
S = sum(Ss);
chi2 = 0;
nu = length(Rs); %currently assuming different different numbers of counts and normalizations
for i = 1:nu
    if Rs(i) ~= 0 || Ss(i) ~= 0
        chi2 = chi2 + (sqrt(S/R)*Rs(i)-sqrt(R/S)*Ss(i))^2/(Rs(i)+Ss(i)); 
    end
end
p = 1- (1 - gammainc(nu/2,chi2/2));