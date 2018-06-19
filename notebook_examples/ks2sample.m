function [D,a] = ks2sample(cdf1,cdf2)
%function that computes the 2-sample Kolmogorov-Smirnov test statistic for
%comparison to whatever significance level you like (see also 
%https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test, the section
%on the 2-sample test
%also returns the value to multiply c(alpha) by such that
%D > a*c(alpha) represents rejection of the null hypothesis and c(alpha) is
%a function that depends on the significance level
%null hypothesis is that the two samples are drawn from the same
%distribution
D = max(abs(cdf1-cdf2));
n = max(size(cdf1,1),size(cdf1,2));
np = max(size(cdf2,1),size(cdf2,2));
a = sqrt((n+np)/(n*np));