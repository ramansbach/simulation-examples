function [l1,l2,linf] = lMoments(lMatrix,cidMatrix)
%takes matrix of the lengths and cluster IDs of the occupying molecules
%over time and computes the first, second, and infinite moment averages 
sz = size(lMatrix);
ttotal = sz(1);
l1 = zeros(ttotal,1);
l2 = zeros(ttotal,1);
linf = zeros(ttotal,1);
for t = 1:ttotal
    lid = clDistrib(lMatrix(t,:),cidMatrix(t,:));
    if sum(~isnan(lid)) > 0
        l1(t) = mean(lid(~isnan(lid)));
        l2(t) = sum(lid(~isnan(lid)).^2)/sum(lid(~isnan(lid)));
        
        linf(t) = max(lid(~isnan(lid)));

    else
       l1(t) = NaN;
       l2(t) = NaN;
       linf(t) = NaN;
    end
end