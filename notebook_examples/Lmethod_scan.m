function [splts,rses] = Lmethod_scan(x,y)
mx = max(y);
steps = length(mx:-1:1);
splts = zeros(1,steps);
rses = zeros(1,steps);
for i=mx:-1:3
    [~,~,splt,rse,rse2] = Lmethod(log(x(y<=i)),log(y(y<=i)));
    splts(i) = splt;
    rses(i) = rse+rse2;
end