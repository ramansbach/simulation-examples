ind = 1;
for dg1plus1 = dg1plus1s
for dgnplus1 = dgnplus1s
x = douts(ind,:);
mouts(ind,6) = (x(end)/x(end-1)>0.1);
ind = ind + 1;
end
end