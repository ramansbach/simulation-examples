function [sigmaxs,f] = computeSigMax(dims,Ofracs,i)
sigmaxs = zeros(1,12);
for j = 1:12
    ind = find(Ofracs(((j-1)*4+1):(4*j),i)==max(Ofracs(((j-1)*4+1):(4*j),i)));
    sigmaxs(j) = dims((j-1)*4+ind,3);
end
f = figure();
set(gca,'fontsize',20)
hold on
eAs = [2.5 2.5 2.5 2.5 5 5 5 5 7.5 7.5 7.5 7.5];
eSCs = [0.2 2 6 10 0.2 2 6 10 0.2 2 6 10];
f = errorMesh(eAs,eSCs,sigmaxs,zeros(size(sigmaxs)),[0 0 0],f,0.);
xlabel('\epsilon_A')
ylabel('\epsilon_{SC}')
zlabel('\sigma_{SC}')
colorbar