function [l1,l2,splt,lsing,aic1,aic2,fig1,fig2] = Lmethod(n,vals,weights,weightbool,errsy)
%fit two lines to the resulting graph and minimize the error
%allow for weighted fit
%actually input std devs not the weights themselves
%{(1,2),(2,l)}...{(1,l-1),(l-1,l)}
l = length(unique(n));
%n = transpose(1:l);
minerr = Inf;
evecn = 0;
fit1 = [];
fit2 = [];
x1 = [];
x2 = [];
uns = unique(n);
lf = fittype('a*x+b','independent','x');
if size(n,1) < size(n,2)
   n = transpose(n); 
end
if size(vals,1) < size(vals,2)
   vals = transpose(vals); 
end
if exist('weights','var')
   if size(weights,1) < size(weights,2)
      weights = transpose(weights); 
   end
end
if ~exist('errsy','var')
   errsy = 0; 
end
for ind = 2:l-1
    %jinds = find(n==uns(ind));
    %j = jinds(1);
    j = ind;
    %P1 = polyfit(n(1:j),vals(1:j),1);
    %P2 = polyfit(n(j:end),vals(j:end),1);
    if weightbool
        %perform weighted fit
        ws1 = 1./(weights(1:j).^2);
        ws2 = 1./(weights(j:end).^2);

        [P1,gP1] = fit(n(1:j),vals(1:j),lf,'Weights',ws1);
        [P2,gP2] = fit(n(j:end),vals(j:end),lf,'Weights',ws2);
    else
        [P1,gP1] = fit(n(1:j),vals(1:j),lf);
        try
        [P2,gP2] = fit(n(j:end),vals(j:end),lf);
        catch
           j; 
        end
    end
    %P1fit = polyval(P1,n(1:j));
    %P2fit = polyval(P2,n(j:end));
    P1fit = P1.a*n(1:j)+P1.b;
    P2fit = P2.a*n(j:end)+P2.b;
    res1 = P1fit - vals(1:j);
    res2 = P2fit - vals(j:end);
    RMSE1 = sqrt(sum(res1.^2)/length(res1));
    RMSE2 = sqrt(sum(res2.^2)/length(res2));
    RMSE1 = gP1.rmse;
    RMSE2 = gP2.rmse;
    RMSET = ((j-1)/(length(vals)-1))*RMSE1 + ((length(vals)-j)/(length(vals)-1))*RMSE2;
    if RMSET < minerr
        minerr = RMSET;
        evecn = j;
        fit1 = P1fit;
        l1 = P1;
        l2 = P2;
        gl1 = gP1;
        gl2 = gP2;
        fit2 = P2fit;
        x1 = n(1:j);
        x2 = n(j:end);
        splt = j;
    end
end
%P = polyfit(n,vals,1);
if weightbool
    [P,gP] = fit(n,vals,lf,'Weights',1./(weights.^2));
else
    [P,gP] = fit(n,vals,lf);
end
%Pfit = polyval(P,n);
lsing = P;
Pfit = lsing.a*n+lsing.b;
rse = Pfit - vals;
P1fit = l1.a*n(1:splt)+l1.b;
P2fit = l2.a*n(splt:end)+l2.b;
res1 = P1fit-vals(1:splt);
res2 = P2fit-vals(splt:end);
rse1 = sum(res1.^2);
rse2 = sum(res2.^2);
%rmse1 = sqrt(sum(res1.^2)/length(res1))
%rmse2 = sqrt(sum(res2.^2)/length(res2))
%rmse = sqrt(sum(rse.^2)/length(rse))
rmse1 = gl1.rmse;
rmse2 = gl2.rmse;
rmse = gP.rmse;
rsesing = sum(rse.^2);
%lsing = P;
aic1 = akaike(3,l,rsesing);
aic2 = akaike(5,l,rse1+rse2);
fig1 = figure();
hold on
ax1 = gca;
fig2 = figure();
hold on
ax2 = gca;
cint1 = confint(l1);
cint2 = confint(l2);
cint = confint(lsing);
clow = cint(1,1)*n+cint(1,2);
chigh = cint(2,1)*n+cint(2,2);
clow1 = cint1(1,1)*x1+cint1(1,2);
chigh1 = cint1(2,1)*x1+cint1(2,2);
clow2 = cint2(1,1)*x2+cint2(1,2);
chigh2 = cint2(2,1)*x2+cint2(2,2);
if weightbool
    if errsy
    errorbar(ax1,n,vals,weights,'.','MarkerSize',15)
    errorbar(ax2,n,vals,weights,'.','MarkerSize',15)
    else
    errorbarxy(ax1,n,vals,weights,[],{'.',[0 0.4470 0.7410],'r'})
    errorbarxy(ax2,n,vals,weights,[],{'.',[0 0.4470 0.7410],'r'})
    end
else
    plot(ax1,n,vals,'.','MarkerSize',15)
    plot(ax2,n,vals,'.','MarkerSize',15)
end
    h1 = plot(ax2,x1,fit1,'LineWidth',2);
    h2 = plot(ax2,x2,fit2,'LineWidth',2);
    patch([transpose(x1), fliplr(transpose(x1))],[transpose(clow1),fliplr(transpose(chigh1))],h1.Color,'FaceAlpha',0.5,'EdgeColor','none');
    patch([transpose(x2), fliplr(transpose(x2))],[transpose(clow2),fliplr(transpose(chigh2))],h2.Color,'FaceAlpha',0.5,'EdgeColor','none');
    h = plot(ax1,n,Pfit,':','LineWidth',2);
    figure(fig1);
    patch([transpose(n), fliplr(transpose(n))],[transpose(clow),fliplr(transpose(chigh))],h.Color,'FaceAlpha',0.5,'EdgeColor','none');

end