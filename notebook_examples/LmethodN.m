function [lNs,gPmins,splts,lsing,gP,aic1,aicN,n,vals] = LmethodN(n,vals,weights,N,weightbool)
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
%uns = unique(n);
[n,ninds] = sort(n);
vals = vals(ninds);
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

possplts = nchoosek(2:(l-1),N-1);
possplts = [ones([size(possplts,1),1]) possplts l*ones([size(possplts,1),1])];
%possplts = possplts(2:end-1,:);
% fid = fopen('poss6.dat','w');
 %size(possplts,1)*size(possplts,2)
 %fprintf(fid,'%d \n',transpose(reshape(transpose(possplts),1,size(possplts,1)*size(possplts,2))));
 %fclose(fid);
for p = 1:size(possplts,1)
    [num2str(p) ' out of ' num2str(size(possplts,1))]
    Ps = cell(size(possplts,2)-1,1);
    gPs = cell(size(possplts,2)-1,1);
    if weightbool
        %perform weighted fit
        for s = 1:(size(possplts,2)-1)
            j = possplts(p,s);
            jn = possplts(p,s+1);
            [P,gP] = fit(n(j:jn),vals(j:jn),lf,'Weights',1./(weights(1:j).^2));
            Ps{s} = P;
            gPs{s} = gP;
        end
    else
        for s = 1:(size(possplts,2)-1)
            j = possplts(p,s);
            jn = possplts(p,s+1);
            try
            [P,gP] = fit(n(j:jn),vals(j:jn),lf);
            catch
               2; 
            end
            Ps{s} = P;
            gPs{s} = gP;
        end
    end
    RMSET = 0;
    for s = 1:(size(possplts,2)-1)
        j = possplts(p,s);
        jn = possplts(p,s+1);
        RMSE = gPs{s}.rmse;
        RMSET = RMSET + ((jn-j+1)/(length(vals)))*RMSE;
    end
    
    if RMSET < minerr
        minerr = RMSET;
        lNs = Ps;
        gPmins = gPs;
        splts = possplts(p,:);
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
rsesing = sum(rse.^2);

rseN = 0;
for s = 1:(length(splts)-1)
   j = splts(s);
   jn = splts(s+1);
   Pf =  lNs{s}.a*n(j:jn)+lNs{s}.b;
   res = Pf - vals(j:jn);
   rseN = rseN+sum(res.^2);
end

aic1 = akaike(3,l,rsesing);
aicN = akaike(3+2*(length(splts)-2),l,rseN);

