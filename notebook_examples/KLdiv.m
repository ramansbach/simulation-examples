function k = KLdiv(P,Q)
% calculates the Kulback-Liebler divergence of Q from P Dkl(P||Q) (which is the
% measurement of how much information is lost when Q is used to approximate
% P).  P and Q are (not necessarily normalized) histogram vectors of the
% same length with counts in each entry
% not defined for Q = 0, P!= 0, need to remove these (since these are a
% binning problem more than anything)
NQ = sum(Q);
lambda = 1;
Q = Q + lambda/(NQ+1); %regularize
P = P./sum(P); %normalize
Q = Q./sum(Q);
%defMat = ~((P~=0).*(Q==0)); %indices for which KL div is defined
%P = P(defMat);
%Q = Q(defMat);

d = P./Q;
d(isnan(d))=0; %deal with 0/0
ld = log(d);
ld(ld==-Inf)=-realmax;
ld(ld==Inf)=realmax;
k = sum(P.*ld);