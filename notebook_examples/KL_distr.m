N = 10000;%number of distribution generated
S1 = 180001*2;%number of sampling points--this is the same as the number of frames looked at generally
S2 = 60001*2;
KLdistr = zeros(1,N);
for i=1:N
    %draw two random vectors of the same size from the same distribution
    a = randn(S1,1);
    b = randn(S2,1);
    %compute their KL divergence
    KLdistr(i) = bondedKL(a,b,101,0,'-');
end
figure;
%create a histogram of the different KL divergences returned for each
%random sample
%so then our histogram is the probability of seeing that KLD given the
%total number of distrubutions we have generated
hist(KLdistr,100);
xlabel('KL','FontSize',30);
ylabel('counts','FontSize',30);