function [clmeans,olmeans,sclmeans,solmeans,cllambdas,ollambdas,cfracs,ofracs,scfracs,sofracs,clamb,sclamb,olamb,solamb] = analyzeLengthDistribution(fbase,runs,molno,cutoff,tsteps,tslices,binno,saveloc,figname,dataname,computeMu,spl1,spl2,spl3,splmid,legendson)
%This function reads in a length distribution from fbase_i where i = 1:runs
%It then computes the first, second, and infinite moments, and plots with
%stddev and saves those figures
%It also computes and plots the change in the distribution over 10
%different slices 
%%%%%%%%%%%%%
%Parameters:
%%%%%%%%%%%%%%
% fbase: string, the base name of the file where the data is to be found
% runs: int, the number of independent simulations
% molno: the total number of molecules in the system
% cutoff: cutoff for contact clusters
% saveloc: the folder to save the figures and data in
% figname: the basic name of the figures
% dataname: string, the basic name of the data
% tsteps: int, the total number of timesteps
% binno: int, number of bins for length distribution computations
% tslices: number of slices for which to compute length distributions
% computeMu: bool, whether to compute the different moments or try to load
% them
% spl1, spl2, spl3: floats, denoting the approximate location of regime
% changeovers from fractal dimension analysis
% splmid: float, denoting the approximate range of the middle
% regime changeover from fractal dimension analysis
% legendson: bool, whether to plot with legends or not

addpath('/home/rachael/Analysis_and_run_code/data_visualization')
%read in data

if exist([saveloc dataname 'distribCIDs.mat'],'file')==2
    load([saveloc dataname 'distribCIDs.mat'])
else
    cDistribs = zeros(tsteps,molno,runs);
    oDistribs = zeros(tsteps,molno,runs);
    oCIDs = zeros(tsteps,molno,runs);
    cCIDs = zeros(tsteps,molno,runs);

    for run = 1:runs
        cDistribs(:,:,run) = importdata([fbase '_runldistrib' int2str(run) '_C.dat']);
        oDistribs(:,:,run) = importdata([fbase '_runldistrib' int2str(run) '_O.dat']);
        oCIDs(:,:,run) = importdata([fbase '_runcut0.35' int2str(run) 'optical-CIDs.dat']);
        cCIDs(:,:,run) = importdata([fbase '_runcut' cutoff int2str(run) 'contact-CIDs.dat']);
    end
    cDistribs(cDistribs==0) = NaN;
    oDistribs(oDistribs==0) = NaN;

    %save variables

    save([saveloc dataname 'distribCIDs.mat'],'cDistribs','oDistribs','oCIDs','cCIDs')
end
%eliminate values > spl3 as spanning clusters
cDistribs(isnan(cDistribs)) = 0;
oDistribs(isnan(oDistribs)) = 0;
cDistribs(cDistribs>spl3) = spl3;
oDistribs(oDistribs>spl3) = spl3;
%eliminate values before time steps of 30 for equilibration
cDistribs = cDistribs(30:end,:,:);
oDistribs = oDistribs(30:end,:,:);
cCIDs = cCIDs(30:end,:,:);
oCIDs = oCIDs(30:end,:,:);
tsteps = tsteps - 30;
%compute and plot length distribution slices
figure()
hold on;
set(gca,'fontsize',20)
ldata = [];
lnames = {};
lind = 1;
mx = 0;
microconvert = (2*10^7/tsteps)*0.001*33/1000;

for t = [1:round(tsteps/tslices):tsteps tsteps]
    try
   [mids,rcounts] = computeDistributionT(t,runs,binno,cDistribs,cCIDs,Inf);
    catch
        2
    end
   if t==tsteps
   terror = errorbar(mids,mean(rcounts),std(rcounts)/sqrt(5),'.');
   
        p = plot(mids,mean(rcounts),'linewidth',2,'color',terror.Color);
   end
   %ldata = [ldata p];
   lnames{lind} = ['t = ' num2str(t*(2*10^7/tsteps)*0.001*33/1000) ' \mu s'];
   lind = lind + 1;
   mx = max(mx,max(mean(rcounts)+std(rcounts)));
end
%compute the fraction of sizes that falls into each of the regimes denoted
%by the fractal dimension analysis
%frac = [frac1 frac2a frac2b frac3]
splmida = mean(splmid);
mrcounts = mean(rcounts);
srcounts = std(rcounts);
cfrac1 = sum(mrcounts(mids<spl1));
scfrac1 = sqrt(sum(srcounts(mids<spl1).^2));
cfrac2a = sum(mrcounts(mids >spl1 & mids < splmida));
scfrac2a = sqrt(sum(srcounts(mids>spl1 & mids < splmida).^2));
cfrac2b = sum(mrcounts(mids >splmida & mids < spl2));
scfrac2b = sqrt(sum(srcounts(mids>splmida & mids < spl2).^2));
cfrac3 = sum(mrcounts(mids >spl2 & mids < spl3));
scfrac3 = sqrt(sum(srcounts(mids>spl2 & mids < spl3).^2));
cfracs = [cfrac1 cfrac2a cfrac2b cfrac3];
scfracs = [scfrac1 scfrac2a scfrac2b scfrac3];

if legendson
    legend(ldata,lnames);
end
%mx = 0.5;
ylim([0 mx]);
plot([spl1 spl1],[0 mx],'k--','linewidth',2)
plot([spl2 spl2],[0 mx],'k--','linewidth',2)
plot([spl3 spl3],[0 mx],'k--','linewidth',2)
plot([mean(splmid) mean(splmid)],[0 mx],'r--','linewidth',2)
%patch([splmid(1) splmid(2) splmid(2) splmid(1)],[0 0 mx mx],'r','edgecolor','None','facealpha',0.3);
savefig([saveloc figname 'cDistribs'])
saveas(gcf,[saveloc figname,'cDistribs'],'pdf')
close()
efit = fittype([num2str(spl3) '*(1-exp(-lambda*x))']);
%efit = fittype('lambda*log(x+1)');

figure()
hold on;
set(gca,'fontsize',20)
ldata = [];
lnames = {};
lind = 1;
mx = 0;
for t = [1:round(tsteps/tslices):tsteps tsteps]
   [mids,rcounts] = computeDistributionT(t,runs,binno,oDistribs,oCIDs,Inf); 
   if t==tsteps
   terror = errorbar(mids,mean(rcounts),std(rcounts)/sqrt(5),'.');
  
   p = plot(mids,mean(rcounts),'linewidth',2,'color',terror.Color);
   end
   %ldata = [ldata p];
   lnames{lind} = ['t = ' num2str(t*(2*10^7/tsteps)*0.001*33/1000) ' \mu s'];
   lind = lind + 1;
   mx = max(mx,max(mean(rcounts)+std(rcounts)));
end

mrcounts = mean(rcounts);
srcounts = std(rcounts);
ofrac1 = sum(mrcounts(mids<spl1));
sofrac1 = sqrt(sum(srcounts(mids<spl1).^2));
ofrac2a = sum(mrcounts(mids >spl1 & mids < splmida));
sofrac2a = sqrt(sum(srcounts(mids>spl1 & mids < splmida).^2));
ofrac2b = sum(mrcounts(mids >splmida & mids < spl2));
sofrac2b = sqrt(sum(srcounts(mids>splmida & mids < spl2).^2));
ofrac3 = sum(mrcounts(mids >spl2 & mids < spl3));
sofrac3 = sqrt(sum(srcounts(mids>spl2 & mids < spl3).^2));
ofracs = [ofrac1 ofrac2a ofrac2b ofrac3];
sofracs = [sofrac1 sofrac2a sofrac2b sofrac3];

if legendson
    legend(ldata,lnames);
end
%mx = 0.5;
ylim([0 mx]);
plot([spl1 spl1],[0 mx],'k--','linewidth',2)
plot([spl2 spl2],[0 mx],'k--','linewidth',2)
plot([mean(splmid) mean(splmid)],[0 mx],'r--','linewidth',2)
plot([spl3 spl3],[0 mx],'k--','linewidth',2)
%patch([splmid(1) splmid(2) splmid(2) splmid(1)],[0 0 mx mx],'r','edgecolor','None','facealpha',0.3);
savefig([saveloc figname 'oDistribs'])
saveas(gcf,[saveloc figname,'oDistribs'],'pdf')
close()
%compute moments
if exist([saveloc dataname 'lmeans.mat'],'file')==2 && ~computeMu
    load([saveloc dataname 'lmeans.mat']);
    cl1mean = clmeans(1,:);
    cl2mean = clmeans(2,:);
    clinfmean = clmeans(3,:);
    ol1mean = olmeans(1,:);
    ol2mean = olmeans(2,:);
    olinfmean = olmeans(3,:);
    scl1mean = sclmeans(1,:);
    scl2mean = sclmeans(2,:);
    sclinfmean = sclmeans(3,:);
    sol1mean = solmeans(1,:);
    sol2mean = solmeans(2,:);
    solinfmean = solmeans(3,:);
else
    cl1s = zeros(runs,tsteps+1);
    cl2s = zeros(runs,tsteps+1);
    clinfs = zeros(runs,tsteps+1);
    ol1s = zeros(runs,tsteps+1);
    ol2s = zeros(runs,tsteps+1);
    olinfs = zeros(runs,tsteps+1);
    
    clfits = {};
    olfits = {};

    cllambdas = [];
    ollambdas = [];
    for run = 1:runs
        [cl1,cl2,clinf] = lMoments(cDistribs(:,:,run),cCIDs(:,:,run));
        [ol1,ol2,olinf] = lMoments(oDistribs(:,:,run),oCIDs(:,:,run));
        cl1s(run,:) = cl1;
        cl2s(run,:) = cl2;
        clinfs(run,:) = clinf;
        ol1s(run,:) = ol1;
        ol2s(run,:) = ol2;
        olinfs(run,:) = olinf;
    end
    for run = 1:runs
        cl2n0 = cl2s(run,:);
        ol2n0 = ol2s(run,:);
        cshift = mean(cl2s(:,1));
        oshift = mean(ol2s(:,1));
        cl2n0(isnan(cl2n0))=0;
        ol2n0(isnan(ol2n0))=0;
        cfit = fit(microconvert*(0:(tsteps))',cl2n0'-cshift,efit,'Lower',0);
        ofit = fit(microconvert*(0:(tsteps))',ol2n0'-oshift,efit,'Lower',0);
        clfits{run} = cfit;
        olfits{run} = ofit;
        cllambdas(run) = cfit.lambda;
        ollambdas(run) = ofit.lambda;
       
    end
    cl1mean = mean(cl1s);
    cl2mean = mean(cl2s);
    clinfmean = mean(clinfs);
    ol1mean = mean(ol1s);
    ol2mean = mean(ol2s);
    olinfmean = mean(olinfs);
    scl1mean = std(cl1s);
    scl2mean = std(cl2s);
    sclinfmean = std(clinfs);
    sol1mean = std(ol1s);
    sol2mean = std(ol2s);
    solinfmean = std(olinfs);

    clmeans = [cl1mean;cl2mean;clinfmean];
    olmeans = [ol1mean;ol2mean;olinfmean];
    sclmeans = [scl1mean;scl2mean;sclinfmean];
    solmeans = [sol1mean;sol2mean;solinfmean];

save([saveloc dataname 'lmeans.mat'],'clfits','olfits','cllambdas','ollambdas','cl1s','cl2s','clinfs','ol1s','ol2s','olinfs','clmeans','olmeans','sclmeans','solmeans');
end
%plot moments

cols = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840],[0 0 0]};
lineProps = struct;
lineProps.col = cols;
lineProps.style = {'-','-','-','-','-','-','-','-'};

figure()
hold on
set(gca,'fontsize',20)
mseb(1:(tsteps+1),[cl1mean; ol1mean],[scl1mean;sol1mean],lineProps);
if legendson
    legend('contact clusters','optical clusters')
end
plot([0,tsteps],[spl1 spl1],'k--','linewidth',2)
plot([0,tsteps],[spl2 spl2],'k--','linewidth',2)
plot([0,tsteps],[spl3 spl3],'k--','linewidth',2)
%plot([0 tsteps],[splmid splmid],'r--','linewidth',2)
cshift = mean(cl2s(:,1));
oshift = mean(ol2s(:,1));
patch([0 tsteps tsteps 0],[splmid(1) splmid(1) splmid(2) splmid(2)],'m','edgecolor','None','facealpha',0.3);
ylim([0 max(cl1mean+scl1mean+5)]);
savefig([saveloc figname 'len1'])
saveas(gcf,[saveloc figname 'len1'],'pdf')
close()

figure()
hold on
set(gca,'fontsize',20)
x = microconvert*(0:(tsteps));
mseb(x,[cl2mean; ol2mean],[scl2mean;sol2mean],lineProps);
if legendson
    legend('contact clusters','optical clusters')
end

clamb = mean(cllambdas);
sclamb = std(cllambdas);
olamb = mean(ollambdas);
solamb = std(ollambdas);
plot(x,spl3*(1-exp(-clamb*x))+cshift,'b--','linewidth',2)
plot(x,spl3*(1-exp(-olamb*x))+oshift,'r--','linewidth',2)
%plot([0,tsteps],[spl1 spl1],'k--','linewidth',2)
%plot([0,tsteps],[spl2 spl2],'k--','linewidth',2)
%plot([0,tsteps],[spl3 spl3],'k--','linewidth',2)
%plot([0 tsteps],[splmid splmid],'r--','linewidth',2)
%patch([0 tsteps tsteps 0],[splmid(1) splmid(1) splmid(2) splmid(2)],'m','edgecolor','None','facealpha',0.3);
ylim([0 max(cl2mean+scl2mean+5)]);
%ylim([0 105]);
savefig([saveloc figname 'len2'])
saveas(gcf,[saveloc figname 'len2'],'pdf')
close()

figure()
hold on
set(gca,'fontsize',20)
mseb(1:(tsteps+1),[clinfmean; olinfmean],[sclinfmean;solinfmean],lineProps);
if legendson
    legend('contact clusters','optical clusters')
end
plot([0,tsteps],[spl1 spl1],'k--','linewidth',2)
plot([0,tsteps],[spl2 spl2],'k--','linewidth',2)
plot([0,tsteps],[spl3 spl3],'k--','linewidth',2)
%plot([0 tsteps],[splmid splmid],'r--','linewidth',2)
patch([0 tsteps tsteps 0],[splmid(1) splmid(1) splmid(2) splmid(2)],'m','edgecolor','None','facealpha',0.3);
ylim([0 max(spl3+5,max(clinfmean+sclinfmean+5))]);
savefig([saveloc figname 'leninf'])
saveas(gcf,[saveloc figname 'leninf'],'pdf')
close()

clear('cDistribs','oDistribs','oCIDs','cCIDs');
end

function [mids,rcounts] = computeDistributionT(i,runs,binno,distribs,cids,maxlen)
    %compute length distribution slices
    mx = 0;
    cldis = {};
    for run = 1:runs
        
        cldi = clDistrib(distribs(i,:,run),cids(i,:,run),maxlen);

        mx = max(mx,max(cldi));
        cldis{run} = cldi;
    end
    edges = 0:(mx/binno):mx;
    mids = (0.5*mx/binno):(mx/binno):mx;
    rcounts = zeros(runs,binno);
    if ~isempty(edges)
        for run = 1:runs
           rcounts(run,:) = histcounts(cldis{run},edges,'normalization','probability'); 
        end
    
    else
        mids = 1:size(rcounts,2);
    end
end
