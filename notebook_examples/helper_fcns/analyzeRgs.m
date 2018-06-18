%script that just grabs the rgs versus mass for all of the different runs
%and rgs versus t and plots (a) asphericity, (b) ratio of r3/r1 and r2/r1
%as well as rinf maybe at each time step?
function [um1,um2,um3,umg,umk,umr1,umr2,umr3,mkap2,mekap2,mr1,mr2,mr3,mRg,mer1,mer2,mer3,meRg,r1,r2,r3,er1,er2,er3] = analyzeRgs(fmids,outfname,sterr,runfigs,throwout)
%fbase = '/home/rachael/cluster/data/mansbac2/coarsegraining/MD/DFAG/';
fbase = '';
%fmids = {'42mols_run1','42mols_run2','42mols_run3'};%,'126mols_run1','126mols_run2','126mols_run3','126mols_run4','200mols_run1','200mols_run2','200mols_run3','200mols_runIC1','200mols_runIC2','200mols_runIC3','252mols_run1','252mols_run2','252mols_run3','252mols_run4','252mols_run5','378mols_run2','378mols_run2','378mols_run3','378mols_run4'};
fendm = '/clustdata/mass_v_Rg.dat';
fend = '/clustdata/clustScalars.dat';
fout = '/home/rachael/cluster/home/mansbac2/coarsegraining/code/hoomd-helpers/analysis-output/';
maxRg = Inf;%approximate maximum "trustworthy" Rg in a box of size 23.5
%preallocate array
rtotm = 0;
rtoth = 0;
%ctot = 0;
if ~exist('sterr','var')
   sterr = 0; 
end
if ~exist('runfigs','var')
   runfigs = 0; 
end
if ~exist('throwout','var')
   throwout = 5; 
end
for f = 1:length(fmids)
   [rm,ctotm] = matsize([fbase fmids{f} fend],'%f');
   [r,ctot] = matsize([fbase fmids{f} fendm],'%f');
   rtotm = rtotm+rm;
   rtoth = rtoth + r;
   %ctot = ctot+c;
end
Rgs = zeros(rtotm,ctotm);
Rhs = zeros(rtoth,ctot);
ind1 = 1;
ind2 = 1;
%read in data
for f = 1:length(fmids)
    [mat,rm,~] = readlinebyline([fbase fmids{f} fend],'%f',[]);
    [mat2,r,~] = readlinebyline([fbase fmids{f} fendm],'%f',[]);
    
    
    try
    Rgs(ind1:ind1+rm-1,:) = mat;
    Rhs(ind2:ind2+r-1,:) = mat2;
    catch
           size(mat)
    end
    ind1 = ind1+rm;
    ind2 = ind2+r;
    
end
%only consider Rgs less than the max Rg
%average Rgs at each time step
%can do them as standard error instead of standard deviation if desired
% [r1,er1,ut] = avg_by_var(Rgs(:,1),Rgs(:,2),5,1,@mean,Inf,sterr);
% [r2,er2,~] = avg_by_var(Rgs(:,1),Rgs(:,3),5,1,@mean,Inf,sterr);
% [r3,er3,~] = avg_by_var(Rgs(:,1),Rgs(:,4),5,1,@mean,Inf,sterr);
% [as,eas,~] = avg_by_var(Rgs(:,1),Rgs(:,6),5,1,@mean,Inf,sterr);
% [Rg,eRg,~] = avg_by_var(Rgs(:,1),Rgs(:,5),5,1,@mean,Inf,sterr);
% ut = ut(Rg<=maxRg);
% 
% eRg = eRg(Rg<=maxRg);
% as = as(Rg<=maxRg);
% eas = eas(Rg<=maxRg);
% r3 = r3(Rg<=maxRg);
% er3 = er3(Rg<=maxRg);
% r2 = r2(Rg<=maxRg);
% er2 = er2(Rg<=maxRg);
% r1 = r1(Rg<=maxRg);
% er1 = er1(Rg<=maxRg);
% Rg = Rg(Rg<=maxRg);

%average Rgs at each cluster size
kap2 = 3/2*((Rgs(:,2).^2+Rgs(:,3).^2+Rgs(:,4).^2)./(Rgs(:,2)+Rgs(:,3)+Rgs(:,4)).^2)-1/2;
[mr1,mer1,um1] = avg_by_var(Rhs(:,1),sqrt(Rgs(:,2)),throwout,1,@mean,Inf,sterr);
[mr2,mer2,um2] = avg_by_var(Rhs(:,1),sqrt(Rgs(:,3)),throwout,1,@mean,Inf,sterr);
[mr3,mer3,um3] = avg_by_var(Rhs(:,1),sqrt(Rgs(:,4)),throwout,1,@mean,Inf,sterr);
[mas,meas,uma] = avg_by_var(Rhs(:,1),Rgs(:,6),throwout,1,@mean,Inf,sterr);
[mRg,meRg,umg] = avg_by_var(Rhs(:,1),sqrt(Rgs(:,5)),throwout,1,@mean,Inf,sterr);
[mkap2,mekap2,umk] = avg_by_var(Rhs(:,1),kap2,throwout,1,@mean,Inf,sterr);
meRg = meRg(mRg<=maxRg);
mas = mas(mRg<=maxRg);
meas = meas(mRg<=maxRg);
mr3 = mr3(mRg<=maxRg);
mer3 = mer3(mRg<=maxRg);
mr2 = mr2(mRg<=maxRg);
mer2 = mer2(mRg<=maxRg);
mr1 = mr1(mRg<=maxRg);
mer1 = mer1(mRg<=maxRg);
mkap2 = mkap2(mRg<=maxRg);
mekap2 = mekap2(mRg<=maxRg);
um1 = um1(mRg<=maxRg);
um2 = um2(mRg<=maxRg);
um3 = um3(mRg<=maxRg);
uma = uma(mRg<=maxRg);
umk = umk(mRg<=maxRg);
umg = umg(mRg<=maxRg);
mRg = mRg(mRg<=maxRg);

r1 = Rgs(:,4)./Rgs(:,2);
r2 = Rgs(:,4)./Rgs(:,3);
r3 = Rgs(:,3)./Rgs(:,2);
[r1,er1,umr1] = avg_by_var(Rhs(:,1),sqrt(Rgs(:,4))./sqrt(Rgs(:,2)),throwout,1,@mean,Inf,sterr);
[r2,er2,umr2] = avg_by_var(Rhs(:,1),sqrt(Rgs(:,4))./sqrt(Rgs(:,3)),throwout,1,@mean,Inf,sterr);
[r3,er3,umr3] = avg_by_var(Rhs(:,1),sqrt(Rgs(:,3))./sqrt(Rgs(:,2)),throwout,1,@mean,Inf,sterr);
%r1 = mr3./mr1;
%r2 = mr3./mr2;
%r3 = mr2./mr1;
%er1 = sqrt(mer1.^2+mer3.^2);
%er2 = sqrt(mer2.^2+mer3.^2);
%er3 = sqrt(mer2.^2+mer1.^2);
%ut = ut/1000;

% figure()
% errorbar(ut,as,eas,'.')
% ylim([-0.5 3])
% xlabel('t (ns)')
% title('asphericity')
% set(gca,'FontSize',15)
% 
% savefig([fout num2str(nmols) 'mols_aspher.fig'])
% saveas(gcf,[fout num2str(nmols) 'mols_aspher'],'epsc')
% saveas(gcf,[fout num2str(nmols) 'mols_aspher.jpg'])

% figure()
% hold on
% errorbar(ut,r3,er3,'.')
% errorbar(ut,r2,er2,'.')
% errorbar(ut,r1,er1,'.')
% title('moments of gyration tensor')
% xlabel('t (ns)')
% legend('r_1','r_2','r_3')
% ylim([0 3.5])
% set(gca,'FontSize',15)
% 
% savefig([fout num2str(nmols) 'mols_moments.fig'])
% saveas(gcf,[fout num2str(nmols) 'mols_moments'],'epsc')
% saveas(gcf,[fout num2str(nmols) 'mols_moments.jpg'])
% 
% figure()
% hold on
% plot(ut,r1,'-^')
% plot(ut,r2,'-o')
% plot(ut,r3,'-s')
% plot(ut,Rg,'-x')
% xlabel('t (ns)')
% ylim([0 3])
% title('moments of gyration tensor')
% legend('r_3','r_2','r_1','R_g')
% set(gca,'FontSize',15)
% 
% savefig([fout num2str(nmols) 'mols_momentsRg.fig'])
% saveas(gcf,[fout num2str(nmols) 'mols_momentsRg'],'epsc')
% saveas(gcf,[fout num2str(nmols) 'mols_momentsRg.jpg'])
% 
% figure()
% hold on
% errorbar(ut,r3./r1,sqrt(er1.^2+er3.^2),'^','linewidth',2,'markersize',20)
% errorbar(ut,r3./r2,sqrt(er2.^2+er3.^2),'o','linewidth',2,'markersize',20)
% errorbar(ut,r2./r1,sqrt(er1.^2+er2.^2),'s','linewidth',2,'markersize',20)
% xlabel('t (ns)')
% legend('r_1/r_3','r_1/r_2','r_2/r_3')
% ylim([0 30])
% set(gca,'FontSize',24)
% 
% savefig([fout num2str(nmols) 'mols_gyrrats.fig'])
% saveas(gcf,[fout num2str(nmols) 'mols_gyrrats'],'epsc')
% saveas(gcf,[fout num2str(nmols) 'mols_gyrrats.jpg'])
% 
if runfigs
figure()
errorbar(uma,mas,meas,'-x','linewidth',2)
xlabel('mass')
%title('asphericity')
set(gca,'FontSize',24)

savefig([fout outfname 'mols_aspherM.fig'])
saveas(gcf,[fout outfname 'mols_aspherM'],'epsc')
saveas(gcf,[fout outfname 'mols_aspherM.jpg'])

figure()
errorbar(umk,mkap2,mekap2,'-x','linewidth',2)
set(gca,'FontSize',24)
xlim([-10 340])
ylim([0 1.1])
savefig([fout outfname 'mols_anisotM.fig'])
saveas(gcf,[fout outfname,'mols_anisotM'],'epsc')
saveas(gcf,[fout outfname,'mols_anisotM.jpg'])


figure()
mseb(transpose(umk),transpose(mkap2),transpose(mekap2));
set(gca,'FontSize',24)
xlim([-10 340])
ylim([0 1.1])

savefig([fout outfname 'mols_anisotMt.fig'])
saveas(gcf,[fout outfname,'mols_anisotMt'],'epsc')
saveas(gcf,[fout outfname,'mols_anisotMt.jpg'])


figure()
hold on
errorbar(um1,mr1,mer1,'-','linewidth',2)
errorbar(um2,mr2,mer2,'--','linewidth',2)
errorbar(um3,mr3,mer3,':','linewidth',2)
errorbar(umg,mRg,meRg,'-.','linewidth',2)
%title('moments of gyration tensor')

%legend('r_3','r_2','r_1','R_g','location','NorthWest')
set(gca,'FontSize',24)
xlim([-10 340])
savefig([fout outfname 'mols_momentsRgM.fig'])
saveas(gcf,[fout outfname 'mols_momentsRgM'],'epsc')
saveas(gcf,[fout outfname 'mols_momentsRgM.jpg'])

lineProps1 = struct;
lineProps2 = struct;
lineProps3 = struct;
lineProps4 = struct;
lineProps1.col = {[0 0.4470 0.7410]};
lineProps2.col = {[0.85 0.325 0.098]};
lineProps3.col = {[0.929 0.694 0.125]};
lineProps4.col = {[0.494 0.184 0.556]};
lineProps1.style = '--';
lineProps2.style = ':';
lineProps3.style = '-';
lineProps4.style = '-.';

figure()
hold on

mseb(transpose(um1),transpose(mr1),transpose(mer1),lineProps1,1);
mseb(transpose(um2),transpose(mr2),transpose(mer2),lineProps2,1);
mseb(transpose(um3),transpose(mr3),transpose(mer3),lineProps3,1);
mseb(transpose(umg),transpose(mRg),transpose(meRg),lineProps4,1);
%title('moments of gyration tensor')

%legend('r_3','r_2','r_1','R_g','location','NorthWest')
set(gca,'FontSize',24)
xlim([-10 340])
savefig([fout outfname 'mols_momentsRgMt.fig'])
saveas(gcf,[fout outfname 'mols_momentsRgMt'],'epsc')
saveas(gcf,[fout outfname 'mols_momentsRgMt.jpg'])
% 
% figure()
% hold on
% plot(um1,mr1,'-','linewidth',2)
% plot(um2,mr2,'--','linewidth',2)
% plot(um3,mr3,'-.','linewidth',2)
% plot(umg,mRg,':','linewidth',2)
% %title('moments of gyration tensor')
% xlabel('mass')
% %legend('r_3','r_2','r_1','R_g','location','NorthWest')
% set(gca,'FontSize',24)
% 
% savefig([fout outfname 'mols_momentsRgplot.fig'])
% saveas(gcf,[fout outfname 'mols_momentsRgplot'],'epsc')
% saveas(gcf,[fout outfname 'mols_momentsRgplot.jpg'])

figure()
hold on
errorbar(um1,mr3./mr1,sqrt(mer1.^2+mer3.^2),'-','linewidth',2)
errorbar(um1,mr3./mr2,sqrt(mer2.^2+mer3.^2),'--','linewidth',2)
errorbar(um1,mr2./mr1,sqrt(mer1.^2+mer2.^2),'-.','linewidth',2)
xlabel('mass')
%legend('r_1/r_3','r_1/r_2','r_2/r_3','location','North')
set(gca,'FontSize',24)

savefig([fout outfname 'mols_gyrratsM.fig'])
saveas(gcf,[fout outfname 'mols_gyrratsM'],'epsc')
saveas(gcf,[fout outfname 'mols_gyrratsM.jpg'])




figure()
hold on
mseb(transpose(um1),transpose(mr3./mr1),transpose(sqrt(mer1.^2+mer3.^2)),lineProps1,1);
mseb(transpose(um1),transpose(mr3./mr2),transpose(sqrt(mer2.^2+mer3.^2)),lineProps2,1);
mseb(transpose(um1),transpose(mr2./mr1),transpose(sqrt(mer2.^2+mer1.^2)),lineProps3,1);
xlim([-10 340])
ylim([-5 25])
%legend('r_1/r_3','r_1/r_2','r_2/r_3','location','North')
set(gca,'FontSize',24)

savefig([fout outfname 'mols_gyrratsMt.fig'])
saveas(gcf,[fout outfname 'mols_gyrratsMt'],'epsc')
saveas(gcf,[fout outfname 'mols_gyrratsMt.jpg'])
end
% figure()
% hold on
% plot(um1,mr3./mr1,'-','linewidth',2)
% plot(um1,mr3./mr2,'--','linewidth',2)
% plot(um1,mr2./mr1,':','linewidth',2)
% xlabel('mass')
% %legend('r_1/r_3','r_1/r_2','r_2/r_3','location','North')
% set(gca,'FontSize',24)

% savefig([fout outfname 'mols_gyrratsplot.fig'])
% saveas(gcf,[fout outfname 'mols_gyrratsplot'],'epsc')
% saveas(gcf,[fout outfname 'mols_gyrratsplot.jpg'])
%savefig([fout num2str(nmols) 'mols_Rgrats.fig'])
%saveas(gcf,[fout num2str(nmols) 'mols_Rgrats'],'epsc')
%saveas(gcf,[fout num2str(nmols) 'mols_Rgrats.jpg'])

%Rgs = Rgs(Rgs(:,5)<maxRg,:);
%masses = Rhs(:,1);
 
