glabels = {'reg. 1','reg. 2a', 'reg. 2b', 'reg. 3'};
gfiglabels = {'1','2a','2b','3'};
for reg = 1:length(glabels)
    for i = 1:length(AAnames)
        icurr = find((dims(:,1)==(str2double(AAnames{i})/100)));
        
        f = figure();
        hold on
        f = errorMesh(dims(icurr,2),dims(icurr,3),Cfracs(icurr,reg),Cfracs(icurr,reg+4)/sqrt(runs),[0 0 0],f,0.3);
        figure(f);
        colorbar
        caxis([min(Cfracs(:,reg)) max(Cfracs(:,reg))])
        xlim([0 10.5])
        ylim([0.98 1.8])
        zlim([min(Cfracs(:,reg)-Cfracs(:,reg+4)) max(Cfracs(:,reg)+Cfracs(:,reg+4))])
        view([-120.7,30])
        set(gca,'fontsize',20);
        xlabel('SC eps')
        ylabel('SC Sig')
        zlabel(['C frac, ' glabels{reg}])
        title(['AA = ' AAnames{i}]);
        savefig([saveloc '/lfracr' gfiglabels{reg} 'A' AAnames{i} '3D'])
        saveas(gcf,[saveloc '/lfracr' gfiglabels{reg} 'A' AAnames{i} '3D'],'pdf')
        
%         f = figure();
%         hold on
%         f = errorMesh(dims(icurr,2),dims(icurr,3),log(Cfracs(icurr,reg)),zeros(size(Cfracs(icurr,reg+4))),[0 0 0],f,0.3);
%         figure(f);
%         colorbar
%         lCfracs = log(Cfracs(:,reg));
%         lCfracs = lCfracs(~isinf(lCfracs));
%         caxis([min(lCfracs) max(lCfracs)])
%         xlim([0 10.5])
%         ylim([0.98 1.8])
%         zlim([min(lCfracs) max(lCfracs)])
%         view([-46.7,7.6])
%         set(gca,'fontsize',20);
%         xlabel('SC eps')
%         ylabel('SC Sig')
%         zlabel(['C frac, log(' glabels{reg} ')'])
%         title(['AA = ' AAnames{i}]);
%         savefig([saveloc '/lfracr' gfiglabels{reg} 'A' AAnames{i} '3D'])
%        saveas(gcf,[saveloc '/lfracr' gfiglabels{reg} 'A' AAnames{i} '3D'],'pdf')
    end
end