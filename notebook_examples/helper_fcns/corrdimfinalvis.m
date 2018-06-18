f = figure()
set(gca,'fontsize',20)
hold on;
runs = 5;
colormatrix = zeros(48,3);
colormatrix(:,1) = ones(48,1);
colormatrix(:,2) = 0:0.8/48:0.79;
colormatrix(:,3) = 0:0.8/48:0.79;
ind = 1;
for a = 1:length(AAnames)
    for sc = 1:length(SCnames)
        for sig = 1:length(Rnames)
            dataMatrix = dataCellFinal{a}{sc}{sig};
            if ~isempty(dataMatrix)
                x = dataMatrix(:,1);
                %y = mean(dataMatrix(:,2:end),2);
                %yMeans(:,yind) = y;
                
                %yind = yind + 1;
                %s = std(dataMatrix(:,2:end),0,2);
                %errorbar(x,y,s,'linewidth',2)
                for run = 1:runs
                    y = dataMatrix(:,1+run);
                    plot(log(x),log(y),'-','linewidth',1.5,'color',colormatrix(ind,:))
                end
                figure(f)
                ind = ind + 1;
                %xlim([log(1.5) log(73.3)])
                %xlim([0 4.5])
                
                %ylim([-15 0])
            end
        end
    end
end