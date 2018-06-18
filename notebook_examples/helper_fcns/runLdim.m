 dataMatrix = dataCellFinal{a}{sc}{sig};
            if ~isempty(dataMatrix)
                x = log(dataMatrix(:,1));
                
                
      
                imin = find(abs(x-x1)==min(abs(x-x1)));
                imax = find(abs(x-x2)==min(abs(x-x2)));
                d1s = zeros([1 runs]);
                d2s = zeros([1 runs]);
                splts = zeros([1 runs]);
                for run = 2:(runs+1)
                    y = log(dataMatrix(:,run));
                    [l1,l2,splt,lsing,aic1,aic2,fig1,fig2] = Lmethod(x(imin:imax),y(imin:imax),[],0,0);
                    savefig(fig1,[saveloc '/' AAnames{a} '-' SCnames{sc} '-' Signames{sig} '-Lmeth1-' int2str(run)]);
                    saveas(fig1,[saveloc '/' AAnames{a} '-' SCnames{sc} '-' Signames{sig} '-Lmeth1-' int2str(run)],'pdf');
                    savefig(fig2,[saveloc '/' AAnames{a} '-' SCnames{sc} '-' Signames{sig} '-Lmeth2-' int2str(run)]);
                    saveas(fig2,[saveloc '/' AAnames{a} '-' SCnames{sc} '-' Signames{sig} '-Lmeth2-' int2str(run)],'pdf');
                    close(fig1);
                    close(fig2);
                    AA = str2double(AAnames{a})/100;
                    if strcmp(SCnames{sc},'02')
                       SC = 0.2; 
                    else
                       SC = str2double(SCnames{sc}); 
                    end
                    if strcmp(Rnames{sig},'001')
                        Sig = 0.01;
                    else
                        Sig = str2double(Rnames{sig});
                    end
                    if aic1 < aic2
                       d1 = lsing.a;
                       d2 = lsing.a;
                       xsplit = NaN;
                    else
                        d1 = l1.a;
                        d2 = l2.a;
                        xsplit = x(splt);
                    end
                    d1s(run-1) = d1;
                    d2s(run-1) = d2;
                    splts(run-1) = xsplit;
                end
                D1s(ind,:) = d1s;% = [D1s; d1s];
                D2s(ind,:) = d2s;% = [D2s; d2s];
                d1 = mean(d1s);
                sd1 = std(d1s);
                d2 = mean(d2s);
                sd2 = std(d2s);
                xsplit = mean(exp(splts));
                sxsplit = std(exp(splts));
                dimvec = [AA SC Sig d1 d2 xsplit sd1 sd2 sxsplit];
                dims(ind,:) = dimvec;%[dims;dimvec];
                ind = ind + 1;
            end
                