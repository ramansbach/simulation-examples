function dataCell = readFractalDim(tframes,AAnames,SCnames,Rnames,fileLocation)

fbase = 'mols10000_';

dataCell = {};


for i = 1:length(AAnames)
    for sc = 1:length(SCnames)
        for sig = 1:length(Rnames)
            foldname = [fileLocation AAnames{i} '-' SCnames{sc} '-' Rnames{sig} '/'];
            fname = [foldname fbase AAnames{i} '-' SCnames{sc} '-' Rnames{sig} '_short_run-corrcalc'];
            
            data = importdata([fname num2str(tframes(1)) '.dat']);

            szdata = size(data);
            dataMatrix = zeros([szdata(1) szdata(2) length(tframes)]);
            tind = 1;
               for t = tframes
                  
                   data = importdata([fname num2str(t) '.dat']);
                   
                   dataMatrix(:,:,tind) = data;
                   tind = tind + 1;
               end
               

            dataCell{i}{sc}{sig} = dataMatrix;
        end
    end
end
