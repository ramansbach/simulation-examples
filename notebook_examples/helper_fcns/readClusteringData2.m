function [contactMatrix,opticalMatrix,alignedMatrix,ctypemat] = readClusteringData2(fileLocation,fbase)
sigs = 4;
aas = 3;
scs = 3;
%ctypesC = [0;0;0;1;2;0;0;2;2;2;2;2;3;3;3;1;1;3;3;3;1;4;3;3;3;0;0;1;2;2;0;0;4;3;3;2;1;2;3;3;1;2;2;3;3;1;1;2;3;3;1;1;1;1;1];
ctypemat = zeros(sigs*aas*scs,6);
contactMatrix = nan*zeros(length(ctypemat),9);
opticalMatrix = nan*zeros(length(ctypemat),9);
alignedMatrix = nan*zeros(length(ctypemat),9);
ctypes = ones(sigs*aas*scs,1);
%ctypes(1:length(ctypesC)) = ctypesC;
i = 1;
for sig = [0.25, 0.5,0.75,1.25]
    for aa = [2.5,5.0,7.5]
        for sc = [0.2,2,4]
            AA = num2str(100*aa);
            if sc == 0.2
               SC = '02'; 
            else
               SC = num2str(sc); 
            end
            if sig == 0.01
                SIG = '001';
            elseif sig == 0.6
                SIG = '06';
            elseif sig == 1.25;
                SIG = '125';
            else
                SIG = ['0' num2str(100*sig)];
            end
            
            basename = [AA '-' SC '-' SIG];
            foldname = [fileLocation basename '/'];
            fnameC = [foldname fbase basename '_short_runsmol-data-contact.dat'];
            if exist(fnameC,'file') == 2
                fnameC = [foldname fbase basename '_short_runsmol-data-contact.dat'];
                fnameO = [foldname fbase basename '_short_runsmol-data-optical.dat'];
                fnameA = [foldname fbase basename '_short_runsmol-data-aligned.dat'];
                dataC = importdata(fnameC,' ',1);
                dataO = importdata(fnameO,' ',1);
                dataA = importdata(fnameA,' ',1);
                try
                ctypemat(i,:) = [aa sc 1.0 ctypes(i) sig 0.0027];
                catch
                    aa
                end
            
                for j = 1:9
                   contactMatrix(i,j) = dataC.data(j);
                   opticalMatrix(i,j) = dataO.data(j);
                   alignedMatrix(i,j) = dataA.data(j);
                end
                i = i+1;
            end
    
        end
    end
end
%additional things to read in by hand if necessary
% AA = '185';
% SC = '2';
% sig = '0';
% basename = [AA '-' SC '-' sig];
% foldname = [fileLocation 'diffuse/' basename '/'];
% fnameC = [foldname fbase basename '_runsmol-data-contact.dat'];
% fnameO = [foldname fbase basename '_runsmol-data-optical.dat'];
% fnameA = [foldname fbase basename '_runsmol-data-aligned.dat'];
% dataC = importdata(fnameC,' ',1);
% dataO = importdata(fnameO,' ',1);
% dataA = importdata(fnameA,' ',1);
% 
% for j = 1:9
%     contactMatrix(i,j) = dataC.data(j);
%     opticalMatrix(i,j) = dataO.data(j);
%     alignedMatrix(i,j) = dataA.data(j);
% end
% i = i+1;
% 
% AA = '0';
% SC = '0';
% sig = '4';
% basename = [AA '-' SC '-' sig];
% foldname = [fileLocation basename '/'];
% fnameC = [foldname fbase basename '_runsmol-data-contact.dat'];
% fnameO = [foldname fbase basename '_runsmol-data-optical.dat'];
% fnameA = [foldname fbase basename '_runsmol-data-aligned.dat'];
% dataC = importdata(fnameC,' ',1);
% dataO = importdata(fnameO,' ',1);
% dataA = importdata(fnameA,' ',1);
% for j = 1:9
%     contactMatrix(i,j) = dataC.data(j);
%     opticalMatrix(i,j) = dataO.data(j);
%     alignedMatrix(i,j) = dataA.data(j);
% end
% i = i+1;
% 
% AA = '185';
% SC = '2';
% sig = '1';
% basename = [AA '-' SC '-' sig];
% foldname = [fileLocation basename '/'];
% fnameC = [foldname fbase basename '_runsmol-data-contact.dat'];
% fnameO = [foldname fbase basename '_runsmol-data-optical.dat'];
% fnameA = [foldname fbase basename '_runsmol-data-aligned.dat'];
% dataC = importdata(fnameC,' ',1);
% dataO = importdata(fnameO,' ',1);
% dataA = importdata(fnameA,' ',1);
% for j = 1:9
%     contactMatrix(i,j) = dataC.data(j);
%     opticalMatrix(i,j) = dataO.data(j);
%     alignedMatrix(i,j) = dataA.data(j);
% end
% i = i+1;
% 
% AA = '185';
% SC = '2';
% sig = '4';
% basename = [AA '-' SC '-' sig];
% foldname = [fileLocation basename '/'];
% fnameC = [foldname fbase basename '_runsmol-data-contact.dat'];
% fnameO = [foldname fbase basename '_runsmol-data-optical.dat'];
% fnameA = [foldname fbase basename '_runsmol-data-aligned.dat'];
% dataC = importdata(fnameC,' ',1);
% dataO = importdata(fnameO,' ',1);
% dataA = importdata(fnameA,' ',1);
% for j = 1:9
%     contactMatrix(i,j) = dataC.data(j);
%     opticalMatrix(i,j) = dataO.data(j);
%     alignedMatrix(i,j) = dataA.data(j);
% end
% i = i+1;
% 
% AA = '185';
% SC = '2';
% sig = '0';
% basename = [AA '-' SC '-' sig];
% foldname = [fileLocation 'rad025/' basename '/'];
% fnameC = [foldname fbase basename '_runsmol-data-contact.dat'];
% fnameO = [foldname fbase basename '_runsmol-data-optical.dat'];
% fnameA = [foldname fbase basename '_runsmol-data-aligned.dat'];
% dataC = importdata(fnameC,' ',1);
% dataO = importdata(fnameO,' ',1);
% dataA = importdata(fnameA,' ',1);
% for j = 1:9
%     contactMatrix(i,j) = dataC.data(j);
%     opticalMatrix(i,j) = dataO.data(j);
%     alignedMatrix(i,j) = dataA.data(j);
% end
% i = i+1;
