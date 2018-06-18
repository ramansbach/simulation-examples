function [contactMatrix,opticalMatrix,alignedMatrix,ctypemat] = readClusteringData(fileLocation,fbase)
bbs = 3;
aas = 5;
scs = 5;
ctypesC = [0;0;0;1;2;0;0;2;2;2;2;2;3;3;3;1;1;3;3;3;1;4;3;3;3;0;0;1;2;2;0;0;4;3;3;2;1;2;3;3;1;2;2;3;3;1;1;2;3;3;1;1;1;1;1];
ctypemat = zeros(bbs*aas*scs,6);
contactMatrix = nan*zeros(length(ctypemat),9);
opticalMatrix = nan*zeros(length(ctypemat),9);
alignedMatrix = nan*zeros(length(ctypemat),9);
ctypes = ones(bbs*aas*scs,1);
ctypes(1:length(ctypesC)) = ctypesC;
i = 2;
for bb = [0, 0.01,0.6]
    for aa = [0, 0.37,1.85,3.7,7.4]
        for sc = [0,0.2,2,4,10]
            AA = num2str(100*aa);
            if sc == 0.2
               SC = '02'; 
            else
               SC = num2str(sc); 
            end
            if bb == 0.01
                BB = '001';
            elseif bb == 0.6
                BB = '06';
            else
                BB = num2str(bb);
            end
            basename = [AA '-' SC '-' BB];
            foldname = [fileLocation basename '/'];
            fnameC = [foldname fbase basename '_runsmol-data-contact.dat'];
            if exist(fnameC,'file') == 2
                fnameC = [foldname fbase basename '_runsmol-data-contact.dat'];
                fnameO = [foldname fbase basename '_runsmol-data-optical.dat'];
                fnameA = [foldname fbase basename '_runsmol-data-aligned.dat'];
                dataC = importdata(fnameC,' ',1);
                dataO = importdata(fnameO,' ',1);
                dataA = importdata(fnameA,' ',1);
                ctypemat(i,:) = [aa sc bb ctypes(i) 0.5 0.027];
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
% BB = '0';
% basename = [AA '-' SC '-' BB];
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
% BB = '4';
% basename = [AA '-' SC '-' BB];
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
% BB = '1';
% basename = [AA '-' SC '-' BB];
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
% BB = '4';
% basename = [AA '-' SC '-' BB];
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
% BB = '0';
% basename = [AA '-' SC '-' BB];
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
