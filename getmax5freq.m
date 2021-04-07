num = [2 3 4];colle = [];
colle = [];
for i = 1:length(num)
    
    dirname = ['結果20200' int2str(num(i)) filesep];
    disp(i);
    D = dir([dirname sprintf('*-features.mat')]);
    
    for k = 1 : length(D)
            tmp=[];
            load([dirname D(k).name]);
            J = (N_original/120000)>0.667;
            [~,IX] = sort(intfreqX,'descend');
            [~,IY] = sort(intfreqY,'descend');
            [~,IZ] = sort(intfreqZ,'descend');
            IX = IX(:,J);
            IY = IY(:,J);
            IZ = IZ(:,J);
            tmp(1:5,:) = IX(1:5,:);
            tmp(6:10,:) = IY(1:5,:);
            tmp(11:15,:) = IZ(1:5,:);
            tmp(16,:) = N_original(J)';
            tmp(17,:) = EX(J)';
            tmp(18,:) = EY(J)';
            tmp(19,:) = EZ(J)';
            colle = [colle,tmp];
            clearvars tmp
    end
    %pause(20)
end