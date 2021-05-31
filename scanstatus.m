%Let the notations:
% D be the kyma original date stamp
% d be the kyma date num
% 5 minutes in matlab to num is 0.0035, denoted by fivemin
% logi is the logical matirx that in the interval of two datelist elements
% V is the data value matrix with operation status
% S be the array that record each intervals have what status
clear

load('kymadata_withsign.mat','D','V');
load('datadatelist.mat');

d = datenum(D);
fivemin = 0.0035;
S = zeros(length(datelist),1);
v1 = V(:,1);
v2 = V(:,2);
for i = 1:length(datelist)
    if(i==length(datelist))
        timetmp = max(datelist)+fivemin;
        logi = ((d>=max(datelist))&(d<=timetmp));
        tmp = v1(logi);
        S(i) = mean(tmp);
        disp('The Final');
        continue
    end
    logi = ((d>=datelist(i))&(d<=datelist(i+1)));
    tmp = v1(logi);
    S(i) = mean(tmp);
    modtmp = round(length(datelist)/4);
    if (mod(i,modtmp)==0)
        disp('Flag')
    end
end
save('RPMs.mat','S')
