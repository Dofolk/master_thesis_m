% D is the date string in kyma
% d is the date num from D
% V is the RPM value in kyma with status sign in columnn 2


load('kymadata_withsign.mat')
statustmp = [];
status = {};

findneg = find(V<0);
for i = 1:length(findneg)
    if( V(findneg(i)-1)==0 | V(findneg(i)+1)==0 )
        V(findneg(i)) = 0;
    end
end

fivemin = 0.0035;
d = datenum(D);
v1 = V(:,1);
clear T D

for i = 1:length(tt)
    logi = ((d>=tt(i))&(d<=tt(i)+fivemin));
    tmp = mean(v1(logi));
    statustmp(i) = tmp;
end

nantmp = (isnan(statustmp)|statustmp<0);
F(nantmp,:) = [];
tt(nantmp) = [];
statustmp(nantmp) = [];

for i = 1:length(statustmp)
    v = statustmp(i);
    if(v>=58)
        status{i} = "high";
    elseif(v>0)
        status{i} = "mid";
    elseif(v==0)
        status{i} = "stop";
    end
end
status = status';
status = string(status);
clear tmp statustmp v logi findneg