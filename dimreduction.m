% F is the processed data
% t is also
clear

load('feature.mat');
feature(:,end) = [];
freq = 1:1:600;
freq = [freq 1000000];

Z1 = normalize(feature(:,1:601));
[id,score,u] = OD_onlinePCA(Z1,0.2);
outs = sum(score>(mean(score)+3*var(score)));
F = feature;
F(id(1:outs),:)=[];
tt(id(1:outs)) = [];

Z = F(:,1:600)-mean(F(:,1:600));
[C,S] = pca(Z);
outf = C(:,1)<0.01;
freq(outf) = [];
F(:,outf) = [];

clearvars -except F feature freq tt nantmp fivemin
