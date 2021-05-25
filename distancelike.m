%%
% load data file and save each info
Fi = dir([pwd filesep '主機*.mat']);
intx = [];inty = [];intz = [];N = [];
for i = 1:length(Fi)
    load(Fi(i).name);
    intx = [intx,intfreqX(:,T1)];
    inty = [inty,intfreqY(:,T1)];
    intz = [intz,intfreqZ(:,T1)];
    N = [N,N_original];
end
feature(1:200,:) = intx;
feature(201:400,:) = inty;
feature(401:600,:) = intz;
feature(601,:) = N;
clearvars -except feature
disp('Part1')
pause(3)
%%
% do kmeans for 4 clustering
rng(1);
[idx,cx] = kmeans(feature(1:200,:),4);
[idy,cy] = kmeans(feature(201:400,:),4);
[idz,cz] = kmeans(feature(401:600,:),4);
[id,c] = kmeans(feature,4);
disp('Part2')

%%
% build the distance matrix(D)
%{
%wasted fk= =
l = length(feature);
D = zeros(l);
for i =1:l
    for j = i:l
        D(i,j) = norm(feature(:,i)-feature(:,j));
    end
end
%}
D = pdist(feature');
D = squareform(D);
disp('Part3')

%%
% build edge matrix
r = linspace(0,0.1,1000);
r1 = linspace(0,0.2685,1000);
E = zeros(length(r),l);
for i = 1:length(r)
    T = double((D>0) & (D<r(i)));
    T = T + T';
    E(i,:) = sum(T);
end
clearvars T

%%
%knnsearch epsilon
A = feature;
datanum = size(feature,1);
D = zeros(datanum,1);
[~,dist] = knnsearch(A(2:datanum,:),A(1,:));
D(1) = dist;
for i = 2:datanum
    [~,dist] = knnsearch(A([1:i-1,i+1:datanum],:),A(i,:));
    D(i) = dist;
end
[sortD,sortDid] = sort(D,'descend');
plot( (1:1:datanum)',sortD,'r+-','Linewidth',2);

%%
%use dbscan to find clusterings
[id,c] = dbscan(feature',0.009,50);

%%
%test dbscan min. pts.
%find that can set min pt = 190~250
% 0.0079<eps<0.01 , not equal 0.01
clearvars -except feature
M = zeros(100,1);
m = zeros(100,1);
for i=1:100
    [id,c] = dbscan(feature,0.009,250-i);
    M(i) = max(id);
    m(i) = sum(id==-1);
end

%%