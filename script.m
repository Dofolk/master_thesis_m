%%
% load data file and save each info
Fi = dir([pwd filesep '*.mat']);
intx = [];inty = [];intz = [];
for i = 1:length(Fi)
    load(Fi(i).name);
    intx = [intx,intfreqX];
    inty = [inty,intfreqY];
    intz = [intz,intfreqZ];
end
feature(1:200,:) = intx;
feature(201:400,:) = inty;
feature(401:600,:) = intz;
clearvars -except feature
disp('Part1')

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
l = length(feature);
D = zeros(l);
for i =1:l
    for j = i:l
        D(i,j) = norm(feature(:,i)-feature(:,j));
    end
end
disp('Part3')

%%
% build edge matrix
r = linspace(0,0.2685,1000);
E = zeros(r,l);
for i = 1:length(r)
    T = double((D>0) & (D<r(i)));
    T = T + T';
    E(i,:) = sum(T);
end
clearvars T
%%
for k = 1:l
    pk = zeros(1,1000);
    for i = 1:length(E)
        T = E(i).X;
        pk(i) = T(k);
    end
    figure;
    plot(r,pk);
    title(sprintf('Point %d',k));
    print(sprintf('point%d.jpg',k),'-djpeg');
    close;
end