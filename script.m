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
clearvars -except intx inty intz
disp('Part1')

%%
% do kmeans for 4 clustering
rng(1);
[idx,cx] = kmeans(intx,4);
[idy,cy] = kmeans(inty,4);
[idz,cz] = kmeans(intz,4);
disp('Part2')

%%
% build the distance matrix(D)
l = length(intx);
D = zeros(l);
for i =1:l
    for j = i:l
        D(i,j) = norm(intx(:,i)-intx(:,j));
    end
end
disp('Part3')

%%
% build edge matrix
r = linspace(0,0.2685,1000);
E = {};
for i = 1:length(r)
    T = double((D>0) & (D<r(i)));
    T = T + T';
    E(i).X = sum(T);
end

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