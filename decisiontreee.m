a = randperm(length(tt));
s = round(length(tt)*0.7);
Train = F(a(1:s),:);
TrainP = status(a(1:s));
Test = F(a(s:end),:);
TestP = status(a(s:end));
clear a s
ctree = fitctree(Train,TrainP);

%method 1
leafs = logspace(1,2,10);
rng('default')
N = numel(leafs);
err = zeros(N,1);
for n=1:N
    t = fitctree(Train,TrainP,'CrossVal','On',...
        'MinLeaf',leafs(n));
    err(n) = kfoldLoss(t);
end
plot(leafs,err);
xlabel('Min Leaf Size');
ylabel('cross-validated error');
% OptimalTree = fitctree(F,status,'minleaf',40);

%method 2
[~,~,~,bestlevel] = cvLoss(ctree,'SubTrees','All','TreeSize','min');
tree = prune(ctree,'Level',bestlevel);
