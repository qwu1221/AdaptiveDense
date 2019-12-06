%% Replicated simulation results
% Clique structure
% tpr and tnr are compared with greedy, adaptive greedy and goldberg's alg


rng(20);
n1 = 30;n2 = 30;N=100;m=15;
mu1 = 0.8;mu2 = 0;sigma = 1;
rate1=0.9;rate0=0.05;
iter_time = 20;
result_list = zeros(iter_time,6);
result_list_ref = zeros(iter_time,6); % use lambda=1 for comparison
result_list_gold = zeros(iter_time,6);
addpath('/Users/qwu/Downloads/Don/dense')
lambda0=0.8:0.025:1.6;

A_iter = zeros(N*(N-1)/2,n1,iter_time);
B_iter = zeros(N*(N-1)/2,n2,iter_time);
p_iter = zeros(N*(N-1)/2,iter_time);
node_perm_iter = zeros(n,iter_time);
error_idx  = zeros(iter_time,1);
for iter=1:iter_time
    iter
    try
G = binornd(1,rate0,N);
G(1:m,1:m)=binornd(1,rate1,m);
for i=1:N
    G(i,i)=0;
end
G_vec = squareform(G);
A = zeros(size(G_vec,2),n1);
B = zeros(size(G_vec,2),n2);
p_vec = zeros(size(G_vec));
for i=1:n1
    A(:,i)=normrnd(mu1*G_vec+mu2,sigma,size(G_vec));
end
A_iter(:,:,iter)=A;
for i=1:n2
    B(:,i)=normrnd(mu2,sigma,size(G_vec));
end
B_iter(:,:,iter)=B;
for i=1:size(G_vec,2)
    [h,p_vec(i)]=ttest2(A(i,:),B(i,:));
end
P=squareform(p_vec);
p_iter(:,iter)=p_vec;
for i=1:N
    P(i,i)=1;
end
nlogp=-log(P);
%figure;imagesc(nlogp);colormap jet;colorbar;snapnow

% do permutation for the observed matrix of p-values
n=size(nlogp,1);
perm_matrix = squareform(1:(n*(n-1)/2));
node_perm_idx = randperm(n);
perm_matrix = perm_matrix(node_perm_idx,node_perm_idx);
perm_vec = squareform(perm_matrix);
[result ID]=sort(perm_vec);
nlogp_vec = squareform(nlogp);
nlogp1_vec = nlogp_vec(perm_vec);
nlogp1 = squareform(nlogp1_vec);

node_perm_iter(:,iter) = node_perm_idx;

%figure;imagesc(nlogp1);colormap jet;colorbar;snapnow
%% Select the lambda with binarized likelihood, then use greedy_ad

rpmf = [0.2 0.3 0.3 0.2];
r=[0.025 0.01 0.005 0.001];
p0=0.05;
[inlist_blik outlist_blik max_lambda] = greedy_ad_blik(nlogp1,lambda0, r, rpmf, p0);
true_ones = 1:m;
true_node = arrayfun(@(x) find(node_perm_idx==x),true_ones);
true_zeros = (m+1):N;
false_node = arrayfun(@(x) find(node_perm_idx==x),true_zeros);
%tpr = sum(ismember(inlist,true_node))/size(inlist,2);
tpr = sum(ismember(true_node,inlist_blik))/size(true_node,2);
fpr = sum(ismember(false_node,inlist_blik))/size(false_node,2);
tnr = sum(ismember(false_node,outlist_blik))/size(false_node,2);
fnr = sum(ismember(true_node,outlist_blik))/size(true_node,2);
%tpr
%inclusion
result_list(iter,1)=max_lambda(1);
result_list(iter,2)=size(inlist_blik,2);
result_list(iter,3)=tpr;
result_list(iter,4)=fpr;
result_list(iter,5)=tnr;
result_list(iter,6)=fnr;

%% use lambda=1 for reference
[inlist outlist]=greedy_ad(nlogp1,p0,1);
%tpr = sum(ismember(inlist,true_node))/size(inlist,2);
tpr = sum(ismember(true_node,inlist))/size(true_node,2);
fpr = sum(ismember(false_node,inlist))/size(false_node,2);
tnr = sum(ismember(false_node,outlist))/size(false_node,2);
fnr = sum(ismember(true_node,outlist))/size(true_node,2);
result_list_ref(iter,1)=1;
result_list_ref(iter,2)=size(inlist,2);
result_list_ref(iter,3)=tpr;
result_list_ref(iter,4)=fpr;
result_list_ref(iter,5)=tnr;
result_list_ref(iter,6)=fnr;



%% Use Goldberg's Algorithm

[inlist_gold outlist_gold] = goldberg(nlogp1,p0);
tpr = sum(ismember(true_node,inlist_gold))/size(true_node,2);
fpr = sum(ismember(false_node,inlist_gold))/size(false_node,2);
tnr = sum(ismember(false_node,outlist_gold))/size(false_node,2);
fnr = sum(ismember(true_node,outlist_gold))/size(true_node,2);
result_list_gold(iter,1)=1;
result_list_gold(iter,2)=size(inlist,2);
result_list_gold(iter,3)=tpr;
result_list_gold(iter,4)=fpr;
result_list_gold(iter,5)=tnr;
result_list_gold(iter,6)=fnr;



  catch
       error_idx(iter,1)=1;
        
        continue
  end
end

mean(result_list(~error_idx,:),1)
mean(result_list_ref(~error_idx,:),1)
mean(result_list_gold(~error_idx,:),1)
sqrt(var(result_list(~error_idx,:),1))
sqrt(var(result_list_ref(~error_idx,:),1))
sqrt(var(result_list_gold(~error_idx,:),1))


save('res2_09005_15_08.mat','result_list','result_list_ref','result_list_gold','error_idx');
