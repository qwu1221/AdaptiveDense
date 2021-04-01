function [p T0 Tm_vec] = perm_label(A, B, W, inlist,lambda_vec,r,rpmf,p0,V)

N = size(W);
n1 = size(A,2);
n2 = size(B,2);
C = [A B];


T0 = cal_loglik(W,inlist,r,rpmf);


Tm_vec=[];
error_idx  = zeros(V,1);
for i=1:V
    try
    label_perm = randperm(n1+n2);
    A_new = C(:,label_perm(1:n1));
    B_new = C(:,label_perm((n1+1):(n1+n2)));
    for j=1:size(A,1)
    [~,p_perm(j)]=ttest2(A_new(j,:),B_new(j,:));
    end
    P_perm=squareform(p_perm);
    P_perm = P_perm+eye(N);
    Wperm=-log(P_perm);
    
    
    [inlist_perm outlist_perm max_lambda max_lik] = greedy_ad_blik(Wperm,lambda_vec, r, rpmf, p0);

    
    Tm = cal_loglik(Wperm,inlist_perm,r,rpmf);
    Tm_vec(i) = Tm;
    catch
       error_idx(i)=1;
       Tm_vec(i)=0;
        continue
  end
end

p = sum(Tm_vec>T0)/sum(error_idx==0);
Tmax_5prct = prctile(Tm_vec(Tm_vec>0),95);
end
