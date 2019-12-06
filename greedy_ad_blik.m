function [inlist outlist max_lambda] = greedy_ad_blik(nlogp1,lambda0, r, rpmf, p0)
%%%%  Adaptive greedy algorithm with lambda selected by integrated likelihood
%%%%  Input: 
%%%%        nlogp1: matrix of negative p values
%%%%        lambda0: vector of possible lambda values
%%%%        r: vector of possible cutoffs
%%%%        rpmf: pmf of possible cutoffs
%%%%        p0: screening cutoff in subgraph extraction
%%%%  Output: 
%%%%         inlist,outlist: vectors of nodes inside and outside the
%%%%         selected subgraph
%%%%         max_lambda: selected lambda values, may be a vector

lik_r0 = zeros(size(lambda0,2),size(r,2));
for i=1:size(r,2)
    r0=r(i);
    para_selection = zeros(size(lambda0,2),5);
    for k=1:size(lambda0,2)
    [inlist outlist]=greedy_ad(nlogp1,p0,lambda0(k));
    W_return = nlogp1([inlist outlist],[inlist outlist]);
    %figure;imagesc(W_return);colormap jet;colorbar;
    W_in = nlogp1(inlist ,inlist);
    %figure;imagesc(W_in);colormap jet;colorbar;
    b_in = W_in;
    b_in(W_in<-log(r0))=0;
    b_in(W_in>-log(r0))=1;
    b_in_vec = squareform(b_in);
    rho1 = mean(b_in_vec);
    loglik_in = sum(b_in_vec)*log(rho1)+(size(b_in_vec,2)-sum(b_in_vec))*log(1-rho1);
    if size(outlist,2)>1
    W_out_vec = [squareform(nlogp1(outlist,outlist)) reshape(nlogp1(inlist,outlist),1,[])];
    b_out_vec = W_out_vec;
    b_out_vec(W_out_vec<-log(r0))=0;
    b_out_vec(W_out_vec>-log(r0))=1;
    rho0 = mean(b_out_vec);
    loglik_out = sum(b_out_vec)*log(rho0)+(size(b_out_vec,2)-sum(b_out_vec))*log(1-rho0);
    elseif size(outlist,2)==1
        W_out_vec = [nlogp1(outlist,outlist) nlogp1(outlist,inlist)];
        b_out_vec = W_out_vec;
        b_out_vec(W_out_vec<-log(r0))=0;
        b_out_vec(W_out_vec>-log(r0))=1;
        rho0 = mean(b_out_vec);
        loglik_out = sum(b_out_vec)*log(rho0)+(size(b_out_vec,2)-sum(b_out_vec))*log(1-rho0);
    else
        loglik_out=0;  
    end
   
    loglik=loglik_in+loglik_out;
    para_selection(k,:)=[lambda0(k) size(W_in,1) loglik_in loglik_out loglik];
    %[lambda0(k) rho1 rho0 size(W_in,1) loglik_in loglik_out loglik]
    
    %max_idx = find(para_selection(:,5)==max(para_selection(:,5)));
    %max_lambda=lambda0(max_idx);
    end
    lik_r0(:,i)=para_selection(:,5);
    
end
lik_int = rpmf*(lik_r0)';
max_idx = find(lik_int==max(lik_int));
max_lambda=lambda0(max_idx);

[inlist outlist]=greedy_ad(nlogp1,p0,max_lambda(1));
end



