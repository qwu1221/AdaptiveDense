function [inlist outlist]=greedy_ad(nlogp,p0,lambda)
%%%% Adaptive greedy algorithm with given lambda
%%%%  nlogp: matrix of negative p values
%%%%  p0: screening cutoff in subgraph extraction
%%%%  lambda: parameter in the density definition: |A(S)|/|S|^lambda

W=nlogp;
W1=W;
W(W1<-log(p0))=0;
%figure;imagesc(W)

z1=find(sum(W)>0);
W=W(z1,z1);

density_list = zeros(size(W,1)-2,4);
W_org = W;
for i=1:(size(W,1)-2)
    degs = sum(W, 2);
    posi_degs = find(degs);
%     if size(posi_degs,1)>0
    min_degs=min(degs(posi_degs)); 
    idx = find(degs==min_degs);
    W(idx(1),:)=zeros(1,size(W,1));W(:,idx(1))=zeros(size(W,1),1);
    density_list(i,1)=size(W,1)-i;
    density_list(i,2)=idx(1);
    density_list(i,3)=min_degs;
    density_list(i,4)=(sum(sum(W,2)))/(density_list(i,1))^lambda;
%     else
%         W(1,:)=zeros(1,size(W,1));W(:,1)=zeros(size(W,1),1);
%         density_list(i,1)=size(W,1)-i;
%         density_list(i,2)=1;
%         density_list(i,3)=min_degs;
%         density_list(i,4)=(sum(sum(W,2)))/(density_list(i,1))^lambda;
%     end
    
end

[max_density,idx]=max(density_list(:,4));
density_raw = (sum(sum(W_org,2))/(size(W_org,1))^lambda);
whole = 1:size(W1,1);
q2 = ismember(whole,z1);
z0 = whole(~q2);
if max_density<density_raw
    inlist = z1;
    outlist = z0;
else
outlist = density_list(1:idx,2);
wholeW = 1:size(W,1);
q=ismember(wholeW,outlist);
inlist=z1(~q);
outlist = sort([z1(outlist') z0]);
end
end
