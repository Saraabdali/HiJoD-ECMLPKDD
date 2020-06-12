function [R]=EstRank_xi(X,alpha,nperm)
r = rank(X);
[~,S,~]=svds(X,r);
lambda=diag(S);
[m,n]=size(X);

R=1;
for i=1:r
    for j=1:nperm
           [m,n]=size(X); 
           X=num2cell(X,2);
           perm_X=cellfun(@(x) x(randperm(n)),X,'UniformOutput',false);
           X_perm=cell2mat(perm_X);
           [~,S,~]=svds(X_perm,r); 
           size(diag(S))
           lambda_perm(j,:)=diag(S);
    end
    clear X_perm
    lambda_perm(i)=prctile(lambda_perm(:,i),100*(1-alpha));
    if(lambda(i)>lambda_perm(i)&& i>R)
      R=i;
    end
end
end
