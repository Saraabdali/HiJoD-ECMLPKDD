function [R]=EstRank(C,alpha,nperm)
X=vertcat(C{:});
r = rank(X);
[~,S,~]=svds(X,r);
lambda=diag(S);
lambda=abs(lambda);
[a,b]=size(X);

R=1;
for i=1:r
    for j=1:nperm
        for k=1:length(C)
           X=C{k};
           [m,n]=size(X); 
           X=num2cell(X,2);
           perm_X=cellfun(@(x) x(randperm(n)),X,'UniformOutput',false);
           X_perm{k}=cell2mat(perm_X);
         end
         X_perm=vertcat(X_perm{:}); 
         [~,S,~]=svds(X_perm,r); 
         lambda_perm(j,:)=abs(diag(S));
         clear X_perm
    end
    lambda_perm(i)=prctile(lambda_perm(:,i),100*(1-alpha));
    if(lambda(i)>lambda_perm(i)&& i>R)
      R=i;
    end
end
end
