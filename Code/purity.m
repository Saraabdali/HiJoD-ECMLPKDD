function [n,v]=purity(label, result)
I=size(label,1);
n=0;
for i = 1:I
    if(label(i)==result(i)) 
        n=n+1;
    end
end
    v=n/I;
end
