function linknews = knnInd(nn)
    ln=size(nn,2);
    lnews=size(nn,1);
    nn
    linknews = zeros((ln-1)*lnews,2);
    l=1;
    for i=1:lnews
        for j=2:ln
            linknews(l,:)=[i,nn(i,j)];
            l=l+1;
        end
    end
end
