function [label,idxnl,labelInd] = createSampleLabel3(label3,plabel)
    lnews=size(label3,1);
    labelInd=label3; 
    labelInd(labelInd==0)=-1;
    vectory=[];
    vectoridx=[];
    labelInd_real=find(labelInd==1);
    labelInd_fake=find(labelInd==-1);
    k_real=ceil(plabel*size(labelInd_real,1))
    k_fake=ceil(plabel*size(labelInd_fake,1))
    [y_real,~] = datasample(labelInd_real,k_real,'Replace',false);
    [y_fake,~] = datasample(labelInd_fake,k_fake,'Replace',false);
    label = zeros(lnews,1);
    label(y_real)=1;
    label(y_fake)=-1;
    idxnl=find(label==0);
end
