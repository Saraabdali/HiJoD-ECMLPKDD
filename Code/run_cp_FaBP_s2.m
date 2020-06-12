function [avgl,prec_avgl,f1_avgl, rec_avgl]= run_cp_FaBP(C,maxnn,it,pl,label3)
[m,v,~,prec,rec,f1,~,vl,~,precl,recl,f1l,specl,spec]=test_FaBP(C,maxnn,it,label3,pl);
disp('FaBP done!');
spec_avg=mean(spec);
prec_avg=mean(prec);
rec_avg=mean(rec);
f1_avg=mean(f1);
avg=mean(v);
vvar = var(v,1); %normalized by N
sd=std(v,1);
prec_avgl=mean(precl);
rec_avgl=mean(recl);
spec_avgl=mean(specl);
f1_avgl=mean(f1l);
avgl=mean(vl);
vvarl = var(vl,1); %normalized by N
sdl=std(vl,1);
clear C maxnn it pl label3 filename;
end

