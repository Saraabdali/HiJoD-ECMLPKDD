%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ECML/PKDD 2020%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Semi-Supervised Multi-aspect Detection of Misinformation using Hierarchical Joint Decomposition%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sara Abdali, Neil Shah, Evangelos E.Papalexakis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%R_Tags1,R_TTA1,R_HTA1 are initial ranks for TAGS, TTA and HTA respectively and nn is the number of neighbor

function HiJoD(R_Tags1,R_TTA1,R_HTA1,nn)
addpath('./tensor_toolbox_2.6');
addpath('./tensor_toolbox_2.6/met');
pl=[0.1;0.2;0.3;0.4];
for k=1:length(pl)
    for h=1:100
        load(strcat('./newsTagsMatrix',num2str(h),'.mat')); %loading TAGS matrix
        label3=tagsHtmlDomain.label3;
        load(strcat('./HTA_',num2str(h),'.mat'));  %loading HTA (Hashtag-Term-Article) Tensor
        hxwxn=wxwxn;
        load(strcat('./TTA_',num2str(h),'.mat'));  %loading TTA (Term-Term-Article) Tensor
        wxwxn=wxwxn;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%Level-1 Decomposition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic
        hxwxn=cp_als(hxwxn,R_HTA1);
        wxwxn=cp_als(wxwxn,R_TTA1);
        [Tags,~,~]=svds(mTags,R_Tags1);
        toc
        Tol=1e-12;
        wxwxn=wxwxn.U{3};
        wxwxn=wxwxn';
        hxwxn=hxwxn.U{3};
        hxwxn=hxwxn';
        Tags=Tags';
        X=[Tags;wxwxn;hxwxn];
        C{1}=Tags;
        C{2}=wxwxn;
        C{3}=hxwxn;
        [wht,n]=size(X);
        Res=zeros(wht,n);
        Err=Inf;
        R_final=rank(X);
        R_final_prev=1;
        tic
        %%%%%%%%%%%%%%%%%%%%%%%%%%Level-1 Decomposition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        while(R_final~=R_final_prev || Err>Tol)
              %calculating J
              R_final_prev=R_final;
              R_final=EstRank(C,0.05,25)
              [U,S,V] = svds(X,R_final);
              I = eye(n,n);
              Proj=I-V*V';
              J=U*S*V';
              J_Tags=J(1:R_Tags1,:);
              Individual_Tags=Tags-J_Tags;
              J_wxwxn=J(R_Tags1+1:R_Tags1+R_TTA1,:);
              Individual_wxwxn=wxwxn-J_wxwxn;
              J_hxwxn=J(R_Tags1+R_TTA1+1:R_Tags1+R_TTA1+R_HTA1,:);
              Individual_hxwxn=hxwxn-J_hxwxn;
              C{1}=Individual_Tags*Proj;
              R_Tags=EstRank(C,0.05,40); 
              [U,S,V] = svds(Individual_Tags*Proj,R_Tags); 
              A_Tags=U*S*V';
              Tags_new=Tags-A_Tags;
              C{1}=Individual_wxwxn*Proj;
              R_TTA=EstRank(C,0.01,40);
              [U,S,V] = svds(Individual_wxwxn*Proj,R_TTA);
              A_wxwxn=U*S*V';
              wxwxn_new=wxwxn-A_wxwxn;
              C{1}=Individual_hxwxn*Proj;
              R_HTA=EstRank(C,0.05,40); 
              [U,S,V] = svds(Individual_hxwxn*Proj,R_HTA);
              A_hxwxn=U*S*V';  
              hxwxn_new=hxwxn-A_hxwxn;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              X_new=[Tags_new;wxwxn_new;hxwxn_new];
              C{1}=Tags_new;
              C{2}=wxwxn_new;
              C{3}=hxwxn_new;
              Res=X_new-X;
              Err=norm(Res, 2);
              X=X_new;
              Tags=X(1:R_Tags1,:);
              wxwxn=X(R_Tags1+1:R_Tags1+R_TTA1,:);
              hxwxn=X(R_Tags1+R_TTA1+1:R_Tags1+R_TTA1+R_HTA1,:);
      end
   C_final=J';
   toc
   save('X.mat','C_final');
   label3(isnan(label3))=0;
  [avgl(h),prec_avgl(h),f1_avgl(h), rec_avgl(h)]= run_cp_FaBP_s2(C_final,nn,100,pl(k),label3);                        
end
 file_name= strcat('./JIVE_',num2str(pl(k)*100),'rank_',num2str(R_final),'.mat');
 save(file_name,'-v7.3','avgl','prec_avgl','f1_avgl', 'rec_avgl');
end
