function [cluster_index, icl, degs2, index_gamma2, thre, sout] = diffnum(sout,S_g, percent,cluster_index,K, icl2, degs2, index_gamma2, thre, m_iter,absmeanr,stdr)
% Initialize the thre, 81 corresponding to thre=mean+std.
f = Inf;
mean_num2 = [linspace(2,max(cluster_index),81)';linspace(max(cluster_index),24,41)'];
mean_num2(81)=[];
mean_num2 = mean_num2/mean_num2(81);
%oo = absmeanr+ [-1:0.025:2].*stdr;

%% sout
try 
    sout.K;
    if ~isempty(find(sout.K==K))
        m_iter = 0;
        degs2 = sout.para{K}.degs2;
        icl2 = sout.para{K}.icl;
        index_gamma2 = sout.para{K}.indexgamma2;
        thre = sout.para{K}.thre;
        cluster_index = sout.Kclu{K};
    end
end

%% adaptive thre
remain = K-length(icl2);
thre_new = thre;
iternum = 10; iterec=2;
cluster_index_new = cluster_index;
recremain = abs(remain); 
thre_n2 = thre_new;
icl = icl2; current_ind = 1;
while (m_iter~=0 && remain~=0 && iternum>0 )
    mc = length(icl);
    if (iternum==10)
        opind = [-1:0.025:2]';
        [~,opindmin] = min(abs(K/mc-mean_num2));
        thre_new = absmeanr+opind((opindmin(1)))*stdr+(current_ind-1)*0.025*stdr;
        if (thre_new >0.99 || thre_new<-0.99)
            break;
        end
        [L_alg, degs_temp, icl_temp, indexgamma_temp] = THRE_ALGnew(S_g,percent,thre_new);
        cluster_index = L_alg;
        if ~isempty(sout)
            if(isempty(find(sout.K==max(L_alg))))
                sout.K = [sout.K, max(L_alg)];
                sout.Kclu{max(L_alg)} = L_alg;
                sout.para{max(L_alg)}.degs2 = degs_temp;
                sout.para{max(L_alg)}.icl = icl_temp;
                sout.para{max(L_alg)}.indexgamma2 = indexgamma_temp;
                sout.para{max(L_alg)}.thre = thre_new;
            end
        end
    else
        thre_new2(1) = thre_n2+(K-mc)*current_ind*0.025*stdr;
        thre_new2(2) = thre_n2-(K-mc)*current_ind*0.025*stdr;
        for rec=1:2
            if (thre_new2(rec) >0.99 || thre_new2(rec)<-0.99)
                break;
            end
            [L_alg2{rec}, degs_temp2{rec}, icl_temp2{rec}, indexgamma_temp2{rec}] = THRE_ALGnew(S_g,percent,thre_new2(rec));
            if ~isempty(sout)
                if(isempty(find(sout.K==max(L_alg2{rec}))))
                    sout.K = [sout.K, max(L_alg2{rec})];
                    sout.Kclu{max(L_alg2{rec})} = L_alg2{rec};
                    sout.para{max(L_alg2{rec})}.degs2 = degs_temp2{rec};
                    sout.para{max(L_alg2{rec})}.icl = icl_temp2{rec};
                    sout.para{max(L_alg2{rec})}.indexgamma2 = indexgamma_temp2{rec};
                    sout.para{max(L_alg2{rec})}.thre = thre_new2(rec);
                end
            end
        end
        [~,rind] = min([abs(K-length(icl_temp2{1})), abs(K-length(icl_temp2{2}))]);
        cluster_index = L_alg2{rind};
        degs_temp =degs_temp2{rind};
        icl_temp =icl_temp2{rind};
        indexgamma_temp = indexgamma_temp2{rind};
        thre_new = thre_new2(rind);
    end
    
    if isempty(cluster_index)
        break;
    end
    
    remain =K-length(icl_temp);
    if (abs(remain)<=recremain)
        recremain = abs(remain);
        degs2 = degs_temp;
        icl = icl_temp;
        index_gamma2 = indexgamma_temp;
        thre_n2 = thre_new;
        cluster_index_new = cluster_index;
        current_ind = 1;
    elseif (abs(remain)>recremain)
        current_ind = current_ind+1;
    end
    iternum = iternum-1;
    iterec = iterec+1;
end

recremain = K-length(icl);
cluster_index = cluster_index_new;

if recremain~=0
    [cluster_index, icl, degs2, index_gamma2] = Split_Merge( S_g, cluster_index,K, icl, degs2, index_gamma2, thre_n2);
end

% correct with last clust
sout = last_clust(sout, S_g, degs2, K);
if ~isempty(sout)
    if (~isempty(sout.thre))
        sout.si = silhouette_coef(sout.cluster_index,S_g);
        si = silhouette_coef(cluster_index,S_g);
        if sout.si>si
            icl = sout.icl;
            cluster_index = sout.cluster_index;
            degs2 = sout.degs2 ;
            index_gamma2 =  sout.index_gamma2;
        end
    end
end
thre  = thre_n2;

