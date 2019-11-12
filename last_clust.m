function sout = last_clust(sout,  S_g, degs2, K)
if ~isempty(sout)
    if (~isempty( sout.thre))  %split cluster_old
        if max(sout.cluster_index)~= K-1
            return;
        end
        afa = max(S_g(sout.icl ,:));
        afa = (afa-min(afa))./(max(afa)-min(afa));
        gamma_new = sout.degs2./afa;
        gamma_new(isnan(gamma_new)) = 0;
        [~,sort_ga] = sort(gamma_new, 'descend');
        for rec = 1:length(sout.icl )
            sort_ga(sort_ga==sout.icl (rec)) = [];
        end
        indga = find(gamma_new==Inf);
        [~,indga2] = sort(degs2(indga),'descend');
        sort_ga(1:length(indga)) = indga(indga2);
        sout.icl  = [sout.icl ; sort_ga(1)'];
        
        count_icl  = cell(length(sout.icl ),1);
        for i=1:length(count_icl)
            tw = sout.icl (i);
            count_icl{i} = [count_icl{i}, tw];
        end
        
        ord = sout.index_gamma2 ;
        for orec=1:length(sout.icl )
            ord(ord==sout.icl (orec)) = [];
        end
        for i = ord
            temp_check = find(S_g(i,sout.icl ) > sout.thre  );
            if ~isempty(temp_check)
                [~,max_temp] = max(S_g(i,sout.icl ));
                if length(count_icl {max_temp})> length( sout.cluster_index )*0.01
                    continue;
                end
                count_icl{max_temp} = [count_icl{max_temp}, i];
            end
        end
        sr = [];   icl_del  = [];
        for lir2=1:length(count_icl)
            if length(count_icl{lir2})==1
                icl_del = [icl_del,lir2];
                continue;
            end
            sr = [sr, mean(S_g(count_icl{lir2},:))'];
        end
        sout.icl(icl_del) = [];
        [~, sout.cluster_index ] = max(sr,[],2);
        
    end
end