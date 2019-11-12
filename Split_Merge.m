function [cluster_index, icl, degs2, index_gamma2] = Split_Merge( S_g, cluster_index,K, icl, degs2, index_gamma2, thre_n2)
remain = K-length(icl);
if (remain>0)  %split
    for gi=1:remain
        afa = max(S_g(icl,:));
        ifa = find(afa>thre_n2);
        afa = (afa-min(afa))./(max(afa)-min(afa));
        gamma_new = degs2./afa;
        gamma_new(isnan(gamma_new)) = 0;
        gamma_new(ifa) = 0;
        [~,sort_ga] = sort(gamma_new, 'descend');
        indga = find(gamma_new==Inf);
        [~,indga2] = sort(degs2(indga),'descend');
        sort_ga(1:length(indga)) = indga(indga2);
        for rec = 1:length(icl)
            sort_ga(sort_ga==icl(rec)) = [];
        end
        icl = [icl; sort_ga(1)];
    end
end

%merge
remain =  K-length(icl);
if (remain<0)
    icl = icl(1:K);
end

count_icl = cell(length(icl),1);
for i=1:length(count_icl)
    count_icl{i} = [count_icl{i}, icl(i)];
end
ord = index_gamma2;
for orec=1:length(icl)
    ord(ord==icl(orec)) = [];
end
for i = ord
    temp_check = find(S_g(i,icl) > thre_n2);
    if ~isempty(temp_check)
        [~,max_temp] = max(S_g(i,icl));
        if length(count_icl{max_temp})> length(cluster_index)*0.01
            continue;
        end
        count_icl{max_temp} = [count_icl{max_temp}, i];
    end
end
sr = [];  icl_del  = [];
for lir2=1:length(count_icl)
    if length(count_icl{lir2})==1
        icl_del = [icl_del,lir2];
        continue;
    end
    sr = [sr, mean(S_g(count_icl{lir2},:))'];
end
icl(icl_del) = [];

[~,cluster_index] = max(sr,[],2);

