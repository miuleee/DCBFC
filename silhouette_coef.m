function [SI, mean_inter, mean_intra, SI_cluster] = silhouette_coef(cluster,r)
SI = [];
ssind = find(cluster~=0);
cluster = cluster(ssind);
r = r(ssind,ssind);
n_samples_per_label = zeros(1,max(cluster));
r = 1-r;
for i=1:max(cluster)
    n_samples_per_label(i) = length(find(cluster==i));
end
intra_clust_dists = zeros(length(cluster),1);
inter_clust_dists = inf*ones(length(cluster),1);


for curr_label=1:max(cluster)
    mask = find(cluster==curr_label);
    current_distances = r(mask,:);
    intra_clust_dists(cluster==curr_label) = sum(current_distances(:,mask))/(n_samples_per_label(curr_label)-1);

    for other_label = 1:max(cluster)
        if other_label ~= curr_label
            other_mask = find(cluster == other_label);
            other_distances = mean(current_distances(:, other_mask),2);
            inter_clust_dists(mask) = min([inter_clust_dists(mask),other_distances],[],2);
        end
    end
end


SI = inter_clust_dists - intra_clust_dists;
SI = SI./max([inter_clust_dists , intra_clust_dists],[],2);
del_ind = find(n_samples_per_label==1);
if ~isempty(del_ind)
    SI(cluster(del_ind)) = 0;
end
SI(isnan(SI))=[];

SI_cluster = [];
for i=1:max(cluster)
    SI_cluster(i) = mean(SI(cluster==i));
end

SI = mean(SI);

mean_inter = mean(1-inter_clust_dists);
mean_intra = mean(1-intra_clust_dists);






















































