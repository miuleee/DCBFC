% ---------------------------------
% This script shows how to obtain different number of clusters for DCBFC.
% data.S_g: preprocessed data, T*N (T:frame number, N:pixel number)
% data.I:   raw data, T*N (T:frame number, N:pixel number)
% data.mb:  brain template
% ---------------------------------

load('./S_g(0503-1).mat')
S = data.S_g;
mb = data.mb;
% calculate the similarity matrix
r = corrcoef(S);

% calculate for different cluster number
oo=[4:20];
tic;
tabl = zeros(length(oo),3);
cf = 1;
cluster_index_tabl = zeros(size(r,2),length(oo));
sout = [];
sout.K = [];
sout.Kclu = [];
sout.para = [];
sout.si = [];
sout.thre = [];
ss = cell(length(oo),1);
for juo=oo
    juo
    [cluster_index,~,alg_time, mean_cluster_total, icl2, sout] = DCBFC(sout, S, r, mb,2,juo,0,1);
    ss{juo-3} = rmfield(sout,{'K','Kclu','para'});
    tabl(cf,1) = max(cluster_index);
    tabl(cf,2) = silhouette_coef(cluster_index,r) ; 
    tabl(cf,3) = alg_time;
    cluster_index_tabl(:,cf) = cluster_index;
    cf = cf+1;
end
toc;

