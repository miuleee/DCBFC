% ---------------------------------
% This script shows how to cluster using DCBFC
% data.S_g: preprocessed data, T*N (T:frame number, N:pixel number)
% data.I:   raw data, T*N (T:frame number, N:pixel number)
% data.mb:  brain template
% ---------------------------------

load('./S_g(0503-1).mat')
S = data.S_g;
mb = data.mb;
% pic = zeros(size(mb)); pic(mb~=0)= data.I(1,:);
% figure; imagesc(pic); colormap gray

% calculate the similarity matrix
tsimi = tic;
r = corrcoef(S);
stime=toc(tsimi);

% DCBFC 
[cluster_index,~,alg_time, mean_cluster_total, icl2] = DCBFC([], S, r, mb,2,0,1,0);

alg_time = alg_time+stime;