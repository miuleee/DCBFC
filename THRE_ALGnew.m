function [cluster_index,degs2, icl2, index_gamma2]  = THRE_ALGnew(S_g,percent,thre)
%tic;
     r = S_g;
     
     index_record = -1*ones(size(S_g,2),1);
     cluster_index = -1*ones(size(S_g,2),1);       
     cluster_index2 = -1*ones(size(S_g,2),1);  
     degs2 = -1*ones(size(S_g,2),1); 
     index_gamma2 = -1*ones(size(S_g,2),1); 
     
     gg = 1;  hh = 1;
     loop = 1; u = 0;
     inc = 1;
     r2 =r;   
     r(r<thre & r> -thre)=0;
     
     %% 
     dc = size(r,1)*percent/100;
     
     turni=0;
     icl2 = [];
     while loop
         len = length(find(cluster_index==-1));
         if len<20
             break;
         end
         u = u+1;
         
         lendeg=sum(logical(r));
         degs = zeros(1,len);
         temp = sum(abs(r));
         degs(lendeg>=dc) = temp(lendeg>=dc)./lendeg(lendeg>=dc);
         
         [degs_sorted,orddegs]=sort(degs,'descend');
         if degs_sorted(1)==0
             break;
         end
         
         rankdeg=zeros(1,len);
         rankdeg(orddegs)=len:-1:1;
         rankdeg = bsxfun(@minus,rankdeg,rankdeg(1,:)');
         rankdeg(rankdeg<0) = 0;
         rr=max(r.*logical(rankdeg),[],2);
         clear rankdeg
         
         degs = (degs-min(degs))./(max(degs)-min(degs));
         
         rr = (rr-min(rr))./(max(rr)-min(rr));
         rr(orddegs(1)) = 0;
         
         gamma = zeros(1,len);
         ind = -1*ones(1,size(r,1));
         ind2 = -1*ones(1,size(r,1));
         gamma=degs./rr';
         gamma(isnan(gamma))=0;
         if isnan(gamma)
             break;
         end
         gamma(gamma==-Inf) = Inf;
         [gamma_sort,index_gamma] = sort(gamma,'descend');
         %     figure;stem(gamma_sort(2:end));
         if turni==0
             turni=1;
             index_gamma2=index_gamma;
             orddegs2 = orddegs;
             degs2 = degs;
         end
         
         ga1 = find(gamma_sort~=Inf,1);
         ga2 = find(gamma_sort~=0,1,'last');
         esp = (gamma_sort(ga1)-1)/exp(1)+1;
         ff = find(gamma_sort>esp,1,'last');
         if(ff>dc)
             ff=ga1;
         end
         
         if isempty(ff)
             break;
         end
         NCLUST = ff(1);
         temp_cluster = find(cluster_index == -1);
         clear icl
         icl(1:NCLUST) = index_gamma(1:NCLUST);
         r_central = r2(temp_cluster(icl),temp_cluster(icl)); 
         r_central(1:size(r_central,1)+1:end)=0;
         subcount = length(icl);
         addcount = 1;
         delind = [];
         while(subcount)
             tempind = find(r_central(addcount,:)>thre);
             tempind(tempind<=addcount) = [];
             delind = [delind, tempind];
             addcount = addcount+1;
             subcount = subcount-1;
         end
         delind = unique(delind);
         icl(delind) = [];
         NCLUST = length(icl);
         
         ord = index_gamma;
         for orec=1:length(icl)
             ord(ord==icl(orec)) = [];
         end
         count_icl = cell(length(icl),1);
         for i=1:length(count_icl)
             count_icl{i} = [count_icl{i}, icl(i)];
         end
         for i = ord
             temp_check = find(r2(temp_cluster(i),temp_cluster(icl)) > thre);
             if ~isempty(temp_check)
                 [~,max_temp] = max(r2(temp_cluster(i),temp_cluster(icl)));
                 if length(count_icl{max_temp})> length(ind)*0.01
                     continue;
                 end
                 count_icl{max_temp} = [count_icl{max_temp}, i];
             end
         end
         
         max_index=zeros(size(r,1),1);
         for i=1:length(count_icl)
             count_icl2 = count_icl{i};
             subthre = r2(temp_cluster(count_icl2),cluster_index == -1)-thre;
             subthre(subthre<0) = 0;
             max_index(sum(subthre)>0)=i;
         end
         
         ic = cell(1,length(icl));
         for i=1:NCLUST
             ic{i} = find(max_index==i);
             ind(ic{i}) = gg;
             ind2(icl(i)) = gg;     gg=gg+1;
         end
         
         cluster_index(cluster_index == -1) = ind;
         cluster_index2(temp_cluster) = ind2;
         icl2 = [icl2; temp_cluster(icl)];
         
         r = r(ind==-1,ind==-1);
         temp = index_record(index_record == -1);
         temp(ind~=-1) = u;
         index_record(index_record == -1) = temp;
         
         len2 = length(find(cluster_index==-1));
         if len2==len
             break;
         end
     end
     
     icl = icl2;
     if isempty(icl)
         cluster_index = [];
         return
     end
          
     ord = index_gamma2;
     for orec=1:length(icl)
         ord(ord==icl(orec)) = [];
     end
     count_icl = cell(length(icl),1);
     for i=1:length(count_icl)
         count_icl{i} = [count_icl{i}, icl(i)];
     end
     for i = ord
         temp_check = find(r2(i,icl) > thre);
         if ~isempty(temp_check)
             [~,max_temp] = max(r2(i,icl));
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
     icl2(icl_del) = [];
     
     [~,cluster_index] = max(sr,[],2);
  
  
