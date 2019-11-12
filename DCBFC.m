function [cluster_index,diretime,alg_time, mean_cluster_total, icl2, sout] = DCBFC(sout, S, S_g, mb,percent,K,m_display, m_iter)
    % [cluster_index,~,alg_time, mean_cluster_total, icl2] = DCBFC([],S, r, mb,2,0,1,0);
    tStart = tic;
     r = S_g;
    %% section 2.1.1 Threshold for Determining Cluster Centers
     r=r.*tril(ones(size(r)),-1);
     r=abs(r(r~=0));
     absmeanr = mean(r);
     stdr = std(r);
     thre = stdr+absmeanr;
     
     r = S_g;
     index_record = -1*ones(size(S_g,2),1);
     cluster_index = -1*ones(size(S_g,2),1);  % cluster label     
     gg = 1;  hh = 1;
     loop = 1; u = 0;     
     inc = 1;
     
     r(r<thre & r> -thre)=0;
     
    %% section 2.1.2 Determining Cluster Centers
    dc = size(r,1)*percent/100;
    turni=0;
    icl2 = [];
    diretime = 0;
    while loop
        len = length(find(cluster_index==-1));
        if len<length(cluster_index)*percent/1000
            break;
        end
        u = u+1;
        diret1 = tic;
        % degs: rank degree; the average of the absolute correlation coefficients that are greater than the threshold between the ith pixel and the whole pixels
        lendeg=sum(logical(r));
        degs = zeros(1,len);
        temp = sum(abs(r));
        degs(lendeg>=dc) = temp(lendeg>=dc)./lendeg(lendeg>=dc);
        diretime = diretime+toc(diret1);
        [degs_sorted,orddegs]=sort(degs,'descend');
        if degs_sorted(1)==0
            break;
        end
        rankdeg=zeros(1,len);
        rankdeg(orddegs)=len:-1:1;
        rankdeg = bsxfun(@minus,rankdeg,rankdeg(1,:)');
        rankdeg(rankdeg<0) = 0;
        
        diret1 = tic;        
        % rr: the maximum correlation coefficient between pixel i and the other
        % pixels whose rank degree values are higher than pixel i
        rr=max(r.*logical(rankdeg),[],2);
        diretime = diretime+toc(diret1);
        clear rankdeg
        
        % Normalized to [0,1]
        degs = (degs-min(degs))./(max(degs)-min(degs));
        rr = (rr-min(rr))./(max(rr)-min(rr));
        rr(orddegs(1)) = 0;
        
        % calculate the gamma value for each pixel
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
        if turni==0
            turni=1;
            index_gamma2=index_gamma;
            orddegs2 = orddegs;
            degs2 = degs;
        end
        %     figure;stem(gamma_sort(1:end));
        
        % select the candidate center pixels
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
        
        % select the center pixels from the candidate center pixels
        icl(1:NCLUST) = index_gamma(1:NCLUST);
        r_central = S_g(temp_cluster(icl),temp_cluster(icl)); 
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
        
        % the pixels which are similar to center pixels are removed from current set {D}        
        ord = index_gamma;%
        for orec=1:length(icl)
            ord(ord==icl(orec)) = [];
        end
        count_icl = cell(length(icl),1);
        for i=1:length(count_icl)
            count_icl{i} = [count_icl{i}, icl(i)];
        end
        for i = ord 
            temp_check = find(S_g(temp_cluster(i),temp_cluster(icl)) > thre); 
            if ~isempty(temp_check)
                [~,max_temp] = max(S_g(temp_cluster(i),temp_cluster(icl))); 
                if length(count_icl{max_temp})> length(ind)*0.01
                    continue;
                end
                count_icl{max_temp} = [count_icl{max_temp}, i];
            end
        end
        
        max_index=zeros(size(r,1),1);
        for i=1:length(count_icl)
            count_icl2 = count_icl{i}; 
            subthre = S_g(temp_cluster(count_icl2),cluster_index == -1)-thre;
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
        icl2 = [icl2; temp_cluster(icl)];
        
        % form new similarity matrix
        r = r(ind==-1,ind==-1);
        temp = index_record(index_record == -1);
        temp(ind~=-1) = u;
        index_record(index_record == -1) = temp;
        
        len2 = length(find(cluster_index==-1));
        if len2==len
            break;
        end
    end

    %% section 2.1.3 Clustering All Pixels to Form the Functional Connection Modules
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
        temp_check = find(S_g(i,icl) > thre); 
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
        if isempty(S)
            sr = [sr, mean(S_g(count_icl{lir2},:))'];
        else
            sr = [sr, mean(S(:,count_icl{lir2}),2)];
        end
    end
    icl(icl_del) = [];
    icl2(icl_del) = [];
    
    if ~isempty(S)
        sr = corr(S,sr);
    end
    [~,cluster_index] = max(sr,[],2);
            
    if K>0 % for different cluster number 
        [cluster_index, icl, degs2, index_gamma2, thre, sout] = diffnum(sout, S_g, percent, cluster_index,K, icl2, degs2, index_gamma2, thre, m_iter,absmeanr,stdr);             
    end
    
    alg_time=toc(tStart);
    
    sout.icl = icl;
    sout.cluster_index = cluster_index;
    sout.degs2 = degs2;
    sout.index_gamma2 = index_gamma2;
    sout.thre = thre;
    
    K = max(cluster_index);
    mean_cluster_total = zeros(K,size(S_g,1));
    for i=1:K
        mean_cluster_total(i,:) = mean(S_g(:,cluster_index==i),2);
    end


%% display results of DCBFC
index_remain = 1:K;
a = 1:length(index_remain);
a = a(randperm(length(a)));
if(m_display)
    li=zeros(size(mb));
    li(mb~=0)=cluster_index;
    pic=li;
    pic1=medfilt2(pic,[5 5]);
    temp=pic1;
    temp(temp==pic(1,1))=0;
    temp(temp~=0)=1;
    pic1=pic1.*temp;
    pic1 = round(pic1);
    figure; imagesc(pic1); colormap jet
    %title('DCBFC','Fontsize',14,'FontWeight','bold','FontName','Times New Roman')
    axis equal
    axis([1,size(mb,1),1,size(mb,2)])
    set(gca,'xtick',[],'ytick',[]);
    colormap jet
    [row,col] = find(mb==0);
    hold on
    stem(col,row,'MarkerSize',5,'MarkerFaceColor',[0 0 0],...
        'MarkerEdgeColor',[0 0 0],...
        'Marker','square',...
        'LineStyle','none',...
        'Color',[0 0 0]);       
end

if (m_display)
    li=zeros(128);
    [row,col] = find(mb==0);
    min_m = min(mean_cluster_total(:));
    max_m = max(mean_cluster_total(:));
    for i=1:size(mean_cluster_total,1)
        li(mb~=0)=mean_cluster_total(i,:);
        li = medfilt2(li,[5,5]);
        figure;imagesc(li);%title(num2str(i));
        axis equal
        axis([1,size(mb,1),1,size(mb,2)])
        set(gca,'xtick',[],'ytick',[],'CLim',[-0.64,0.85]);
        %set(gca,'xtick',[],'ytick',[],'CLim',[min_m,max_m]);
        colormap jet
        hold on
        stem(col,row,'MarkerSize',8,'MarkerFaceColor',[0 0 0],...
            'MarkerEdgeColor',[0 0 0],...
            'Marker','square',...
            'LineStyle','none',...
            'Color',[0 0 0]);
    end
end
    