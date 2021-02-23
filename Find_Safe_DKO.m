function [dko_safe,dko_toxic,grRateWT_ctl] = Find_Safe_DKO(DKO_all_inf,model_ctl,thershold)
    %% Find DKO genes in infected model that is safe on healthy model
    % build a 1d gene list of the infected DKO
    [r,c] = size(DKO_all_inf)
    dko_list = reshape(table2array(DKO_all_inf),[r*c,1])
    dko_list = table2cell(unique( DKO_all_inf))
    dko_list = string(dko_list)
    % find the common between the infected DKO and all genes in normal model
    dko_common = intersect(dko_list,model_ctl.genes)
    if numel(dko_common)==0
        dko_unq = DKO_all_inf;
        grRateWT_ctl = 0;
        dko_safe = dko_list;
        dko_toxic =  0;
    else
        [grRatio_ctl, ~, grRateWT_ctl]= doubleGeneDeletion(...
            model_ctl, 'FBA',dko_common');
    
        %create a table of DKO in normal model 
        [r,c]=find(grRatio_ctl <thershold);
        r2=r(r~=c);
        c2=c(r~=c);
        DKO_ctl=[dko_common(r2),dko_common(c2)];
        %DKO_ctl=unique(sort(DKO_ctl')', 'rows');
        %DKO_all_inf_ = [DKO_all_inf;DKO_all_inf(:,1),DKO_all_inf(:,2)]
        DKO_ctl_pairs = strcat(DKO_ctl(:,1),'-',DKO_ctl(:,2))

        DKO_all = table2cell(DKO_all_inf);
        DKO_all_inf_pairs =  strcat(DKO_all(:,1),'-',DKO_all(:,2));
        % Safe drug pairs
        dko_unq_ = setdiff(DKO_all_inf_pairs,DKO_ctl_pairs);
        % Toxic drug pairs
        dko_shared_ = intersect(DKO_all_inf_pairs,DKO_ctl_pairs);
        dko_unq = zeros(numel(dko_unq_),2);
        for i=1:numel(dko_unq_)
            dko_unq(i,:) = strsplit(dko_unq_(i),'-');
        end
        dko_shared = zeros(numel(dko_shared_),2);
        for i=1:numel(dko_shared_)
            dko_shared(i,:) = strsplit(dko_shared_(i),'-');
        end
        dko_toxic = dko_shared;
        dko_safe = dko_unq;
    end
end

