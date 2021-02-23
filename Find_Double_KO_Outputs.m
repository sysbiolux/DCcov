    %(c) Ali Kishk, Maria Pires Pacheco
    % 09 April 2020
    % Luxembourg University
function [DKO_all,DKO_non_ess,DKO_syn,DKO_both] = Find_Double_KO_Outputs(...
    grRatioDble,all_genes,SKO_genes)
%Find_Double_KO_Outputs Extracts 3 outputs of the grRatio in double ko
    %Inputs: 
    % grRatioDble: n*n double grRatio for Double KO
    % all_genes: n*1 string gene names of the model
    % SKO_genes: m*1 cell of Single KO gene names
    
    %% DKO_all: are couples of genes that have a flux of 0 when they are knock out together
    %% DKO_non_ess: are couple of two non essential genes that have a flux below the threshold
    %% DKO_syn are: couple of essential and non essential have a flux below the flux of essential gene alone
    %% DKO_all are included in DKO_non_ess
    %% DKO_both both of non essential and synergistic
    sko_idx = find(ismember(all_genes, SKO_genes));
    non_essential=setdiff(1:numel(all_genes), sko_idx);
    [r,c]=find(grRatioDble(non_essential,non_essential) ==0);
    r2=r(r~=c);
    c2=c(r~=c);

    DKO_all=[non_essential(r2)',non_essential(c2)'];
    DKO_all=unique(sort(DKO_all')', 'rows');
    DKO_all=table(all_genes(DKO_all(:,1)),all_genes(DKO_all(:,2)));
    DKO_all.Properties.VariableNames = [{'c'}    {'r'}]
    [r,c]=find(grRatioDble == 0);
    total_ko_genes=all_genes(r(r==c))

    sko_idx = find(ismember(all_genes, SKO_genes));
    % t3 (essential genes with a non zero flux) in the model genes
    t3 = setdiff(SKO_genes,total_ko_genes);

    % indices of t3 (essential genes with a non zero flux) in the model genes
    t3_idx= find(ismember(all_genes,t3));
    % Get all pairs of Double KO genes
    [r,c] =  find(grRatioDble <= 0.1);
    dbl_aal_pairs = table(all_genes(r),all_genes(c));

    %%Creating a table of non essential genes
    not_an_essential=setdiff(1:numel(all_genes),sko_idx);
    [c,r]=find(grRatioDble(not_an_essential,not_an_essential)< 0.1);
    % non essential genes from non essential index
    c = all_genes(not_an_essential(c)');
    r = all_genes(not_an_essential(r)');
    DKO_non_ess=table(c,r);
    % remove dubulicated pairs
    DKO_non_ess = table2cell(DKO_non_ess)
    DKO_non_ess = string(DKO_non_ess)
    DKO_non_ess=unique(sort(DKO_non_ess')', 'rows');
    % convert the table to char table
    c = DKO_non_ess(:,1);
    r = DKO_non_ess(:,2);
    DKO_non_ess=table(cellstr(c),cellstr(r));
    DKO_non_ess.Properties.VariableNames = [{'c'}    {'r'}]
    
    %% Finding synergistic gene pairs %%%
    % Getting the index of Synergistic DKO genes pairs
    DKO_syn=zeros(10000, 2);
    match = 0;
    k = 0;
    total_ko_genes_Id=find(ismember(all_genes,total_ko_genes));
    for i=1:numel(t3_idx)
        % check if the KO of the gene put he flux to 0
        if ~isempty(setdiff(t3_idx(i),total_ko_genes_Id))

            rate_single_del=grRatioDble(t3_idx(i), t3_idx(i));

            [match, ~]=find(grRatioDble(:, t3_idx(i))< rate_single_del*0.9);
            k = k+numel(match);
            match2= intersect(match, sko_idx);
            match = setdiff(match,sko_idx);
            if ~ isempty(match)
                DKO_syn(k:k+numel(match)-1,1)=match;
                DKO_syn(k:k+numel(match)-1,2)= t3_idx(i);
            end
            % check if two essential genes can bring the flux down
            if ~ isempty(match2)
                a=setdiff(match2, total_ko_genes_Id);
                if ~isempty(a)
                    clc
                end
            end

        end
    end
    
    DKO_syn = DKO_syn(all(DKO_syn,2),:);
    c1 = DKO_syn(:,1);
    c2 = DKO_syn(:,2);
    c = all_genes(c1');
    r = all_genes(c2');
    DKO_syn=table(c,r);    
    % Creating DKO_both for both the non essenial and synergestic
    %DKO_both = array2table([DKO_non_ess;DKO_syn]);
    
    r = size(DKO_all,1);
    if ~isempty(DKO_all) & ~isempty(DKO_non_ess) &  ~isempty(DKO_syn)
        DKO_both = [DKO_all;DKO_non_ess;DKO_syn]
        DKO_both = DKO_both(r+1:end,:)
    elseif ~isempty(DKO_all) & ~isempty(DKO_non_ess)
         DKO_both = [DKO_all;DKO_non_ess]
         DKO_both = DKO_both(r+1:end,:)
    elseif ~isempty(DKO_all) & ~isempty(DKO_syn)
          DKO_both = [DKO_all;DKO_syn]
          DKO_both = DKO_both(r+1:end,:)
    elseif ~isempty(DKO_non_ess) &  ~isempty(DKO_syn)
        DKO_both = [DKO_non_ess;DKO_syn]
    elseif  ~isempty(DKO_non_ess)
        DKO_both = DKO_non_ess
    else
        DKO_both = DKO_syn
    end
end