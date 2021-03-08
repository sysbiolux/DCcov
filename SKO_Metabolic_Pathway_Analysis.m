%% Metabolic Pathway Analysis of the SKO genes

%Load Recon3D model
load('./Recon3DModel_301.mat')
%Removing gene version from the model
genes = Recon3DModel.genes;
for i=1:numel(genes)
    x = strsplit(table2array(genes(i,1)),'.');
    genes(i,1) = cellstr(x(1,1));
end
Recon3DModel.genes = genes;

recon_model = Recon3DModel;
load('dico_recon.mat')
idx = find(contains(string(recon_model.subSystems),"Pyrimidine"));
[string(recon_model.subSystems(idx)),recon_model.subSystems(idx)]
%% Lung study 1
sko_df = readtable('KO_data/severity_SKO_All.csv');

grps= unique(sko_df.Series);
% define the phenotype
T_SKO = cell2table(cell(0,3));
T_SKO.Properties.VariableNames = [{'Pathway','NB','Phenotype'}];


for i=1:numel(grps)
    grp = grps{i};
    %convert the 2d table to 1d list of entreez ids
    sko = sko_df(find(ismember(sko_df.Series,grp)),3:end);
    sko = string(reshape(table2cell(sko),[1,size(sko,1)*size(sko,2)]));
    sko = rmmissing(sko);
    sko = strjoin(sko,',');
    sko = unique(strsplit(sko,','));
    sko(find(ismember(sko,""))) = [];
    % Find shared genes between RECON model genes and SKO
    genes_metabolic=find(ismember(recon_model.genes,sko));
    [r_genes,~]=find(recon_model.rxnGeneMat(:, genes_metabolic));
    genes_ss=recon_model.subSystems(r_genes);

    [pathways_genes, ~, ub] = unique(string(genes_ss));
    path_genes = histc(ub, 1:length(pathways_genes))*1;
    grp = repmat(string(grp),size(path_genes,1),1);
    T_genes = table(pathways_genes, path_genes,grp);
    T_genes.Properties.VariableNames={'Pathway', 'NB','Phenotype'};
    T_SKO = [T_SKO;T_genes];
end

writetable(T_SKO,'KO_data/severity_SKO_RECON_Pathways.csv');

%% Lung study 2
sko_df = readtable('KO_data/timeseries_SKO_All.csv');

grps= unique(sko_df.Series)
% define the phenotype
T_SKO = cell2table(cell(0,3));
T_SKO.Properties.VariableNames = [{'Pathway','NB','Phenotype'}];


for i=1:numel(grps)
    grp = grps{i};
    %convert the 2d table to 1d list of entreez ids
    sko = sko_df(find(ismember(sko_df.Series,grp)),3:end);
    sko = string(reshape(table2cell(sko),[1,size(sko,1)*size(sko,2)]));
    sko = rmmissing(sko);
    sko = strjoin(sko,',');
    sko = unique(strsplit(sko,','));
    sko(find(ismember(sko,""))) = [];
    % Find shared genes between RECON model genes and DEGs
    genes_metabolic=find(ismember(recon_model.genes,sko));
    [r_genes,~]=find(recon_model.rxnGeneMat(:, genes_metabolic));
    genes_ss=recon_model.subSystems(r_genes);

    %
    [pathways_genes, ~, ub] = unique(string(genes_ss));
    path_genes = histc(ub, 1:length(pathways_genes))*1;
    grp = repmat(string(grp),size(path_genes,1),1);
    T_genes = table(pathways_genes, path_genes,grp);
    T_genes.Properties.VariableNames={'Pathway', 'NB','Phenotype'};
    T_SKO = [T_SKO;T_genes];
end

writetable(T_SKO,'KO_data/timeseries_RECON_SKO_Pathways.csv');


