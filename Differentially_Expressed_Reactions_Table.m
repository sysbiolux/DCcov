%% Build a table of differentially expressed reaction in GSE147507

%Load Recon2 model
load('consistRecon2_4.mat')
%Removing gene version from the model
Recon2 = model
genes = Recon2.genes;
for i=1:numel(genes)
    x = strsplit(table2array(genes(i,1)),'.');
    genes(i,1) = cellstr(x(1,1));
end
Recon2.genes = genes;

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

% List DEGs files
degs_files = dir('./data/severity_study/DEGs/');
degs_files = {degs_files.name};
degs_files = degs_files(find(contains(degs_files,'.txt')));

T_DER = cell2table(cell(0,4));
T_DER.Properties.VariableNames = [{'Pathway','Down_regulated','Up_regulated','Phenotype'}];

unq_rxns = cell2table(cell(numel(degs_files),3));
unq_rxns.Properties.VariableNames = [{'Phenotype','Down_regulated','Up_regulated'}];

for i=1:numel(degs_files)
    %degs_df = readtable('BioJupies/DEGS_all_genes.csv');
    filepath  =string(strcat('./data/severity_study/DEGs/',  degs_files(i)));
    opts = detectImportOptions(filepath,'NumHeaderLines',0);
    degs_df = readtable(filepath,opts);
    colnames = degs_df.Properties.VariableNames;
    colnames = colnames(1:6);
    colnames = ['Genes',colnames];
    degs_df.Properties.VariableNames = colnames
    %degs_df.Properties.RowNames = degs_df.Var1;
    degs_df.Properties.RowNames = degs_df.Genes;
    %up = strsplit(string(degs_df{3,2}),'; ')' ;
    %down = strsplit(string(degs_df{3,3}),'; ')' ;
    up = degs_df.Genes(degs_df.padj<0.05 & degs_df.log2FoldChange>=1);
    down = degs_df.Genes(degs_df.padj<0.05 & degs_df.log2FoldChange<=-1);

    %mapping DEGs symbols to entrez
    up_idx = find(ismember(string(up),dico_RECON.SYMBOL))
    up_entrez = dico_RECON{up_idx,2};
    down_idx = find(ismember(string(down),dico_RECON.SYMBOL))
    down_entrez = dico_RECON{down_idx,2};
    % Find shared genes between RECON model genes and DEGs
    up_metabolic=find(ismember(recon_model.genes,up_entrez));
    down_metabolic=find(ismember(recon_model.genes,down_entrez));
    [r_down,~]=find(recon_model.rxnGeneMat(:, down_metabolic));
    [r_up,~]=find(recon_model.rxnGeneMat(:, up_metabolic));
    % Map unique rxns names for each phenotype to calculate the p-value
    unq_rxns(i,2) = {strjoin( unique(recon_model.rxnNames(r_down)),';')}
    unq_rxns(i,3) = {strjoin(unique(recon_model.rxnNames(r_up)),';')};
    unq_rxns(i,1) = {degs_files{i}(1:end-4)};
    
    both=intersect(r_up,r_down);
    up_ss=recon_model.subSystems(setdiff(r_up, both));
    down_ss=recon_model.subSystems(setdiff(r_down, both));
    [pathways_down, ~, ub] = unique(string(down_ss));
    path_down = histc(ub, 1:length(pathways_down))*-1;
    Tdown = table(pathways_down, path_down);
    %
    [pathways_up, ~, ub] = unique(string(up_ss));
    path_up = histc(ub, 1:length(pathways_up))*1;
    Tup = table(pathways_up, path_up);
    Tup.Properties.VariableNames={'Pathways', 'NB'};
    pathway_diff=union(pathways_up, pathways_down);
    match=find(ismember(pathway_diff, pathways_down));
    u=zeros(numel(pathway_diff),1);
    u(match)=path_down;
    match=find(ismember(pathway_diff, pathways_up));
    v=zeros(numel(pathway_diff),1);
    v(match)=path_up;
    grp = repmat(string(degs_files{i}(1:end-4)),size(pathway_diff,1),1);
    T_diff2=table(pathway_diff,u,v,grp);
    T_diff2.Properties.VariableNames =  [{'Pathway','Down_regulated','Up_regulated','Phenotype'}];
    T_DER = [T_DER;T_diff2];
end
writetable(T_DER,'./data/severity_study/Differentially_Expressed_Reactions.csv');
writetable(unq_rxns,'./data/severity_study/DEGs_Unique_Reactions.csv');

% Pivoting T_DER to 2 sparse arrays m*n (n: phenotypes, m: reactions)
uniq_phenotypes = unique(T_DER.Phenotype)';
uniq_rxns = unique(T_DER.Pathway)';
% ccreate a table of rxn pathways
T_pathway = table(recon_model.subSystems,recon_model.rxnNames)
writetable(T_pathway,'./RECON_Pathway_Reactions.csv');
