
% Calculate median, minimum and maximum for all models across the 2 studies

% List models mat files for the severity study  

model_files = dir('./models/severity/');
model_files = {model_files.name};
model_files2 = dir('./models/timeseries/');
model_files2 = {model_files2.name};
model_files = [model_files';model_files2']';
model_files = model_files(find(contains(model_files,'.mat')));

%count the number of rxns, mets and genes
%Creating a table of model, model_path, #mets, #rxns, #genes,
%Study,Condition, State, Generic_Model
T = cell(numel(model_files),9);
T(:,1) = model_files;
severity_conds= {'A549_0.02', 'A549_0.2_ACE2', 'NHBE_2', 'A549_2','A549_2_ACE2', 'Calu3_2'};
severity_conds_={'S2','S6','S1','S5','S16','S7'};
for i=1:size(T,1)
    %Study
    if i<=24
        dir = "./models/severity/";
        T(i,6) = {"Severity"};
    else
        dir = './models/timeseries/';
        T(i,6) = {"Timeseries"};
    end
    T(i,2) = {dir+string(model_files{i})};
    % Condition 
    grp = strsplit(string(T(i,1)),'_');

    % State (Infected / Mock)
    if contains(T(i,1), 'SARS')
        T(i,8) = {"Infected"};
    else
        T(i,8) = {"Mock"};
    end
    
    % Generic_Model 
    generic_model = strsplit(grp(end),'.');
    T(i,9) = {generic_model(1)};
    
    % Condition 
    grp = strsplit(string(T(i,1)),'_');
    if contains(string(T(i,2)), 'severity')
        cond =grp(1);
        cond_ = severity_conds(find(ismember(severity_conds_,cond)));
        T(i,7) = {cond_}; % change severity condition name
        T(i,1) = {replace(string(T(i,1)),cond,cond_)} ;% change severity model name
    else
        T(i,7) = {strjoin(grp(1:2),'_')};
    end
    T(i,7) = {replace(string(T(i,7)),'SARS','')}; % remove SARS from some condition names
       
    %

    X = load(string(T(i,2)));

    var_name = fieldnames(X);
    num_rxns = numel(unique(X.(string(var_name)).rxns));
    num_mets = numel(unique(X.(string(var_name)).mets));
    num_genes = numel(unique(X.(string(var_name)).genes));
    T(i,3) = {num_rxns};
    T(i,4) = {num_mets};
    T(i,5) = {num_genes};
end

T = cell2table(T);

T.Properties.VariableNames =[{'Model', 'model_path', 'Metabolites', 'Reactions', 'Genes','Study','Condition', 'State', 'Generic_Model'}];
writetable(T,'KO_data/Model_statistics.csv')

% Median number of metabolites
[median(T.Metabolites),min(T.Metabolites),max(T.Metabolites)]
% Median number of reactions
[median(T.Reactions),min(T.Reactions),max(T.Reactions)]
% Median number of genes
[median(T.Genes),min(T.Genes),max(T.Genes)]

%% Jaccard similirity of built models

% Recon3 models
T_recon3d= T(T.Generic_Model=='3D',:);
load('Recon3DModel_301.mat')
model_keep = zeros(numel(T_recon3d.Condition),numel(Recon3DModel.rxns));
for i=1:size(T_recon3d,1)
    X = load(T_recon3d.model_path(i));
    var_name = fieldnames(X);
    model_rxns = unique(X.(string(var_name)).rxns);
    model_keep(i,:) = ismember(Recon3DModel.rxns,model_rxns);
end

colnames = T_recon3d.Condition +'_'+T_recon3d.State;
colnames = replace(colnames,'_','\_');

J = squareform(pdist(model_keep,'jaccard'));

altcolor= [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
    255 0 0; 204 0 0; 152 0 0; 102 0 0; 51 0 0]/255;

%get and save the default size
defaultPosition = get(0,'DefaultFigurePosition');
%get the current screen size
screensize = get( groot, 'Screensize' );
%screensize = get(0, 'Screensize'); %for earlier Matlab versions (e.g. Matlab 2010)
%set default figure position to full screen
screensize = [1,1,800,800];
set(0, 'DefaultFigurePosition', screensize);

cgo_J = clustergram(1-J,...
    'RowLabels', colnames,...
    'ColumnLabels', colnames,...
    'ColumnLabelsRotate',315, ...
    'Cluster', 'all', ...
    'symmetric','False',...
    'Colormap', altcolor)
addTitle(cgo_J,{"Model similarity using Jaccard distance"," based on the Recon3D reconstructed models' reactions"})

%adding colorbar
cbButton = findall(gcf,'tag','HMInsertColorbar');
ccb = get(cbButton,'ClickedCallback')
set(cbButton,'State','on')
ccb{1}(cbButton,[],ccb{2})
%save plot
saveas(gcf,['Figs\Models_Jaccard_Similarity_Recon3D.png']);


% Recon2 models
T_recon2= T(T.Generic_Model=='2',:);
load('consistRecon2_4.mat')
model_keep = zeros(numel(T_recon2.Condition),numel(model.rxns));
for i=1:size(T_recon2,1)
    X = load(T_recon2.model_path(i));
    var_name = fieldnames(X);
    model_rxns = unique(X.(string(var_name)).rxns);
    model_keep(i,:) = ismember(model.rxns,model_rxns);
end

colnames = T_recon2.Condition +'_'+T_recon2.State;
colnames = replace(colnames,'_','\_');

J = squareform(pdist(model_keep,'jaccard'));

altcolor= [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
    255 0 0; 204 0 0; 152 0 0; 102 0 0; 51 0 0]/255;

%get and save the default size
defaultPosition = get(0,'DefaultFigurePosition');
%get the current screen size
screensize = get( groot, 'Screensize' );
%screensize = get(0, 'Screensize'); %for earlier Matlab versions (e.g. Matlab 2010)
%set default figure position to full screen
screensize = [1,1,800,800];
set(0, 'DefaultFigurePosition', screensize);

cgo_J = clustergram(1-J,...
    'RowLabels', colnames,...
    'ColumnLabels', colnames,...
    'ColumnLabelsRotate',315, ...
    'Cluster', 'all', ...
    'symmetric','False',...
    'Colormap', altcolor)
addTitle(cgo_J,{"Model similarity using Jaccard distance"," based on the Recon2 reconstructed models' reactions"})

%adding colorbar
cbButton = findall(gcf,'tag','HMInsertColorbar');
ccb = get(cbButton,'ClickedCallback')
set(cbButton,'State','on')
ccb{1}(cbButton,[],ccb{2})

saveas(gcf,['Figs\Models_Jaccard_Similarity_Recon2.png']);


% cgfig = findall(0,'type','figure', 'tag', 'Clustergram'); % Figure handle
% cgax = findobj(cgfig, 'type','axes','tag','HeatMapAxes'); % main axis handle
% % Change fontsize
% cgax.Fontsize = 8;

