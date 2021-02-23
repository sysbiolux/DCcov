
% Calculate median, minimum and maximum for all models across the 2 studies

%Creating a table of model_path, #_mock_mets, #_infected_mets, #_mock_rxns, #_infected_rxns
T = cell2table(cell(28,5));

% List models mat files for the severity study  

model_files = dir('./output');
model_files = {model_files.name};
model_files = model_files(find(contains(model_files,'final_recon')));
% Remove S15 , lung biobsy & S16 with drug pertubation
model_files = model_files(~contains(model_files,'s15'))
model_files = model_files(~contains(model_files,'s16d'))
T(1:12,1)  = strcat('./output/',model_files)';

% List models mat files for the time-series study
model_files = dir('./output_');
model_files = {model_files.name};
model_files = model_files(find(startsWith(model_files,'SKO')));
T(13:28,1)  = strcat('./output_/',model_files)';

for i=1:size(T,1)
    filepath = string(T{i,1});
    X = load(filepath) ;
    models = fieldnames(X) ;
    %finds if the codition have ctl model based on rpkm size
    has_ctl_model  = models(find(startsWith(models,'rpkm_ctl')))    
    
    models = models(find(startsWith(models,'model_')));
    % determine the variable names of the models
    model_ctl = models(find(contains(models,'_ctl')));
    model_cov = models(find(contains(models,'_cov')));
    
    %count the number of rxns and mets
    rxns_cov = numel(X.(string(model_cov)).rxns);
    mets_cov = numel(X.(string(model_cov)).mets);
    
    if  size(has_ctl_model,1) ==0
        rxns_ctl = numel(X.(string(model_ctl)).rxns);
        mets_ctl = numel(X.(string(model_ctl)).mets);
    elseif numel(X.(string(has_ctl_model)))>=1
        rxns_ctl = numel(X.(string(model_ctl)).rxns);
        mets_ctl = numel(X.(string(model_ctl)).mets);        
    else 
        rxns_ctl = 0;
        mets_ctl = 0;
    end
    T(i,2:5) = table({mets_ctl},{mets_cov},{rxns_ctl},{rxns_cov});
end

T.Properties.VariableNames = [{'model_path'},{'mock_mets'},{'infected_mets'},{'mock_rxns'},{'infected_rxns'}];

% Number of mets without conditions with no healthy samples
x =  cell2mat([T{:,2};T{:,3}]);
x = x(x>1);
median(x)
min(x)
max(x)

% Number of rxns without conditions with no healthy samples
x =  cell2mat([T{:,4};T{:,5}]);
x = x(x>1);
median(x)
min(x)
max(x)