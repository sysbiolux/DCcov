function [up_keep,down_keep] = doubleRxnDeletion(model,rxns_to_check)
%DoubleRxnsDeletion Apply double Reaction Deletion
% returns up and down regulated rxns
    changeCobraSolver('ibm_cplex');
    solverOK=changeCobraSolver('ibm_cplex','all')
    model_keep=model;
    up_keep=[];
    down_keep=[];
    for i=1:numel(rxns_to_check)
        %% reinitialize the model
        model= model_keep;
        model.ub(find(ismember(model.rxns, rxns_to_check(i))))=0;
        model.lb(find(ismember(model.rxns, rxns_to_check(i))))=0;
        [grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(model);
        down= model.rxns(grRatio<grRateWT*0.99);
        up= model.rxns(grRatio>grRateWT*1.01);
        down(:,2)=repmat(rxns_to_check(i), size(down,1),1);
        up(:,2)=repmat(rxns_to_check(i), size(up,1),1);
        up_keep=[up_keep;up];
        down_keep=[down_keep;up];
        p = i/numel(rxns_to_check)*100;
        fprintf('Double Reaction Deletion in progress ... %.4f \n', p)
    end
end