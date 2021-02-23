function [SKO_dn_safe,SKO_dn_toxic,SKO_up_safe,SKO_up_toxic] = Find_SKO_Rxns(model_inf,model_ctl)
%Find_SKO_Rxns Apply Single Rxn deletion on infected and healthy models
%   Returns safe and toxic rxns
    changeCobraSolver('ibm_cplex');
    solverOK=changeCobraSolver('ibm_cplex','all')

    threshold_down = 0.7;
    threshold_up = 1.3;

    [grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(model_inf);
    rxns_down =  model_inf.rxns(grRatio<= threshold_down);
    rxns_up =  model_inf.rxns(grRatio >= threshold_up);

    [grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(model_ctl);
    rxns_down_ =  model_ctl.rxns(grRatio<= threshold_down);
    rxns_up_ =  model_ctl.rxns(grRatio >= threshold_up);

    %% Find SKO that is safe to target on healthy models
    SKO_dn_safe = setdiff(rxns_down,rxns_down_);
    SKO_dn_toxic = intersect(rxns_down,rxns_down_);
    SKO_up_safe = setdiff(rxns_up,rxns_up_);
    SKO_up_toxic = intersect(rxns_up,rxns_up_);
end

