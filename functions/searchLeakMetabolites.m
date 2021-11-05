function LeakRxns = searchLeakMetabolites(model)

modelClosed = modelHarmonization(model);
[~, selExc] = findExchangeReactions(modelClosed);
modelExchanges = findExcRxns(modelClosed);
% Set lower bound of biomass reaction to 0 
modelClosed.lb(ismember(modelClosed.rxns,'biomass_reaction'))=0;
modelClosed.lb(ismember(modelClosed.rxns,'biomass_maintenance_noTrTr'))=0;
modelClosed.lb(ismember(modelClosed.rxns,'biomass_maintenance'))=0;
% Set exchange reactions
modelClosed.lb(ismember(modelClosed.rxns,modelClosed.rxns(modelExchanges)))=0;
modelClosed.ub(ismember(modelClosed.rxns,modelClosed.rxns(modelExchanges)))=1000;
% find leaking metabolites
[LeakRxns,modelTested,LeakRxnsFluxVector] = fastLeakTest(modelClosed, modelClosed.rxns(selExc), 'false');

end