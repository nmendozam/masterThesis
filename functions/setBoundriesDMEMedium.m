function model = setBoundriesDMEMedium(model)

modelexchanges1 = strmatch('Ex_',model.rxns);
modelexchanges4 = strmatch('EX_',model.rxns);
modelExchanges = unique([modelexchanges1;modelexchanges4]);
% modelExchanges = findExchangeReactions(model);
model.lb(ismember(model.rxns,model.rxns(modelExchanges)))=0;
model.ub(ismember(model.rxns,model.rxns(modelExchanges)))=0;

% Open just the DMEM's exchange reactions
model = changeRxnBounds(model,'EX_gly[e]', -0.030,'l');
model = changeRxnBounds(model,'EX_gln_L[e]', -0.584,'l');
model = changeRxnBounds(model,'EX_ile_L[e]', -0.105,'l');
model = changeRxnBounds(model,'EX_leu_L[e]', -0.105,'l');
model = changeRxnBounds(model,'EX_met_L[e]', -0.030,'l');
model = changeRxnBounds(model,'EX_phe_L[e]', -0.066,'l');
model = changeRxnBounds(model,'EX_ser_L[e]', -0.042,'l');
model = changeRxnBounds(model,'EX_thr_L[e]', -0.095,'l');
model = changeRxnBounds(model,'EX_trp_L[e]', -0.016,'l');
model = changeRxnBounds(model,'EX_val_L[e]', -0.094,'l');
% this is a derivate of Folic acid which is the original metabolite in the DMEM
model = changeRxnBounds(model,'EX_dhf[e]', -0.004,'l'); 
model = changeRxnBounds(model,'EX_nad[e]', -0.004,'l');
model = changeRxnBounds(model,'EX_ribfiv[e]', -0.0004,'l');
model = changeRxnBounds(model,'EX_inost[e]', -0.0072,'l');
model = changeRxnBounds(model,'EX_glc_D[e]', -4.5,'l');

model = changeRxnBounds(model,'EX_udpglcur[e]', -0.05, 'l');
model = changeRxnBounds(model,'EX_M02047[e]', -0.05, 'l');
model = changeRxnBounds(model,'EX_cytd[e]', -0.016, 'l');
% model = changeRxnBounds(model,'EX_dcyt[e]', -0.005, 'l');

end