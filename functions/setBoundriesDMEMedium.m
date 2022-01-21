function model = setBoundriesDMEMedium(model)

modelexchanges1 = strmatch('Ex_',model.rxns);
modelexchanges4 = strmatch('EX_',model.rxns);
modelExchanges = unique([modelexchanges1;modelexchanges4]);
% modelExchanges = findExchangeReactions(model);
% model.lb(ismember(model.rxns,model.rxns(modelExchanges)))=-1000;
model.lb(ismember(model.rxns,model.rxns(modelExchanges)))=0;
% model.ub(ismember(model.rxns,model.rxns(modelExchanges)))=0;


% Open just the DMEM's exchange reactions
model = changeRxnBounds(model,'EX_gly[e]', -0.030,'l');
model = changeRxnBounds(model,'EX_gln_L[e]', -0.584,'l');
model = changeRxnBounds(model,'EX_ile_L[e]', -0.105,'l');
model = changeRxnBounds(model,'EX_leu_L[e]', -0.105,'l');
model = changeRxnBounds(model,'EX_met_L[e]', -0.030,'l');
model = changeRxnBounds(model,'EX_phe_L[e]', -0.066,'l');
model = changeRxnBounds(model,'EX_ser_L[e]', -0.042,'l');
model = changeRxnBounds(model,'EX_thr_L[e]', -0.095,'l');
% model = changeRxnBounds(model,'EX_trp_L[e]', -0.016,'l');
model = changeRxnBounds(model,'EX_val_L[e]', -0.094,'l');
% this is a derivate of Folic acid which is the original metabolite in the DMEM
model = changeRxnBounds(model,'EX_dhf[e]', -0.004,'l'); 
model = changeRxnBounds(model,'EX_nad[e]', -0.004,'l');
% model = changeRxnBounds(model,'EX_ribfiv[e]', -0.0004,'l');
model = changeRxnBounds(model,'EX_inost[e]', -0.0072,'l');
model = changeRxnBounds(model,'EX_glc_D[e]', -4.5,'l');

% Open the FBS' exchange realctions
model = changeRxnBounds(model,'EX_inost[e]', -69.9 * 10^-3, 'l');
model = changeRxnBounds(model,'EX_thymd[e]', -1.5 * 10^-3, 'l');
model = changeRxnBounds(model,'EX_ribflv[e]', -581.9 * 10^-6, 'l');
model = changeRxnBounds(model,'EX_pydx[e]', -150.7 * 10^-6, 'l');
model = changeRxnBounds(model,'EX_cys_L[e]', -100.2 * 10^-3, 'l');
model = changeRxnBounds(model,'EX_asn_L[e]', -56.8 * 10^-3, 'l');
model = changeRxnBounds(model,'EX_trp_L[e]', -44.2 * 10^-3, 'l');
model = changeRxnBounds(model,'EX_ser_L[e]', -249.8 * 10^-3,'l');
model = changeRxnBounds(model,'EX_ala_L[e]', -50 * 10^-3,'l');
model = changeRxnBounds(model,'EX_arg_L[e]', -846.7 * 10^-3,'l');
model = changeRxnBounds(model,'EX_val_L[e]', -451.3 * 10^-3,'l');
model = changeRxnBounds(model,'EX_ile_L[e]', -415.2 * 10^-3,'l');
model = changeRxnBounds(model,'EX_leu_L[e]', -450.1 * 10^-3,'l');
model = changeRxnBounds(model,'EX_lys_L[e]', -624.1 * 10^-3,'l');
model = changeRxnBounds(model,'EX_pro_L[e]', -149.9 * 10^-3,'l');
model = changeRxnBounds(model,'EX_thr_L[e]', -448.8 * 10^-3,'l');
model = changeRxnBounds(model,'EX_Lcystin[e]', -130.2 * 10^-3,'l');
model = changeRxnBounds(model,'EX_gln_L[e]', -2.5,'l');
model = changeRxnBounds(model,'EX_hxan[e]', -15.4 * 10^-3,'l');
model = changeRxnBounds(model,'EX_lipoate[e]', -508.9 * 10^-3,'l');

%% FROM: https://www.researchgate.net/figure/Lipid-composition-of-FBS-and-FBSLess-mM-sera_fig13_49820996
% FBS' Triglycerides
model = changeRxnBounds(model, 'EX_tag_hs[e]', -1.26, 'l');
% FBS' phospholipids
% model = changeRxnBounds(model, 'EX_pe_hs[e]', -0.72, 'l');
model = changeRxnBounds(model, 'EX_ps_hs[e]', -0.72, 'l');
% FBS' cholesterol
model = changeRxnBounds(model, 'EX_chsterol[e]', -1.1, 'l'); % HDL


% model = changeRxnBounds(model, 'EX_56dura[e]', 0, 'l');     % 5,6-Dihydrouracil
% model = changeRxnBounds(model, 'EX_cbasp[e]', 0, 'l');      % N-Carbamoyl-L-Aspartate
% model = changeRxnBounds(model, 'EX_csn[e]', 0, 'l');        % Cytosine
% model = changeRxnBounds(model, 'EX_duri[e]', 0, 'l');       % Deoxyuridine
% model = changeRxnBounds(model, 'EX_hmcr[e]', 0, 'l');       % Homocitrulline
% model = changeRxnBounds(model, 'EX_orot[e]', 0, 'l');   % Orotidine
% model = changeRxnBounds(model, 'EX_orot5p[e]', 0, 'l');     % Orotate5'-Phosphate
% model = changeRxnBounds(model, 'EX_udpglcur[e]', 0, 'l');   % Uridine-5'-Diphosphate-Glucuronate
% model = changeRxnBounds(model, 'EX_ura[e]', 0, 'l');        % Uracil
% model = changeRxnBounds(model, 'EX_utp[e]', 0, 'l');        % Uridine-5'-TrIphosphate
% model = changeRxnBounds(model, 'EX_utp[e]', -0.1, 'l');        % Uridine-5'-TrIphosphate

%% Exchange gases
% O2, CO2, H2, NH4 and H2S b
model = changeRxnBounds(model,'EX_o2[e]', -1000,'l');
model = changeRxnBounds(model,'EX_co2[e]', -1000,'l');
% model = changeRxnBounds(model,'EX_h2[e]', -1000,'l');
model = changeRxnBounds(model,'EX_nh4[e]', -1000,'l');
% model = changeRxnBounds(model,'EX_h2s[e]', -1000,'l');
model = changeRxnBounds(model,'EX_no[e]', -1000,'l');
%% Water
model = changeRxnBounds(model,'EX_h2o[e]', -1000,'l');

%% SOURCE: https://www.researchgate.net/publication/15827024_Purine_base_and_nucleoside_cytidine_and_uridine_concentrations_in_foetal_calf_and_other_sera
model = changeRxnBounds(model,'EX_cytd[e]', -1.3 * 10^-3, 'l');        % Cytidine
model = changeRxnBounds(model, 'EX_uri[e]', -5.1 * 10^-3, 'l');        % Uridine
model = changeRxnBounds(model, 'EX_urate[e]', -130 * 10^-3, 'l');      % Urate
model = changeRxnBounds(model, 'EX_hxan[e]', -74.7 * 10^-3, 'l');      % Hypoxanthine


model = changeRxnBounds(model, 'EX_hmcr[e]', -30.1 * 10^-3, 'l');        % Homocitrulline


end