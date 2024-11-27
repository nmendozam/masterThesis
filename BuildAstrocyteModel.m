% This file contains the code to build the astrocyte model
% using proteomic and transcriptomic data located in 'data/' folder

close all force; clear variables; clc

%% Startup the COBRA Toolbox
addpath('cobratoolbox')
addpath('functions')
initCobraToolbox(false) % Don't update the toolbox
changeCobraSolver ('gurobi', 'all', 1); % For large models
clc

%% Load integrated omic data
load('out/omicIntegratedData.mat') %loads 'omicIntegratedData

%% Open and standardize models
load('models/Recon3D_301.mat')

% Standardize Recon3D rxnECNumber rules
Recon3D = harmonizeReconECnumbers(Recon3D);
Recon3D = makePrules(Recon3D);
% load base model
Osorio_model=readCbModel('models/Astrocyte_Osorio2019.xml');
Osorio_model.rxns = strrep(Osorio_model.rxns,"(","[");
Osorio_model.rxns = strrep(Osorio_model.rxns,")","]");

%% Reinforce Exchanges from Osorio
% The exchange reactions from the Osorio model are reinforced, as they were
% manually curated. Hence, by adding a high value to them, we ensure that
% they will be included in the final model.
% OsorioExchanges = Osorio_model.rxns(strmatch('EX_', Osorio_model.rxns));
% ReconExchanges_from_Osorio = ismember(Recon3D.rxns, OsorioExchanges);
% omicIntegratedData(ReconExchanges_from_Osorio) = omicIntegratedData(ReconExchanges_from_Osorio) + 500;

%% Recon3D Model reduction to astrocite specific GEM
disp('Creating Astrocite Specific Model ...')
options.solver = 'iMAT';
options.expressionRxns = omicIntegratedData();
options.threshold_lb = prctile(omicIntegratedData,0.99);
options.threshold_ub = prctile(omicIntegratedData,0.01);

astrociteModel = createTissueSpecificModel(Recon3D, options);
astrociteModel.description = 'MendozaAstrociteModel';
astrociteModel.rxnECNumbers = {astrociteModel.rxnECNumbers{:}}';
save('out/dirtyAstrociteModel.mat', 'astrociteModel');

FBAsolution = optimizeCbModel(astrociteModel,'max');
fprintf("The first FBA flux after reduction is %f\n", FBAsolution.f)

%% Cleaning the model
disp('Searching and removing Leaking metabolites ...')
LeakRxns = searchLeakMetabolites(astrociteModel);
if ~isempty(LeakRxns)
    disp('Seraching for a consistent network');
    [~, ~, ~, ~, modelNoLeak] = findFluxConsistentSubset(astrociteModel);
end
save('out/NoLeakAstrociteModel.mat', 'modelNoLeak');

FBAsolution = optimizeCbModel(modelNoLeak,'max');
fprintf("The FBA flux after removing the leaking metabolites is %f\n", FBAsolution.f)

%% Contextualize the model for all the reactions
minimum = min(omicIntegratedData(~isinf(omicIntegratedData)));
omicIntegratedData(isinf(omicIntegratedData)) = minimum;
boundaries = omicIntegratedData - minimum;

% Assuming Recon3D.rxnNames and model.rxnNames are cell arrays of strings
[~, idx] = ismember(modelNoLeak.rxnNames, Recon3D.rxnNames);
reduced_boundaries = boundaries(idx);

modelNoLeak.lb(modelNoLeak.lb < 0) = -1;
modelNoLeak.ub(modelNoLeak.ub > 0) = 1;
modelNoLeak.lb =  modelNoLeak.lb .* reduced_boundaries;
modelNoLeak.ub = modelNoLeak.ub .* reduced_boundaries ;

FBAsol = optimizeCbModel(modelNoLeak, 'max'); FBAsol.f
fprintf("The FBA flux after reduction to astrocyte %f\n", FBAsol.f)

%% Set growth medium
astrociteModelDMEM = setBoundriesDMEMedium(modelNoLeak);
FBAsolution = optimizeCbModel(astrociteModelDMEM ,'max');
fprintf("The FBA flux after constrain to DMEM %f\n", FBAsolution.f)


%% Contextualize with heatlhy
% Contextualize with the healthy omic data
% to assess the model's behavior in a healthy context before gap-filling
load('out/omicHealthy.mat');
minimum = min(omicHealthy(~isinf(omicHealthy)));
omicHealthy(isinf(omicHealthy)) = minimum;
boundaries = omicHealthy - minimum;

% Assuming Recon3D.rxnNames and model.rxnNames are cell arrays of strings
[~, idx] = ismember(astrociteModelDMEM.rxnNames, Recon3D.rxnNames);
reduced_boundaries = omicHealthy(idx);

healthyAstrocyte = astrociteModelDMEM;

healthyAstrocyte.lb(healthyAstrocyte.lb < 0) = -1;
healthyAstrocyte.ub(healthyAstrocyte.ub > 0) = 1;
healthyAstrocyte.lb =  healthyAstrocyte.lb .* reduced_boundaries;
healthyAstrocyte.ub = healthyAstrocyte.ub .* reduced_boundaries ;


model = healthyAstrocyte;
%% Remove cytidine dependency
model = addReaction(model, 'DCMPDA', 'reactionFormula', 'h2o[c] + h[c] + dcmp[c] <=> nh4[c] + dump[c]');
%% Addition of tibolone reactions
model = addReaction(model, 'EX_Tibolone[e]', 'reactionFormula', 'Tibolone[e] <=>');
model = addReaction(model, 'T2', 'reactionFormula', 'Tibolone[e] <=> a3OHtibolone[e]');
model = addReaction(model, 'T3', 'reactionFormula', 'Tibolone[e] <=> b3OHtibolone[e]');
model = addReaction(model, 'T4', 'reactionFormula', 'Tibolone[e] -> d4tibolone[e]');
model = addReaction(model, 'T5', 'reactionFormula', 'b3OHtibolone[e] -> d4tibolone[e]');
model = addReaction(model, 'T6', 'reactionFormula', 'a3OHtibolone[e] -> estradiol[c]');
model = addReaction(model, 'T7', 'reactionFormula', 'b3OHtibolone[e] -> estradiol[c]');
model = addReaction(model, 'T8', 'reactionFormula', 'd4tibolone[e] -> prgstrn[c] + tststerone[c]');
model = addReaction(model, 'T9', 'reactionFormula', 'a3OHtibolone[e] <=> a3SOtibolone[e]');
model = addReaction(model, 'EX_a3SOtibolone[e]', 'reactionFormula', 'a3SOtibolone[e] ->');
%% Addition of Objective functions
model = addReaction(model, 'Glu2Gln', 'reactionFormula', '1 glu_L[e] + 1 gln_L[c] -> 1 gln_L[e] + 1 glu_L[c]');
model = addReaction(model, 'Gly2SerD', 'reactionFormula', '1 gly[e] + 1 ser_D[c] -> 1 ser_D[e] + 1 gly[c]');
model = addReaction(model, 'Glc2Lac', 'reactionFormula', '1 glc_D[e] + 2 lac_L[c] -> 2 lac_L[e] + 1 glc_D[c]');
model = addReaction(model, 'Glc2ATP', 'reactionFormula', '1 glc_D[e] + 36 atp[c] -> 36 atp[e] + 1 glc_D[c]');
model = addReaction(model, 'Cys2GTHRD', 'reactionFormula', '1 cys_L[e] + 1 glu_L[c] + 1 gly[c] -> 1 gthrd[e]');
model = changeObjective(model, {'biomass_reaction', 'Glu2Gln', 'Gly2SerD', 'Glc2Lac', 'Glc2ATP', 'Cys2GTHRD'});

FBAsol = optimizeCbModel(model, 'max'); FBAsol.f
fprintf("The FBA flux after healthy context %f\n", FBAsol.f)