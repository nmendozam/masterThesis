close all force; clear variables; clc

%% Startup the COBRA Toolbox
initCobraToolbox(false) % Don't update the toolbox
changeCobraSolver ('gurobi', 'all', 1); % For large models
clc

%% Load data and model
addpath('functions')  
% Load omic data
abundance = readtable('data/sup_material_8.xlsx');
expresion = readtable('data/DESeq2_normalised_counts.xls');
% Get list of proteinECNumbers
proteinECNumbers_FileName = 'ProteinECNumbers.mat';
if isfile(proteinECNumbers_FileName)
    load(proteinECNumbers_FileName, 'ProteinECNumbers');
else
    ProteinECNumbers = getProteinCodes(abundance);
    save(proteinECNumbers_FileName, 'ProteinECNumbers');
end

% Load Recon3D model
load('models/Recon3D_301.mat')
% Recon3D stats
N_filled = find(~cellfun(@isempty,Recon3D.rxnECNumbers));
percent_of_filled = length(N_filled)/length(Recon3D.rxnECNumbers) * 100;
fprintf('%.02f%% of the reactions in Recon3D are associated with EC numbers\n', percent_of_filled)
% Standarize Recon3D rxnECNumber rules
Recon3D = harmonizeReconECnumbers(Recon3D);
Recon3D = makePrules(Recon3D);

%% Map omic data to Recon3D reactions
disp('Mapping omic data to Recon3D reactions ...')
profile on
abundanceTable = mapAbundance(Recon3D, abundance, ProteinECNumbers);
expresion = EnsemblToEntrez(expresion, 'data/mart_export.txt');
expressionTable = mapExpression(Recon3D, expresion);
profile viewer

%% Transcriptomic and proteomic data integration
disp('Integrating Transcriptomic and Proteomic data ...')
omicIntegratedData = omicIntegrationPCA(abundanceTable, expressionTable);

%% Recon3D Model reduction to astrocite specific GEM
disp('Creating Astrocite Specific Model ...')
options.solver = 'iMAT';
options.expressionRxns = omicIntegratedData;
options.threshold_lb = 1000;
options.threshold_ub = -1000;

astrociteModel = createTissueSpecificModel(Recon3D, options);
astrociteModel.description = 'MendozaAstrociteModel';
save('out/dirtyAstrociteModel.mat', 'astrociteModel');

%% First FBA to check flux
FBAsolution = optimizeCbModel(astrociteModel,'max');
fprintf("The first FBA flux after reduction is %f\n", FBAsolution.f)

%% Cleaning the model
modelClosed = modelHarmonization(astrociteModel);
[modelExchanges, selExc] = findExchangeReactions(modelClosed);
% Set lower bound of biomass reaction to 0 
modelClosed.lb(ismember(modelClosed.rxns,'biomass_reaction'))=0;
modelClosed.lb(ismember(modelClosed.rxns,'biomass_maintenance_noTrTr'))=0;
modelClosed.lb(ismember(modelClosed.rxns,'biomass_maintenance'))=0;
% Set exchange reactions
modelClosed.lb(ismember(modelClosed.rxns,modelClosed.rxns(modelExchanges)))=0;
modelClosed.ub(ismember(modelClosed.rxns,modelClosed.rxns(modelExchanges)))=1000;
% Remove leaking metabolites
disp('Searching and removing Leaking metabolites ...')
modelNoLeak = modelClosed;
[LeakRxns,modelTested,LeakRxnsFluxVector] = fastLeakTest(modelNoLeak, modelNoLeak.rxns(selExc), 'false');
if ~isempty(LeakRxns)
    disp('Removing leaking metabolites');
    modelNoLeak = removeRxns(modelNoLeak, LeakRxns);
end
save('out/NoLeakAstrociteModel.mat', 'modelNoLeak');

% FBA after removing the leaking metabolites
modelExchanges = findExchangeReactions(modelNoLeak);
modelNoLeak.lb(ismember(modelNoLeak.rxns,modelNoLeak.rxns(modelExchanges)))=-1000;
modelNoLeak.ub(ismember(modelNoLeak.rxns,modelNoLeak.rxns(modelExchanges)))=1000;
FBAsolution = optimizeCbModel(modelNoLeak,'max');
fprintf("The FBA flux after removing the leaking metabolites is %f\n", FBAsolution.f)

%% Set DMEM medium
astrociteModelDMEM = setBoundriesDMEMedium(modelNoLeak);
FBAsolution = optimizeCbModel(astrociteModelDMEM ,'max');
% How many exchange reactions are active
modelexchanges1 = strmatch('Ex_',astrociteModelDMEM.rxns);
modelexchanges4 = strmatch('EX_',astrociteModelDMEM.rxns);
modelExchanges = unique([modelexchanges1;modelexchanges4]);
exRxns = find(ismember(astrociteModelDMEM.rxns,astrociteModelDMEM.rxns(modelExchanges)) & FBAsolution.v ~= 0);
astrociteModelDMEM.rxns(exRxns)

