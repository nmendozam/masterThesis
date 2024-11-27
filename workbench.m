close all force; clear variables; clc
%% Startup the COBRA Toolbox
addpath('cobratoolbox')
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

Osorio_model=readCbModel('models/Astrocyte_Osorio2019.xml');
Osorio_model.rxns = strrep(Osorio_model.rxns,"(","[");
Osorio_model.rxns = strrep(Osorio_model.rxns,")","]");

Martin_model=readCbModel('models/Astrocyte_Martin-Jimenez2017.xml');
Martin_model.rxns = strrep(Martin_model.rxns,"(","[");
Martin_model.rxns = strrep(Martin_model.rxns,")","]");

% Load Recon3D model
load('models/Recon3D_301.mat')
% Recon3D stats
N_filled = find(~cellfun(@isempty,Recon3D.rxnECNumbers));
percent_of_filled = length(N_filled)/length(Recon3D.rxnECNumbers) * 100;
fprintf('%.02f%% of the reactions in Recon3D are associated with EC numbers\n', percent_of_filled)
N_filled = find(~cellfun(@isempty,Recon3D.rules));
percent_of_filled = length(N_filled)/length(Recon3D.rules) * 100;
fprintf('%.02f%% of the reactions in Recon3D are associated with gene rules\n', percent_of_filled)
% Standarize Recon3D rxnECNumber rules
Recon3D = harmonizeReconECnumbers(Recon3D);
Recon3D = makePrules(Recon3D);

%% Map omic data to Recon3D reactions
disp('Mapping omic data to Recon3D reactions ...')

abundanceTable_FileName = 'out/abundanceTable.mat';
expressionTable_FileName = 'out/expressionTable.mat';
if isfile(expressionTable_FileName) && isfile(abundanceTable_FileName)
    load(abundanceTable_FileName, 'abundanceTable');
    load(expressionTable_FileName, 'expressionTable');
else
    profile on
    abundanceTable = mapAbundance(Recon3D, abundance, ProteinECNumbers);
    expresion = EnsemblToEntrez(expresion, 'data/mart_export.txt');
    expressionTable = mapExpression(Recon3D, expresion);
    profile viewer

    save(abundanceTable_FileName, 'abundanceTable');
    save(expressionTable_FileName, 'expressionTable');
end


%% Transcriptomic and proteomic data integration
disp('Integrating Transcriptomic and Proteomic data ...')
omicIntegratedData = omicIntegrationPCA(abundanceTable, expressionTable);
omicIntegratedData = log(abs(min(omicIntegratedData)) + omicIntegratedData);

%% Reinforce Exchanges from Osorio
OsorioExchanges = Osorio_model.rxns(strmatch('EX_', Osorio_model.rxns));
ReconExchanges_from_Osorio = ismember(Recon3D.rxns, OsorioExchanges);
omicIntegratedData(ReconExchanges_from_Osorio) = omicIntegratedData(ReconExchanges_from_Osorio) + 500;

%% Recon3D Model reduction to astrocite specific GEM
disp('Creating Astrocite Specific Model ...')
options.solver = 'iMAT';
options.expressionRxns = omicIntegratedData();
options.threshold_lb = 500;
options.threshold_ub = -500;

astrociteModel = createTissueSpecificModel(Recon3D, options);
astrociteModel.description = 'MendozaAstrociteModel';
astrociteModel.rxnECNumbers = {astrociteModel.rxnECNumbers{:}}';
save('out/dirtyAstrociteModel.mat', 'astrociteModel');

%% First FBA to check flux
FBAsolution = optimizeCbModel(astrociteModel,'max');
fprintf("The first FBA flux after reduction is %f\n", FBAsolution.f)

%% Cleaning the model
disp('Searching and removing Leaking metabolites ...')
LeakRxns = searchLeakMetabolites(astrociteModel);
if ~isempty(LeakRxns)
    disp('Seraching for a consistent network');
    [~, ~, ~, ~, modelNoLeak] = findFluxConsistentSubset(astrociteModel);
%     modelNoLeak = removeRxns(astrociteModel, LeakRxns);
end
save('out/NoLeakAstrociteModel.mat', 'modelNoLeak');

% FBA after removing the leaking metabolites
FBAsolution = optimizeCbModel(modelNoLeak,'max');
fprintf("The FBA flux after removing the leaking metabolites is %f\n", FBAsolution.f)

%% Set DMEM medium
astrociteModelDMEM = setBoundriesDMEMedium(modelNoLeak);
FBAsolution = optimizeCbModel(astrociteModelDMEM ,'max');
fprintf("The FBA flux after constrain to DMEM %f\n", FBAsolution.f)

Exchanges = {'EX_gly[e]', 'EX_gln_L[e]', 'EX_ile_L[e]', 'EX_leu_L[e]', 'EX_met_L[e]', 'EX_phe_L[e]', 'EX_ser_L[e]', 'EX_thr_L[e]', 'EX_trp_L[e]', 'EX_val_L[e]', 'EX_dhf[e]', 'EX_nad[e]', 'EX_ribfiv[e]', 'EX_inost[e]', 'EX_glc_D[e]', 'EX_udpglcur[e]', 'EX_cytd[e]'};
result = {};
result(:,1) = FBAsolution.v(ismember(astrociteModelDMEM.rxns, Exchanges))
result(1,:)=astrociteModelDMEM.rxns(ismember(astrociteModelDMEM.rxns, Exchanges))

FBAsolution = optimizeCbModel(Osorio_model ,'max');
fprintf("The FBA flux for osorio model %f\n", FBAsolution.f)
Osorio_modelDMEM = setBoundriesDMEMedium(Osorio_model);
FBAsolution = optimizeCbModel(Osorio_modelDMEM ,'max');
fprintf("The FBA flux after constrain to DMEM %f\n", FBAsolution.f)


FBAsolution.v(strmatch('EX_', Osorio_model.rxns))
Osorio_model.rxns(strmatch('EX_', Osorio_model.rxns))

only_recon3d = compareModelsRxns(Recon3D, astrociteModelDMEM);
only_astrcyte = compareModelsRxns(astrociteModelDMEM, Recon3D);
fprintf("reactions only recon: %d only astrocyte %d shared %d \n", ...
    length(only_recon3d), length(only_astrcyte), length(Recon3D.rxns) - length(only_recon3d))

only_osorio = compareModelsRxns(Osorio_model, astrociteModelDMEM);
only_astrcyte = compareModelsRxns(astrociteModelDMEM, Osorio_model);
fprintf("reactions only osorio: %d only astrocyte %d shared %d %d\n", ...
    length(only_osorio), length(only_astrcyte), length(Osorio_model.rxns) - length(only_osorio), length(astrociteModelDMEM.rxns) - length(only_astrcyte))
%% Contexto
omicIntegratedData(omicIntegratedData < 0) = 0;

%% Graficas
figure;
subplot(2,3,1);
histogram(table2array(abundanceTable(:,:)))
ylim([0 14000])
title('Abundances in Rxns')
subplot(2,3,2);
histogram(table2array(expressionTable(:,:)))
ylim([0 14000])
title('Expressions in Rxns')
subplot(2,3,3);
histogram(omicIntegratedData)
ylim([0 14000])
title('Integrated data in Rxns')

subplot(2,3,4);
histogram(real(log(table2array(abundanceTable(:,:)))))
ylim([0 8000])
title('Log Abundances in Rxns')
subplot(2,3,5);
histogram(real(log(table2array(expressionTable(:,:)))))
ylim([0 8000])
title('Log Expressions in Rxns')
subplot(2,3,6);
histogram(real(log(omicIntegratedData)))
ylim([0 8000])
title('Log Integrated data in Rxns')

Exchange = modelNoLeak.rxns(strncmp('EX_', modelNoLeak.rxns, 3));
Exchange(contains(Exchange, 'ac'))

% printRxnFormula(modelNoLeak, 'rxnAbbrList', Exchange,'printFlag', true, 'metNameFlag', true)




astrociteModelDMEM = setBoundriesDMEMedium(astrociteModel);
FBAsolution = optimizeCbModel(astrociteModelDMEM ,'max');
fprintf("The FBA flux after constrain to DMEM %f\n", FBAsolution.f)

Exchange = astrociteModelDMEM.rxns(strncmp('EX_', astrociteModelDMEM.rxns, 3));
Exchange = Exchange(~strcmp(Exchange,'EX_dcmp[e]')); % desoxicitidina monofosfato  ADN
Exchange = Exchange(~strcmp(Exchange,'EX_cytd[e]')); % citidina
Exchange = Exchange(~strcmp(Exchange,'EX_dcyt[e]')); % Deoxycytidine
for i=1:length(Exchange)
    Current = Exchange{i};
    astrociteModelDMEM = changeRxnBounds(astrociteModelDMEM, Current, -1000, 'l');
    
    FBAsolution = optimizeCbModel(astrociteModelDMEM ,'max');
    if FBAsolution.f > 0 || i == 806
        fprintf("%s got %f\n", Current, FBAsolution.f)
    else
        for j=1:length(Exchange)
            Current_ = Exchange{j};
            astrociteModelDMEM = changeRxnBounds(astrociteModelDMEM, Current_, -1000, 'l');

            FBAsolution_ = optimizeCbModel(astrociteModelDMEM ,'max');
            if FBAsolution_.f > 0
                fprintf("%s and %s got %f\n", Current, Current_, FBAsolution_.f)
            end
            astrociteModelDMEM = changeRxnBounds(astrociteModelDMEM, Current_, 0, 'l');
        end
    end
    
    astrociteModelDMEM = changeRxnBounds(astrociteModelDMEM, Current, 0, 'l');
end

% Exchange = strncmp('EX_', astrociteModelDMEM.rxns, 3);
% astrociteModelDMEM.rxns(Exchange & FBAsolution.v < 0)
% FBAsolution.v(Exchange & FBAsolution.v < 0)

%% See esential reacions
allModels = {};
newModelName = 'AstrociteDMEM';
allModels.(newModelName) = astrociteModelDMEM;
newModelName = 'Osorio';
allModels.(newModelName) = Osorio_model;
[essentialRxn4Models, dataStruct] = essentialRxn4MultipleModels(allModels , 'biomass_reaction');

x = table2cell(essentialRxn4Models);
names = x(:,1);
essentials = names([x{:,2}] < 0.05);
essentials(strncmp('EX_', essentials, 3))

