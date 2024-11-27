
close all force; clear variables; clc
%% Startup the COBRA Toolbox
addpath('cobratoolbox')
initCobraToolbox(false) % Don't update the toolbox
changeCobraSolver('gurobi', 'all', 1); % For large models
clc

%% Load data and base model
% Load transcriptomic and proteomic data
addpath('functions')  
abundance = readtable('data/sup_material_8.xlsx');
expresion = readtable('data/DESeq2_normalised_counts.xls');
disp("Omic data is loaded successfully!")

% Load Recon3D model
load('models/Recon3D_301.mat')
disp("Recon3D model is loaded successfully!")

%% Clean and standardize Recon3D model

% Standarize Recon3D rxnECNumber rules
Recon3D = harmonizeReconECnumbers(Recon3D);
Recon3D = makePrules(Recon3D);

%% Map omic data to Recon3D reactions
% Get list of proteinECNumbers
proteinECNumbers_FileName = 'ProteinECNumbers.mat';
if isfile(proteinECNumbers_FileName)
    load(proteinECNumbers_FileName, 'ProteinECNumbers');
else
    ProteinECNumbers = getProteinCodes(abundance);
    save(proteinECNumbers_FileName, 'ProteinECNumbers');
end

abundanceTable_FileName = 'out/abundanceTable.mat';
expressionTable_FileName = 'out/expressionTable.mat';
if isfile(expressionTable_FileName) && isfile(abundanceTable_FileName)
    load(abundanceTable_FileName, 'abundanceTable');
    load(expressionTable_FileName, 'expressionTable');
else
    fprintf("The files %s and %s were not found. Mapping omic data to Recon3D reactions...\n", abundanceTable_FileName, expressionTable_FileName);
    disp("This may take a while...")
    % profile on
    abundanceTable = mapAbundance(Recon3D, abundance, ProteinECNumbers);
    expresion = EnsemblToEntrez(expresion, 'data/mart_export.txt');
    expressionTable = mapExpression(Recon3D, expresion);
    % profile viewer

    save(abundanceTable_FileName, 'abundanceTable');
    save(expressionTable_FileName, 'expressionTable');
end


%% Integrate omic data into a single dimension.
% This is done with PCA and will be used to determine the reactions
% that are active in the astrocyte and need to be included in the
% astrocyte specific model.
disp('Integrating Transcriptomic and Proteomic data ...')
omicIntegratedData = omicIntegrationPCA(abundanceTable, expressionTable);
omicIntegratedData = log(abs(min(omicIntegratedData)) + omicIntegratedData);

save('out/omicIntegratedData.mat', 'omicIntegratedData');

%% Context-specific omic data
heal_samples = ["Control_NHA1_veh_tech2", "Control_NHA1_veh_tech1", "Control_NHA2_veh_tech1", "Control_NHA2_veh_tech2", "Control_NHA3_veh_tech1", "Control_NHA3_veh_tech2"];
pal_samples = ["Sample_NHA1_pal_tech2", "Sample_NHA1_pal_tech1", "Sample_NHA2_pal_tech1", "Sample_NHA2_pal_tech2", "Sample_NHA3_pal_tech1", "Sample_NHA3_pal_tech2"];
tib_samples = ["Sample_NHA1_tib_pal_tech2", "Sample_NHA1_tib_pal_tech1", "Sample_NHA2_tib_pal_tech1", "Sample_NHA2_tib_pal_tech2", "Sample_NHA3_tib_pal_tech1", "Sample_NHA3_tib_pal_tech2"];

heal_samples_exp = ["NHA1_DMEM", "NHA1_DMEM_1", "NHA2_DMEM", "NHA1_DMEM_1", "NHA3_DMEM", "NHA3_DMEM_1"];
pal_samples_exp = ["NHA1_palmitato", "NHA1_palmitato_1", "NHA2_palmitato", "NHA2_palmitato_1", "NHA3_palmitato", "NHA3_palmitato_1"];
tib_samples_exp = ["NHA3_tib_pal", "NHA3_tib_pal_1", "NHA2_tib_pal", "NHA2_tib_pal_1", "NHA1_tib_pal", "NHA1_tib_pal_1"];


% Heathy astrocites
heal_abundance = abundanceTable(:, heal_samples);
heal_expression = expressionTable(:, heal_samples_exp);
omicHealthy = omicIntegrationPCA(heal_abundance, heal_expression);
omicHealthy = log(abs(min(omicHealthy)) + omicHealthy);
save('out/omicHealthy.mat', 'omicHealthy');
% Inflamed with palmitate
pal_abundance = abundanceTable(:, pal_samples);
pal_expression = expressionTable(:, pal_samples_exp);
omicInflamed = omicIntegrationPCA(pal_abundance, pal_expression);
omicInflamed = log(abs(min(omicInflamed)) + omicInflamed);
save('out/omicInflamed.mat', 'omicInflamed');
% Recovered with tibolone
tib_abundance = abundanceTable(:, tib_samples);
tib_expression = expressionTable(:, tib_samples_exp);
omicTibolone = omicIntegrationPCA(abundanceTable, expressionTable);
omicTibolone = log(abs(min(omicTibolone)) + omicTibolone);
save('out/omicTibolone.mat', 'omicTibolone');