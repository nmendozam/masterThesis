close all force; clear variables; clc
%% Startup the COBRA Toolbox
addpath('cobratoolbox')
initCobraToolbox(false) % Don't update the toolbox
changeCobraSolver ('gurobi', 'all', 1); % For large models
clc
%% Load data and model
addpath('functions') 
load('out/dirtyAstrociteModel.mat');
model = astrociteModel;
[nMets, nRxns, nCtrs, nVars, nGenes, nComps] = getModelSizes(model);
fprintf('%6s\t%6s\t%6s\t%6s%6s\t%6s\n','#mets','#rxns', '#nCntrs', '#nVars', '#nGenes', '#nComps')
fprintf('%6u\t%6u\t%6u\t%6u\t%6u\t%6u\t%s\n',nMet,nRxn, nCtrs, nVars, nGenes, nComps,' totals.')

properties_model = 'out/model_properties.mat';
if isfile(properties_model )
    load(properties_model);
else
    model = checkModelProperties(model);
    save(properties_model, 'model');
end

[SConsistentMetBool, SConsistentRxnBool, SInConsistentMetBool, ...
    SInConsistentRxnBool, unknownSConsistencyMetBool, ...
    unknownSConsistencyRxnBool, model] = findStoichConsistentSubset(model);


%% Set parameters
printLevel=2;
modelBoundsFlag = 1;
leakParams.epsilon = (getCobraSolverParams('LP', 'feasTol')*100);
% leakParam.eta = feasTol*100;
leakParams.theta = 0.5;
leakParams.method = 'quasiConcave';

%% Find leakage or siphons in heuristically internal part using the bounds given with the model
if 1
    [leakMetBool,leakRxnBool,siphonMetBool,siphonRxnBool,leakY,siphonY,statp,statn]...
        = findMassLeaksAndSiphons(model,model.SIntMetBool,model.SIntRxnBool,modelBoundsFlag,leakParams,printLevel);
end
%% For each leaking metabolite find a minimal cardinality leakage mode
leakParams.epsilon=1e-4;
minLeakParams.eta = getCobraSolverParams('LP', 'feasTol')*100;
leakParams.method='dc';

if 1
    minLeakParams.monoMetMode=1;
    [minLeakMetBool,minLeakRxnBool,minSiphonMetBool,minSiphonRxnBool,leakY,siphonY,statp,statn] =...
        findMinimalLeakageModeMet(model,leakMetBool,model.SIntRxnBool,modelBoundsFlag,minLeakParams,printLevel);
end
%% For each siphon metabolite find a minimal cardinality siphon mode
if 1
    minLeakParams.monoMetMode=1;
    [minLeakMetBool,minLeakRxnBool,minSiphonMetBool,minSiphonRxnBool,leakY,siphonY,statp,statn] =...
        findMinimalLeakageModeMet(model,siphonMetBool,model.SIntRxnBool,modelBoundsFlag,minLeakParams,printLevel);
end

%% Print reactions
LeakMetsInd = find(any(minLeakMetBool'));

%% Test if corrected
MetInd = LeakMetsInd(1);
TestMet = false(size(siphonMetBool));
TestMet(MetInd) = true;
Met = model.mets(MetInd);
[Ematrix, elements] = getElementalComposition(model.metFormulas{MetInd}, model.Elements);
MetElemet = elements(Ematrix > 0);
while true
    % Find minimal leakage model for the metabolite of interest
    [minLeakMetBool,minLeakRxnBool,minSiphonMetBool,minSiphonRxnBool, ~, ~, ~, ~] =...
        findMinimalLeakageModeMet(model, TestMet,model.SIntRxnBool,modelBoundsFlag,minLeakParams,printLevel);
    if sum(minLeakRxnBool)
        leakInd = find(minLeakMetBool(MetInd,:));
        RxnInd = find(minLeakRxnBool(:, leakInd));
    elseif sum(minSiphonRxnBool)
        leakInd = find(minSiphonMetBool(MetInd,:));
        RxnInd = find(minSiphonRxnBool(:, leakInd));
    else
        break;
    end

    % Generate table with stats of the leaking reactions
    Table = array2table(RxnInd);
    Table.Rxns = model.rxns(RxnInd);
    Table.Formula = printRxnFormula(model, 'rxnAbbrList', Table.Rxns, 'printFlag', false);
    for i=1:length(model.Elements)
        element = model.Elements{i};
        balance = checkBalance(model, element, 0);
        Table.(element) = balance(Table.RxnInd);
    end

    % Remove balanced reactions
    Table = Table(sum(Table{:, 4:end}') ~= 0, :);
    % if ~any(sum(Table{:, MetElemet}))
    %     leakInd = find(minSiphonMetBool(MetInd,:));
    %     RxnInd = find(minSiphonRxnBool(:, leakInd));
    %     Table = Table(sum(Table{:, 4:end}') ~= 0, :);
    % end
    while any(sum(Table{:, MetElemet}))
        % Search pair of disbalanced reaction
        difference = sum(Table{:, MetElemet});
        pairs = nchoosek(Table{:, MetElemet}, 2);
        balanced = find(sum(pairs') == 0);
        disRxns = ~ismember(Table{:, MetElemet}, pairs(balanced, :));
        Table = Table(disRxns, :);

        % Onece we have found the disblanced reactions
        % Fix the stoichiometry finding a midle point between those
        stoich_change = ones(1, height(Table));
        balance = -difference/(2 * sum(Ematrix));
        % Balance values for each side of the reaction
        stoich_change(Table.H > 0) = balance * Table.H(Table.H > 0) / sum(Table.H(Table.H > 0));
        stoich_change(Table.H < 0) = balance * Table.H(Table.H < 0) / sum(Table.H(Table.H < 0));
        stoich_change = model.S(MetInd, Table.RxnInd) + stoich_change;

        % Change the stoichiometry in the model
        [model, ModifiedRxns] = changeRxnMets(model, Met, Met, Table.Rxns, stoich_change);
        for i=1:height(Table)
            fprintf('Changing %s with %d\n', Table.Rxns{i}, stoich_change(i))
        end

        % Update table
        Table.Formula = printRxnFormula(model, Table.Rxns);
        for i=1:length(model.Elements)
            element = model.Elements{i};
            balance = checkBalance(model, element, 0);
            Table.(element) = balance(Table.RxnInd);
        end
    end
    disp('Balanced looking if there is other leak or siphon reactions...')
end
RxnInterest = 25;
Form = model.metFormulas(findMetIDs(model, findMetsFromRxns(model, Rxns(RxnInterest))));
mets = findMetsFromRxns(model, Rxns(RxnInterest));
{Form{:};mets{:}}
stoich({Form{:}, 'H2O', 'H'})
stoich(model.metFormulas(findMetIDs(model, findMetsFromRxns(model, Rxns(3)))))
stoich({'C21H36N7O16P3S','C5H9O8P' ,'C21H41O7P' ,'C37H62N7O17P3S', 'H2O', 'H', 'O2', 'CO2', 'Na'})


% Repair for h20[c]
[model, ModifiedRxns] = changeRxnMets(model, {'h2o[r]', 'h[r]'}, {'h2o[r]', 'h[r]'}, 'HMR_6813', [0, 0]);
[model, ModifiedRxns] = changeRxnMets(model, {'h[c]'}, {'h2o[c]'}, 'HMR_7656', [-1]);
[model, ModifiedRxns] = changeRxnMets(model, {'h[c]'}, {'h2o[c]'}, 'AGPAT2', [-1]);
[model, ModifiedRxns] = changeRxnMets(model, {'alpa_hs[c]'}, {'alpa_hs[c]', 'o2[c]', 'co2[c]', 'h[c]', 'pi[c]'}, 'AGPAT2', [-3, -8, 23, 35, 2]');
[model, ModifiedRxns] = changeRxnMets(model, {'h[c]', 'h2o[c]'}, {'h[c]', 'h2o[c]'}, 'HMR_1186', [-2, 1]');
% Repair for2 o2[c]
[model, ModifiedRxns] = changeRxnMets(model, {'h[c]'}, {'h[c]'}, 'AGPAT1', [-2]');
[model, ModifiedRxns] = changeRxnMets(model, {'o2[c]'}, {'o2[c]'}, 'IDL_HSDEG', [-8]');
% [model, ModifiedRxns] = changeRxnMets(model, {'o2[c]'}, {'o2[c]'}, 'HDL_HSDEG', [0]');
[model, ModifiedRxns] = changeRxnMets(model, {'o2[c]'}, {'o2[c]'}, 'HDL_HSSYN', [-2]');
[model, ModifiedRxns] = changeRxnMets(model, {'o2[c]'}, {'o2[c]'}, 'LDL_HSDEG', [-2]');
[model, ModifiedRxns] = changeRxnMets(model, {'o2[c]'}, {'o2[c]'}, 'CHYLO_HSDEG', [-2]');

%% get rigth formulas
model.metFormulas([5450, 5451, 5453, 5455, 5457, 5459, 5461, 5463, 5465, 5467, 5469, 5471, 5473, 5475, 5477, 5479, 5481, 5483, 5485, 5487]) = '';
computeMetFormulae(model)


%% For each leaking metabolite find a minimal cardinality leakage mode
if 1
    minLeakParams.monoMetMode=0;
    [minLeakMetBool,minLeakRxnBool,minSiphonMetBool,minSiphonRxnBool,leakY,siphonY,statp,statn] =...
        findMinimalLeakageModeMet(model,leakMetBool,model.SIntRxnBool,modelBoundsFlag,minLeakParams,printLevel);
end
%% For each siphon metabolite find a minimal cardinality siphon mode
if 1
    minLeakParams.monoMetMode=0;
    [minLeakMetBool,minLeakRxnBool,minSiphonMetBool,minSiphonRxnBool,leakY,siphonY,statp,statn] =...
        findMinimalLeakageModeMet(model,siphonMetBool,model.SIntRxnBool,modelBoundsFlag,minLeakParams,printLevel);
end
%% For each heuristically internal but stoichiometrically inconsistent reaction (one at a time), find the min cardinality leakage mode
if 1
    rxnBool=model.SIntRxnBool & model.SInConsistentRxnBool;
    metBool=true(nMet,1);
    minLeakParams.monoRxnMode=1;
    [minLeakMetBool,minLeakRxnBool,minSiphonMetBool,minSiphonRxnBool,leakY,siphonY,statp,statn] =...
        findMinimalLeakageModeRxn(model,rxnBool,metBool,modelBoundsFlag,minLeakParams,printLevel);
end

%% For each heuristically internal but stoichiometrically inconsistent reaction, find the min cardinality leakage mode
if 1
    rxnBool=model.SIntRxnBool & model.SInConsistentRxnBool;
    metBool=true(nMet,1);
    minLeakParams.monoRxnMode=0;
    [minLeakMetBool,minLeakRxnBool,minSiphonMetBool,minSiphonRxnBool,leakY,siphonY,statp,statn] =...
        findMinimalLeakageModeRxn(model,rxnBool,metBool,modelBoundsFlag,minLeakParams,printLevel);
end
%% For each heuristically internal but unknown stoichiometric consistency reaction (one at a time), find a minimal cardinality leakage mode
if 1
    rxnBool=model.SIntRxnBool & model.unknownSConsistencyRxnBool;
    metBool=true(nMet,1);
    minLeakParams.monoRxnMode=1;
    [minLeakMetBool,minLeakRxnBool,minSiphonMetBool,minSiphonRxnBool,leakY,siphonY,statp,statn] = findMinimalLeakageModeRxn(model,rxnBool,metBool,modelBoundsFlag,minLeakParams,printLevel);
end
%% For each heuristically internal but unknown stoichiometric consistency reaction, find a minimal cardinality leakage mode
if 1
    rxnBool=model.SIntRxnBool & model.unknownSConsistencyRxnBool;
    metBool=true(nMet,1);
    minLeakParams.monoRxnMode=0;
    [minLeakMetBool,minLeakRxnBool,minSiphonMetBool,minSiphonRxnBool,leakY,siphonY,statp,statn] = findMinimalLeakageModeRxn(model,rxnBool,metBool,modelBoundsFlag,minLeakParams,printLevel);
end

if printLevel>0
    fprintf('%6u\t%6u\t%s\n',nnz(~model.SIntMetBool),nnz(~model.SIntRxnBool),' heuristically exchange.')
    fprintf('%6u\t%6u\t%s\n',nnz(model.SIntMetBool),nnz(model.SIntRxnBool),' All internally stoichiometrically consistent.');
    fprintf('%6u\t%6u\t%s\n',nnz(model.fluxConsistentMetBool),nnz(model.fluxConsistentRxnBool),' All flux consistent.');
    fprintf('%s\n','Input model assumed to be stoichiometrically and flux consistent.');
end