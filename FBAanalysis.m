close all force; clear variables; clc
%% Startup the COBRA Toolbox
addpath('cobratoolbox')
addpath('functions')  
initCobraToolbox(false) % Don't update the toolbox
changeCobraSolver ('gurobi', 'all', 1); % For large models
clc

%% Load models
MenHealthyModel = readCbModel('GSMs/Astrocyte_Healty_Mendoza2022.xml');
solFBA = optimizeCbModel(MenHealthyModel, 'max'); 
fprintf("The FBA flux of healthy model %f\n", solFBA.f)

MenInflamedModel = readCbModel('GSMs/Astrocyte_InflamedPalmitate_Mendoza2022.xml');
solFBA = optimizeCbModel(MenInflamedModel, 'max'); 
fprintf("The FBA flux of inflamed model %f\n", solFBA.f)

MenTreatedModel = readCbModel('GSMs/Astrocyte_TreatedTibolone_Mendoza2022.xml');
solFBA = optimizeCbModel(MenTreatedModel, 'max'); 
fprintf("The FBA flux of treated model %f\n", solFBA.f)

%% See esential reacions
allModels = {};
newModelName = 'AstrociteDMEM';
allModels.(newModelName) = astrociteModelDMEM;
newModelName = 'Osorio';
allModels.(newModelName) = Osorio_model;
[essentialRxn4Models, dataStruct] = essentialRxn4MultipleModels(allModels , 'biomass_reaction');


%% Model statistics
n_Exchanges = length(findExchangeReactions(model));
n_Transports = length(findTransRxns(model));
fprintf("Exchanges %d Transports %d comparimetized %d of %d total\n", n_Exchanges, n_Transports, length(model.rxns) - (n_Transports + n_Exchanges), length(model.rxns))

% Reactions enzimes
exp = { '^(?<number>\d)\.[A-Z]',  '^(?<number>\d)\.\d'};
for i=1:length(exp)
    tok = regexp(model.rxnECNumbers,exp{i}, 'tokens');
    x = [tok{~cellfun(@isempty, tok)}];
    x = [x{:}];
    y=unique(x);
    for k=1:length(y)
      freq(k)=sum(ismember(x,y(k)));
    end
    out=[cellstr(y)' num2cell(freq')]
end

% Reactions by compartment
ex = false(length(model.rxns),1); tr = false(length(model.rxns),1);
ex(findExchangeReactions(model)) = true;
tr(findRxnIDs(model, findTransRxns(model))) = true;
compartments = {'[c]','[e]','[g]','[i]','[l]','[m]','[n]','[r]','[x]'};
for i=1:length(compartments)
    x = findRxnFromCompartment(model, compartments{i});
    freq(i) = sum(ismember(model.rxns, x(:, 1)) & ~ex & ~tr);
end

