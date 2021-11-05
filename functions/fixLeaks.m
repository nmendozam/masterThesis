function model = fixLeaks(model)
%fixLeaks - Description
%
% Syntax: model = fixLeaks(input)
%
% Long description

% properties_model = 'out/model_properties.mat';
% if isfile(properties_model )
%     load(properties_model);
% else
%     model = checkModelProperties(model);
%     save(properties_model, 'model');
% end

% printLevel=2;
% massBalanceCheck=1;
% [SConsistentMetBool, SConsistentRxnBool, SInConsistentMetBool, SInConsistentRxnBool, ...
%     unknownSConsistencyMetBool, unknownSConsistencyRxnBool, model, stoichConsistModel] ...
%     = findStoichConsistentSubset(model, massBalanceCheck, printLevel);

% model = removeRxns(model, model.rxns(model.rxnRemoveBool));

load stoich.mat
% save stoich.mat

% load insideFixLeaks.mat
model = removeRxns(model, 'steroids');
[model, ModifiedRxns] = changeRxnMets(model, {'h2o[r]', 'h[r]'}, {'h2o[r]', 'h[r]'}, 'HMR_6813', [0, 0]);
[model, ModifiedRxns] = changeRxnMets(model, {'h[c]'}, {'h2o[c]'}, 'HMR_7656', [-1]);
[model, ModifiedRxns] = changeRxnMets(model, {'h[c]'}, {'h2o[c]'}, 'AGPAT2', [-1]);
[model, ModifiedRxns] = changeRxnMets(model, {'alpa_hs[c]'}, {'alpa_hs[c]', 'o2[c]', 'co2[c]', 'h[c]', 'pi[c]'}, 'AGPAT2', [-3, -8, 23, 35, 2]');
[model, ModifiedRxns] = changeRxnMets(model, {'h[c]', 'h2o[c]'}, {'h[c]', 'h2o[c]'}, 'HMR_1186', [-2, 1]');
[model, ModifiedRxns] = changeRxnMets(model, {'h[m]'}, {'h[m]'}, 'HMR_9818', [4517]');

%% Test if corrected
while true
    %% Find leakage or siphons in heuristically internal part using the bounds given with the model
    modelBoundsFlag = 1;
    leakParams.epsilon = (getCobraSolverParams('LP', 'feasTol')*100);
    leakParams.theta = 0.5;
    leakParams.method = 'quasiConcave';
    printLevel=2;
    [leakMetBool,leakRxnBool,siphonMetBool,siphonRxnBool,leakY,siphonY,statp,statn]...
        = findMassLeaksAndSiphons(model,model.SIntMetBool,model.SIntRxnBool,modelBoundsFlag,leakParams,printLevel);

    %% For each leaking and siphon metabolite find a minimal cardinality leakage mode
    leakParams.epsilon=1e-4;
    minLeakParams.eta = getCobraSolverParams('LP', 'feasTol')*100;
    leakParams.method='dc';
    printLevel=1;
    minLeakParams.monoMetMode=1;
    [minLeakMetBool,minLeakRxnBool,minSiphonMetBool,minSiphonRxnBool,leakY,siphonY,statp,statn] =...
        findMinimalLeakageModeMet(model,(siphonMetBool | leakMetBool),model.SIntRxnBool,modelBoundsFlag,minLeakParams,printLevel);

    %% Balance first metabolite
    LeakMets = model.mets(any(minLeakMetBool'));
    model = balanceMetabolite(model, LeakMets(1));
%     save stoich.mat
end

end

function Table = getStoichTable(model, RxnInd)
%getStoichTable - Gets stoichiometric table with Formula and balance of elements
%
% Syntax: Table = getStoichTable(model)
%
% This function assemples a table with:
% 

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

end

function Tables = getCommonMetRxns(model, metabolite, Table)
    
    % find reactions with common metabolites
    Table.Mets = cellfun(@(x) findMetsFromRxns(model, x), Table.Rxns, 'UniformOutput',false);
    Table.Mets = cellfun(@(x) regexprep(x, '\[.*\]', ''), Table.Mets, 'UniformOutput', false);
    % Remove current metabolite
    metabolite = regexprep(metabolite, '\[.*\]', '');
    Table.Mets = cellfun(@(x) x(~ismember(x, metabolite)), Table.Mets, 'UniformOutput', false);
    
    % Find groups of reactions with common names
    commonMets = cell(1, height(Table));
    for i=1:height(Table)
        commonMets{i} = cellfun(@(x) any(ismember(Table.Mets{i}, x)), Table.Mets);
    end

    m = cat(3, commonMets{:});
    m = reshape(m, [], size(m, 3)).';
    rxn_groups = unique(m, 'rows');

    % Remove reactions that don't have common metabolites
    Tables = cell(1, size(rxn_groups, 1));
    for i=1:size(rxn_groups, 1)
        Tables{i} = Table(rxn_groups(i,:), :);
    end

end

function weights = get_reactions_weights(model, Table)
    Table.Mets = cellfun(@(x) findMetsFromRxns(model, x), Table.Rxns, 'UniformOutput',false);
    Table.Mets = cellfun(@(x) regexprep(x, '\[.*\]', ''), Table.Mets, 'UniformOutput', false);
    for i=1:height(Table)
        commonMets{i} = ismember(Table.Mets{i}, cat(1, Table.Mets{Table.H < 0}));
        Table.Mets{i}(x)
    end
end

function model = balance(model, metabolite, Table)
    MetInd = find(strcmp(metabolite, model.mets));
    [Ematrix, elements] = getElementalComposition(model.metFormulas{MetInd}, model.Elements);
    MetElemet = elements(Ematrix > 0);
    Ematrix = Ematrix(Ematrix > 0);

    while any(sum(Table{:, MetElemet}))
        % Search pair of disbalanced reaction
        disRxns = zeros(length(MetElemet), height(Table));
        difference = zeros(length(MetElemet),1);
        for m=1:length(MetElemet)
            difference(m) = sum(Table{:, MetElemet(m)});
            pairs = nchoosek(Table{:, MetElemet(m)}, 2);
            balanced = sum(pairs') == 0;
            disRxns(m,:) = ~ismember(Table{:, MetElemet(m)}, pairs(balanced, :));
        end
        Table = Table(any(disRxns), :);

        % Once we have found the disblanced reactions
        % Fix the stoichiometry finding a midle point between those
        for m=1:length(MetElemet)
            stoich_change = zeros(height(Table), 1);
            balance = -difference(m)/(2 * Ematrix(m));
            % Balance values for each side of the reaction
            MetColum = Table.(MetElemet{m});
            stoich_change(MetColum > 0) = balance * MetColum(MetColum > 0) / sum(MetColum(MetColum > 0));
            stoich_change(MetColum < 0) = balance * MetColum(MetColum < 0) / sum(MetColum(MetColum < 0));
            stoich_change = model.S(MetInd, Table.RxnInd) + stoich_change';
            
            % Find the metabolite candidate
            met_compartment = regexp(metabolite, '\[.*\]', 'match');
            met_candidates = model.mets(ismember(model.metFormulas, MetElemet(1)));
            met_candidate = met_candidates{contains(met_candidates, met_compartment{1}{1})};

            % Change the stoichiometry in the model
            [model, ModifiedRxns] = changeRxnMets(model, {met_candidate}, {met_candidate}, Table.Rxns, stoich_change);
            for i=1:height(Table)
                fprintf('Changing %s with %d %d\n', Table.Rxns{i}, stoich_change(i), met_candidate)
            end
        end

        % Update table
        Table = getStoichTable(model, Table.RxnInd);
    end
    disp('Balanced looking if there is other leak or siphon reactions...')
end


function model = balanceMetabolite(model, metabolite)
%balanceMetabolite - Description
%
% Syntax: model = balanceMetabolite(model)
%
% Long description

printLevel=1;
modelBoundsFlag = 1;
minLeakParams.eta = getCobraSolverParams('LP', 'feasTol')*100;
minLeakParams.monoMetMode=1;


MetInd = find(strcmp(metabolite, model.mets));
TestMet = false(size(model.mets));
TestMet(MetInd) = true;

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
    if (strcmp(metabolite, 'h[m]') && sum(minSiphonRxnBool) >= 32)
        disp('pare aqui')
    end

    % Generate table with stats of the leaking reactions
    Table = getStoichTable(model, RxnInd);
    if height(Table) > 1
        Tables = getCommonMetRxns(model, metabolite, Table);
        for i=1:length(Tables)
            model = balance(model, metabolite, Tables{i});
        end
    else
        model = balance(model, metabolite, Table);
    end
    % If this are balanced use the sink
%     if ~any(sum(Table{:, MetElemet}))
%         leakInd = find(minSiphonMetBool(MetInd,:));
%         RxnInd = find(minSiphonRxnBool(:, leakInd));
%         Table = getStoichTable(model, RxnInd);
%     end
%     if height(Table) == 1
%         mets = findMetsFromRxns(model, Table.Rxns);
%         Form = model.metFormulas(findMetIDs(model, mets));
%         Rs = array2table(sum(stoich(Form)')');
%         Rs.Formulas = Form;
%         Rs.Names = mets
%     else
        
%     end

end

end