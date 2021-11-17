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

% load stoich.mat
% save stoich.mat

% load insideFixLeaks.mat
model = removeRxns(model, 'steroids');
[model, ModifiedRxns] = changeRxnMets(model, {'h2o[r]', 'h[r]'}, {'h2o[r]', 'h[r]'}, 'HMR_6813', [0, 0]);
[model, ModifiedRxns] = changeRxnMets(model, {'h[c]'}, {'h2o[c]'}, 'HMR_7656', [-1]);
[model, ModifiedRxns] = changeRxnMets(model, {'h[c]'}, {'h2o[c]'}, 'AGPAT2', [-1]);
[model, ModifiedRxns] = changeRxnMets(model, {'alpa_hs[c]'}, {'alpa_hs[c]', 'o2[c]', 'co2[c]', 'h[c]', 'pi[c]'}, 'AGPAT2', [-3, -8, 23, 35, 2]');
[model, ModifiedRxns] = changeRxnMets(model, {'h[c]', 'h2o[c]'}, {'h[c]', 'h2o[c]'}, 'HMR_1186', [-2, 1]');
[model, ModifiedRxns] = changeRxnMets(model, {'h[m]'}, {'h[m]'}, 'HMR_9818', [4517]');
% [model, ModifiedRxns] = changeRxnMets(model, {'o2[c]'}, {'o2[c]'}, 'HMR_5168', [2336/4]');
% [model, ModifiedRxns] = changeRxnMets(model, {'o2[c]'}, {'o2[c]'}, 'HMR_9818', [2336/4]');

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
    Indexes = find(siphonMetBool | leakMetBool);
    TestMet = false(size(model.mets));
    for i=1:length(Indexes)
        TestMet(Indexes(i)) = true;

        leakParams.epsilon=1e-4;
        minLeakParams.eta = getCobraSolverParams('LP', 'feasTol')*100;
        leakParams.method='dc';
        printLevel=1;
        minLeakParams.monoMetMode=1;
        [minLeakMetBool,minLeakRxnBool,minSiphonMetBool,minSiphonRxnBool,leakY,siphonY,statp,statn] =...
            findMinimalLeakageModeMet(model,TestMet,model.SIntRxnBool,modelBoundsFlag,minLeakParams,printLevel);
        
        if any(any(full(minLeakMetBool)))
            break;
        end
    end

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
% This function assembles a table with:
% 

    Table = array2table(RxnInd);
    Table.Rxns = model.rxns(RxnInd);
    model = extractSubNetwork(model, Table.Rxns);
    RxnInd = findRxnIDs(model, Table.Rxns);

    % Remove unknown elements
    % model.Elements = model.Elements(~ismember(model.Elements, {'X', 'R', 'FULLR'}));

    % Table.Formula = printRxnFormula(model, 'rxnAbbrList', Table.Rxns, 'printFlag', false);
    for i=1:length(model.Elements)
        element = model.Elements{i};
        balance = checkBalance(model, element, 0);
        Table.(element) = balance(RxnInd);
    end

    % Remove balanced reactions
    Table = Table(sum(Table{:, 4:end}') ~= 0, :);

end

function Tables = getCommonMetRxns(model, metabolite, Table)
    
    % find reactions with common metabolites
    [Table.Mets, Table.Stoichiometries]  ...
        = cellfun(@(x) findMetsFromRxns(model, x), Table.Rxns, 'UniformOutput',false);
    Table.Mets = cellfun(@(x) regexprep(x{1}, '\[.*\]', ''), Table.Mets, 'UniformOutput', false);
    Table.Stoichiometries = cellfun(@(x) full(x{1}), Table.Stoichiometries, 'UniformOutput', false);
    % Remove current metabolite compartment
    metabolite = regexprep(metabolite, '\[.*\]', '');
    % Remove irrelevant metabolites
    Table.Mets = cellfun(@(x) x(~ismember(x, metabolite)), Table.Mets, 'UniformOutput', false);
    Table.Stoichiometries = cellfun(@(met, stoich) ...
        stoich(~ismember(met, metabolite)), Table.Mets, Table.Stoichiometries, 'UniformOutput', false);
    Table.Mets = cellfun(@(x) x(~ismember(x, {'h2o'})), Table.Mets, 'UniformOutput', false);
    Table.Stoichiometries = cellfun(@(met, stoich) ...
        stoich(~ismember(met, {'h2o'})), Table.Mets, Table.Stoichiometries, 'UniformOutput', false);
    Table.Mets = cellfun(@(x) x(~ismember(x, {'o2'})), Table.Mets, 'UniformOutput', false);
    Table.Stoichiometries = cellfun(@(met, stoich) ...
        stoich(~ismember(met, {'o2'})), Table.Mets, Table.Stoichiometries, 'UniformOutput', false);
    Table.Mets = cellfun(@(x) x(~ismember(x, {'h'})), Table.Mets, 'UniformOutput', false);
    Table.Stoichiometries = cellfun(@(met, stoich) ...
        stoich(~ismember(met, {'h'})), Table.Mets, Table.Stoichiometries, 'UniformOutput', false);

    Tables = {};
    
    % Getr products and reactants
    products = cellfun(@(met, stoich) ...
        met(stoich > 0),Table.Mets, Table.Stoichiometries, 'UniformOutput', false);
    reactants = cellfun(@(met, stoich) ...
        met(stoich < 0),Table.Mets, Table.Stoichiometries, 'UniformOutput', false);

    % Find groups of reactions with common names
    commonMets = zeros(height(Table));
    for i=1:height(Table)
        commonMets(i,:) = cellfun(@(x) sum(ismember(reactants{i}, x)), products);
    end

    while true
        [value, xindex] = max(commonMets);
        [value, yindex] = max(value);
        if any(sum(commonMets == value) > 1) || isempty(Table)
            break;
        end
        xindex = xindex(yindex);
        fprintf('Value %d [%d,%d]\r\n', value, xindex, yindex)
        Tables{end + 1} = Table([xindex, yindex],:);
        commonMets([xindex, yindex], :) = 0;
        commonMets(:, [xindex, yindex]) = 0;
    end

    rxn_groups = unique(commonMets >= 1, 'rows');
    rxn_groups = rxn_groups(sum(rxn_groups, 2) == 2, :);

    % if isempty(rxn_groups)
    %     commonMets = zeros(height(Table));
    %     for i=1:height(Table)
    %         commonMets(i,:) = cellfun(@(x) sum(ismember(Table.Mets{i}, x)), Table.Mets);
    %     end
    %     rxn_groups = unique(commonMets >= 1, 'rows');
    % end
    
    % Remove reactions that don't have common metabolites
    for i=1:size(rxn_groups, 1)
        Tables{end + 1} = Table(rxn_groups(i,:), :);
    end


end

function weights = get_reactions_weights(model, Table)
    % Get metabololites and stoichiometries
    [Table.Mets, Table.Stoichiometries]  ...
        = cellfun(@(x) findMetsFromRxns(model, x), Table.Rxns, 'UniformOutput',false);
    Table.Mets = cellfun(@(x) regexprep(x{1}, '\[.*\]', ''), Table.Mets, 'UniformOutput', false);
    Table.Stoichiometries = cellfun(@(x) full(x{1}), Table.Stoichiometries, 'UniformOutput', false);

    % Split rxns in by reactants and products imbalance
    imbalance = sum(table2array(Table(:,4:end-2)), 2);
    reactants = imbalance < 0;
    products = imbalance > 0;

    % Get reactant and product metabolites
    Table.Mets(reactants) = cellfun(@(met, stoich) ...
        met(stoich < 0),Table.Mets(reactants), Table.Stoichiometries(reactants), 'UniformOutput', false);
    Table.Mets(products) = cellfun(@(met, stoich) ...
        met(stoich > 0),Table.Mets(products), Table.Stoichiometries(products), 'UniformOutput', false);
    Table.Stoichiometries(reactants) = cellfun(@(stoich) ...
        stoich(stoich < 0), Table.Stoichiometries(reactants), 'UniformOutput', false);
    Table.Stoichiometries(products) = cellfun(@(stoich) ...
        stoich(stoich > 0), Table.Stoichiometries(products), 'UniformOutput', false);

    % Keep only common metabolites between reactants and products
    for i=1:height(Table)

        commonReactants = ismember(Table.Mets{i}, cat(1, Table.Mets{reactants}));
        Table.Stoichiometries{i} = Table.Stoichiometries{i}(commonReactants);
        Table.Mets{i} = Table.Mets{i}(commonReactants);

        commonPoducts = ismember(Table.Mets{i}, cat(1, Table.Mets{products}));
        Table.Stoichiometries{i} = Table.Stoichiometries{i}(commonPoducts);
        Table.Mets{i} = Table.Mets{i}(commonPoducts);
    end
end


function met_candidate = findMetCandidate(model, metabolite, compartment)
    met_no_num = regexprep(model.metFormulas, '\d', '');
    
    met_candidates = model.mets(ismember(met_no_num, metabolite));
    if isempty(met_candidates)
        met_candidate = {};
    else
        met_candidate = met_candidates(contains(met_candidates, compartment));
    end
    % Return only one candidate
    if length(met_candidate) > 1
        met_candidate = met_candidate(1);
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
        difference = zeros(1, length(MetElemet));
        for m=1:length(MetElemet)
            difference(m) = sum(Table{:, MetElemet(m)});
            pairs = nchoosek(Table{:, MetElemet(m)}, 2);
            balanced = sum(pairs') == 0;
            disRxns(m,:) = ~ismember(Table{:, MetElemet(m)}, pairs(balanced, :));
        end
        if length(MetElemet) > 1
            Table = Table(any(disRxns), :);
        else
            Table = Table(disRxns >= 1, :);
        end

        if height(Table) <= 1
            break;
        end
        
        met_compartment = regexp(metabolite, '\[.*\]', 'match');
        met_compartment = met_compartment{1}{1};
        balance = -difference./(2 * Ematrix);
        
        % Once we have found the disblanced reactions
        % Fix the stoichiometry finding a midle point between those
        for m=1:length(MetElemet)
            % Find the metabolite candidate
            met_candidate = findMetCandidate(model, MetElemet(m), met_compartment);
            if isempty(met_candidate)
                met_candidate = metabolite;
            end

            met_candidateInd = findMetIDs(model, met_candidate);

            stoich_change = zeros(height(Table), 1) + balance(m);
            stoich_change = model.S(met_candidateInd, Table.RxnInd) + stoich_change';

            % Change the stoichiometry in the model
            [model, ~] = changeRxnMets(model, met_candidate, met_candidate, Table.Rxns, stoich_change);
            for i=1:height(Table)
                fprintf('Changing %s with %d %s\n', Table.Rxns{i}, stoich_change(i), met_candidate{1})
            end
            
            if isempty(findMetCandidate(model, MetElemet(m), met_compartment))
                [candidate_Ematrix, candidate_elements] = getElementalComposition(model.metFormulas{MetInd}, model.Elements);
                candidate_elements = candidate_elements(candidate_Ematrix > 0);
                candidate_Ematrix = candidate_Ematrix(candidate_Ematrix > 0);
                
                otherEmatrix = candidate_Ematrix(1:length(candidate_elements) ~= m);
                otherElements = candidate_elements(1:length(candidate_elements) ~= m);
                otherMet = findMetCandidate(model, otherElements, met_compartment);
                
                met_candidateInd = findMetIDs(model, otherMet);
                stoich_change = (zeros(height(Table), 1) + (balance(m)));
                stoich_change = (model.S(met_candidateInd, Table.RxnInd)) - stoich_change';
                [model, ~] = changeRxnMets(model, otherMet, otherMet, Table.Rxns, stoich_change);
                break;
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

    % Generate table with stats of the leaking reactions
    Table = getStoichTable(model, RxnInd)
    % Remove balanced reactions
    % Table = Table(sum(Table(5:end)),:)
    if height(Table) > 1
        Tables = getCommonMetRxns(model, metabolite, Table);
        if isempty(Tables)
            leakInd = find(minSiphonMetBool(MetInd,:));
            RxnInd = find(minSiphonRxnBool(:, leakInd));
            Table = getStoichTable(model, RxnInd)
            Tables = getCommonMetRxns(model, metabolite, Table);
        end
        
        for i=1:length(Tables)
            model = balance(model, metabolite, Tables{i});
            % model = balance(model, {'h[c]'}, Tables{i});
            % model = balance(model, {'co2[c]'}, Tables{i});
        end
        
    else
        break;
    end
    % RxnBool = false(size(model.rxns));
    % RxnBool(Table.RxnInd) = true;
    % [minLeakMetBool,minLeakRxnBool,minSiphonMetBool,minSiphonRxnBool, ~, ~, ~, ~] =...
    %     findMinimalLeakageModeMet(model, TestMet, RxnBool,modelBoundsFlag,minLeakParams,printLevel);

end

end