function only_in_model1_by_name = compareModelsRxns(model1, model2)

tf = ismember(model1.rxns, model2.rxns);
only_in_model1_by_name_indexes = find(~tf);
only_in_model1_by_name = model1.rxns(~tf);
% 
% % Get reactions in model1 that are missing in model2
% only_in_model1 = {};
% for i=1:length(only_in_model1_by_name)
%     model1Rxn = find(ismember(model1.rxns, only_in_model1_by_name(i)));
%     % Products and reactants from model1
%     model1Products = model1.mets(model1.S(:,model1Rxn) > 0);
%     model1Reactants = model1.mets(model1.S(:,model1Rxn) < 0);
%     
%     % Metabolite names of the reaction in model1
%     metsNames = model1.mets(model1.S(:,model1Rxn) ~= 0);
%     % Inexes of same metabolites but in model 2
%     metsNum = ismember(model2.mets, metsNames);
%     % Get reactions with the metabolites form model2
%     [~, rxnsWithMet] = find(model2.S(metsNum,:) ~= 0);
%     rxnsWithMet = unique(rxnsWithMet);
%     
%     
%     
%     % Compare the reactions that contain the metabolites
%     It_has_other_name =  false;
%     for j=1:length(rxnsWithMet)
%         % Products and reactants from model2
%         model2Products = sort(model2.mets(model2.S(:,rxnsWithMet(j)) > 0));
%         model2Reactants = sort(model2.mets(model2.S(:,rxnsWithMet(j)) < 0));
%         
%            % Equal number of products, Equal number of reactants and the
%            % same products and reactants
%         if ~isempty(model2Products) && length(model2Products) == length(model1Products) ...
%            && ~isempty(model2Reactants) && length(model2Reactants) == length(model1Reactants) ...
%            && all(ismember(model2Products, model1Products)) && all(ismember(model2Reactants, model1Reactants)) ...
%            && (model1.lb(only_in_model1_by_name_indexes(i)) < 0) == (model2.lb(rxnsWithMet(j)) < 0)
%        
%             It_has_other_name = true;
%             printRxnFormula(model1, 'rxnAbbrList', only_in_model1_by_name(i));
%             printRxnFormula(model2, 'rxnAbbrList', model2.rxns(rxnsWithMet(j)));
%         end
%     end
%     if ~It_has_other_name
%         only_in_model1 = [only_in_model1(:)', only_in_model1_by_name{i}];
%     end
% end

end