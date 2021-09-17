function abundanceTable = mapAbundance (base_model, abundance, proteinECNumbers)

abundanceTable = table();
abundance = addvars(abundance, proteinECNumbers, 'After', 1);
for k=3:length(abundance.Properties.VariableNames)
    abundanceToMap =  abundance(:, [2 k]);
    abundanceToMap.Properties.VariableNames = {'id' 'value'};
    
    % Remove missing values
    abundanceToMap = rmmissing(abundanceToMap);
    abundanceToMap = abundanceToMap(~cellfun('isempty', abundanceToMap.id), :);
    
    [abundanceRxns, ~] = mapAbundanceToReactions(base_model, abundanceToMap);
    abundanceTable = addvars(abundanceTable, abundanceRxns);
end

abundanceTable.Properties.VariableNames = abundance.Properties.VariableNames(3:end);

end