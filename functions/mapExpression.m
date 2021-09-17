function expressionTable = mapExpression(base_model, expresion)

base_model.genes = regexprep(base_model.genes,"\.[0-9]*","");

expresion_columns = expresion.Properties.VariableNames;
expressionTable = table();
for k=3:length(expresion_columns)
    expressionToMap =  expresion(:, [2 k]);
    expressionToMap.Properties.VariableNames = {'gene' 'value'};
    expressionToMap.gene(cellfun('isempty', expressionToMap.gene)) = {' '};
    [expressionRxns, ~] = mapExpressionToReactions(base_model, expressionToMap);
    expressionTable = addvars(expressionTable, expressionRxns);
end
expressionTable.Properties.VariableNames = expresion.Properties.VariableNames(3:end);

end