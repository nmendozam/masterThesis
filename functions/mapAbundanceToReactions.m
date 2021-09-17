function [abundanceRxns, parsedPR, proteins_used, signifRxns] = mapAbundanceToReactions(model, abundanceData, minSum)                                          
% Determines the expression data associated to each reaction present in
% the model 
%
% USAGE:
%
%    [expressionRxns parsedGPR, gene_used] = mapExpressionToReactions(model, expressionData) 
%    [expressionRxns, parsedGPR, gene_used, signifRxns] =  mapExpressionToReactions(model, expressionData, minSum)
%
% INPUTS:
%	model                   model strusture
%	expressionData          mRNA expression data structure
%       .gene               	cell array containing GeneIDs in the same
%                               format as model.genes
%       .value                  Vector containing corresponding expression
%                               value (FPKM/RPKM)
%       .sig:               [optional field] Vector containing significance values of
%                           expression corresponding to expression values in
%                           expressionData.value (ex. p-values)
%
% OPTIONAL INPUT:
%    minSum:         instead of using min and max, use min for AND and Sum
%                    for OR (default: false, i.e. use min)
%
% OUTPUTS:
%   expressionRxns:         reaction expression, corresponding to model.rxns.
%   parsedGPR:              cell matrix containing parsed GPR rule
%   gene_used:              gene identifier, corresponding to model.rxns, from GPRs
%                           whose value (expression and/or significance) was chosen for that
%                           reaction
%
% OPTIONAL OUTPUTS:
%   signifRxns:              significance of reaction expression, corresponding to model.rxns.

%
% Authors:
%       - Nicolas Mendoza Mejia, Apr 2020, Abundance mapping for reactions

if ~exist('minSum','var')
    minSum = false;
end

if isfield(abundanceData, 'sig') 
    exprSigFlag = 1; 
else
    exprSigFlag = 0;
end 

% Extracting GPR data from model
parsedPR = PRparser(model,minSum); % This could be integrated with GPRparser in a function called ruleParser


if exprSigFlag == 0
    % Find wich genes in expression data are used in the model
    [protein_id, protein_abun] = findUsedLevels(model.ECNumbers,abundanceData); %this can be unified with findUsedGenesLevels

    % Link the gene to the model reactions
    [abundanceRxns,  proteins_used] = selectGeneFromGPR(model, protein_id, protein_abun, parsedPR, minSum);
else
    
    [protein_id, protein_abun, gene_sig] = findUsedLevels(model, abundanceData);
    [abundanceRxns,  proteins_used, signifRxns] = selectGeneFromGPR(model, protein_id, protein_abun, parsedPR, minSum, gene_sig);
    
end
