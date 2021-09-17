function parsedPR = PRparser(model, getCNFSets)
% Maps the PR rules of the model to a specified format that is used by
% the model extraction methods 
%
% USAGE:
%   parsedPR = PRparser(model)
%
% INPUT:
%   model:       cobra model structure
%
% OPTIONAL INPUT:
%    getCNFSets:    whether to get the CNF sets (true) or DNF sets (false).
%                   DNF sets represent functional enzyme complexes, while
%                   CNF sets represent the possible subunits of a complex.
%                   (default: false , i.e. DNF sets)
%
% OUTPUT:
%   parsedPR:   cell matrix containing parsed PR rule
%
% AUTHORS: Nicolas Mendoza-Mejia, Apr 2020

if ~exist('getCNFSets','var')
    getCNFSets = false;
end

parsedPR = {};
fp = FormulaParser();
for i = 1:numel(model.rxns)
    if ~isempty(model.Prules{i})
        head = fp.parseFormula(model.Prules{i});
        currentSets = head.getFunctionalGeneSets(model.ECNumbers,getCNFSets)';
        parsedPR{i}=currentSets;
    else
        parsedPR{i}={''};
    end
end
parsedPR=parsedPR';
end