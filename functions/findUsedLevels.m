function [id, expr, sig] = findUsedLevels(id, exprData, printLevel)
% Returns vectors of gene identifiers and corresponding gene expression
% levels or protein abundance for each gene present in the model 
% ('model.genes' or 'model.ECNumbers').
%
% USAGE:
%
%    [gene_id, gene_expr] = findUsedGenesLevels(model, exprData)
%    [gene_id, gene_expr, gene_sig] = findUsedGenesLevels(model, exprData)
%
% INPUTS:
%
%   model:               input model (COBRA model structure)
%
%   exprData:            mRNA expression data structure
%       .id                cell array containing Gene or protein IDs in the
%                            same format as model.genes
%       .value               Vector containing corresponding expression value (FPKM)
%       .sig:                [optional field] Vector containing significance values of
%                            expression corresponding to expression values in exprData.value (ex. p-values)
%
% OPTIONAL INPUTS:
%    printLevel:         Printlevel for output (default 0);
%
% OUTPUTS:
%
%   id:             vector of gene identifiers present in the model
%                        that are associated with expression data
%
%   expr:           vector of expression values associated to each
%                        'id'
%
% OPTIONAL OUTPUTS:
%   sig:             vector of significance values associated to each
%                        'id'
%
%   
% Authors: - S. Opdam & A. Richelle May 2017
%       - Chaitra Sarathy, Oct 2019, add significance value as optional input
%       - Nicolas Mendoza-Mejia, Apr 2020, support for protein abundances

if ~exist('printLevel','var')
    printLevel = 0;
end

if isfield(exprData, 'sig') 
    exprSigFlag = 1; 
else
    exprSigFlag = 0;
end 

expr=[];
sig=[];

for i = 1:numel(id)
    
    cur_ID = id{i};
    % Find matches
    if contains(cur_ID, ".-")
        cur_ID = regexprep(cur_ID, '(\.-)+', '.');
        found = cellfun(@(x) contains(x, cur_ID), exprData.id, 'UniformOutput', false);
    else
        found = cellfun(@(x) strcmp(x, cur_ID), exprData.id, 'UniformOutput', false);
    end
    
    dataID = find(cellfun(@(c) any(c(:)), found));
    
    if isempty (dataID)
        expr(i)=-1;        
    elseif length(dataID)==1
        expr(i)=exprData.value(dataID);
        if exprSigFlag == 1 
            sig(i) = exprData.sig(dataID);
        end 
    elseif length(dataID)>1    	
        if printLevel > 0
            disp(['Double for ',num2str(cur_ID)])
        end
        expr(i) = sum(exprData.value(dataID));
        if exprSigFlag == 1 
            sig(i) = sum(exprData.sig(dataID));
        end 
    end    
end
           
end