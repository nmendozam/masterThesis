function model = makePrules(model, gECNumbers)
%makePrules This function takes model.rxnECNumbers and makes a list with
%the protein rules (model.Prules) and a list of unique EC numbers (model.ECNumbers)
%
% USAGE:
%   model = makePrules(model)
%
% INPUT:
%   model:       cobra model structure
%
% OUTPUT:
%   model:   cobra model with aditional properties Prules and ECNumbers
%
% AUTHORS: Nicolas Mendoza-Mejia, Apr 2020

if ~exist('gECNumbers','var')
    gECNumbersFlag = 0;
    gECNumbers = model.rxnECNumbers;
else
    gECNumbersFlag = 1;
end

ecPattern = '\<(?!EC:|^\>)([0-9A-Z]+.(([0-9A-Z]|-)+.)*([0-9A-Z]|-)+)';


if gECNumbersFlag
    parsedECNumbers = regexprep(gECNumbers, ecPattern, 'x($1)');
    [Prules, ECNumbers] = makePrulesWithGrules(model.rules, parsedECNumbers);
else
    parsedECNumbers = regexp(gECNumbers, ecPattern, 'match');
    [Prules, ECNumbers] = makePrulesWithRxnECNumbers(model, parsedECNumbers);
end

model.Prules = Prules;
model.ECNumbers = ECNumbers;

end

function [index, ECList, size] = makeECList(inECNumber, ECList, size)

match = strcmp(inECNumber, ECList);

if isempty(match) || sum(match) == 0 
    size = size + 1;
    index = size;
    % Save in the list of unique ECNumbers
    ECList{size} = inECNumber;
else
    index = find(match);
end

end

function [Prules, ECNumbers] = makePrulesWithRxnECNumbers(model, rxnECNumbers)
nRxns = numel(model.rxns);
Prules = cell (nRxns, 1);
ECNumbers = {};
sizeEcNumbers = 0;

for i=1:nRxns

    nECNum = numel(rxnECNumbers{i});
    
    % This is done to avoid adding non-indexed things from rxnECNumbers
    if (nECNum > 0)
        Prules{i} = model.rxnECNumbers{i};
    end
    
    for j=1:nECNum
        [ECNumberIndex, ECNumbers, sizeEcNumbers] = makeECList(rxnECNumbers{i}{j}, ECNumbers, sizeEcNumbers);
        
        % Modify the protein rules
        pattern = regexptranslate('escape',rxnECNumbers{i}{j});
        replacement = ['x(' num2str(ECNumberIndex) ')'];
        Prules{i} = regexprep(Prules{i}, pattern, replacement, 'once');
    end
   
end

end

function [Prules, ECNumbers] = makePrulesWithGrules(gRules, gECNumbers)
Prules = gRules;
ECNumbers = {};
sizeEcNumbers = 0;

geneInexes = regexp(gRules,'(?<=x\()([0-9]+)(?=\))','match');

for i=1:numel(geneInexes)
    indexes = unique(str2double(geneInexes{i}));
    for j=1:numel(indexes)
        index = indexes(j);
        pattern = ['x(' num2str(index) ')'];
        
        if isempty(gECNumbers{index})
            % We need to delete the entry and one of the logic characters
            % around it
            regExp = regexptranslate('escape', pattern);
            Prules(i) = regexprep(Prules(i), [regExp ' *[&|] *'], '');
            Prules(i) = regexprep(Prules(i), ['( *[&|] *)?' regExp], '');
        else
            replacement = gECNumbers{index};
            match = regexp(gECNumbers{index},'(?<=x\().*?(?=\))','match');
            for k=1:numel(match)
                [ECNumberIndex, ECNumbers, sizeEcNumbers] = makeECList(match{k}, ECNumbers, sizeEcNumbers);
                pat = ['x(' ECNumbers{ECNumberIndex} ')'];
                rep = ['x(' num2str(ECNumberIndex) ')'];
                replacement = strrep(replacement, pat, rep);
            end
            
            Prules(i) = strrep(Prules(i), pattern, replacement);
        end
    end
end

Prules = regexprep(Prules,'((\||&) *\( *\) *)*$',''); % Remove emty parenthesis at the end
Prules = regexprep(Prules,' *\( *\) *(\||&)?',''); % Remove emty parenthesis at the start and in the midle

end


