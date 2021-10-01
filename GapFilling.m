close all force; clear variables; clc

%% Startup the COBRA Toolbox
initCobraToolbox(false) % Don't update the toolbox
changeCobraSolver ('gurobi', 'all', 1); % For large models
clc

%% Load model
filename = 'out/NoLeakAstrociteModel.mat';
load(filename, 'modelNoLeak');
model = modelNoLeak; clear modelNoLeak

%% Identify deadend metabolites
% Detect dead ends
outputMets = detectDeadEnds(model);
% Print the names
DeadEnds = model.mets(outputMets);

% Identify associated reactions
[rxnList, rxnFormulaList] = findRxnsFromMets(model, DeadEnds);
% look at the bounds
model.lb(find(ismember(model.rxns,rxnList)))
model.ub(find(ismember(model.rxns,rxnList)))

% Identify gaps with downstream metabolites
[allGaps, rootGaps, downstreamGaps] = gapFind(model, 'true');

%% Identify bloqued reactions
BlockedReactions = findBlockedReaction(model);


%% statistics
cnt = 1;
Stats{cnt,1} = 'Model name';cnt = cnt+1;
Stats{cnt,1} = 'Size S (original model)';cnt = cnt+1;
Stats{cnt,1} = 'Number of compartments';cnt = cnt+1;
Stats{cnt,1} = 'List of compartments';cnt = cnt+1;
Stats{cnt,1} = 'Number of blocked reactions';cnt = cnt+1;
Stats{cnt,1} = 'Number of solvable blocked reactions';cnt = cnt+1;
Stats{cnt,1} = 'Size S (flux consistent)';cnt = cnt+1;
Stats{cnt,1} = 'Size SUX (including solvable blocked reactions)';cnt = cnt+1;
Stats{cnt,1} = 'Number of added reactions (all)';cnt = cnt+1;
Stats{cnt,1} = 'Number of added metabolic reactions ';cnt = cnt+1;
Stats{cnt,1} = 'Number of added transport reactions ';cnt = cnt+1;
Stats{cnt,1} = 'Number of added exchange reactions ';cnt = cnt+1;
Stats{cnt,1} = 'Time preprocessing';cnt = cnt+1;
Stats{cnt,1} = 'Time fastGapFill';cnt = cnt+1;

col = 1;
RxnList={};
cnt = 1;
% Remove constrains from exchange reactions
EX = strmatch('EX_',model.rxns);
model.lb(EX)=-1000;
model.ub(EX)=1000;
clear EX
% Basic statistics
Stats{cnt,col+1} = filename;cnt = cnt+1;
[a,b] = size(model.S);
Stats{cnt,col+1} = strcat(num2str(a),'x',num2str(b));cnt = cnt+1;
% List of compartments
[tok,rem] = strtok(model.mets,'\[');
rem = unique(rem);
Stats{cnt,col+1} = num2str(length(rem));cnt = cnt+1;
Rem = rem{1};
for j = 2:length(rem)
Rem = strcat(Rem,',',rem{j});
end
Stats{cnt,col+1} = Rem;cnt = cnt+1;
clear Rem tok rem;

%% Pre-Gapfill
% Parameters
weights.MetabolicRxns = 0.1; % Kegg metabolic reactions
weights.ExchangeRxns = 0.5; % Exchange reactions
weights.TransportRxns = 10; % Transport reactions
% Prepare for gapFill (Adds reference from database and transport
% reactions)
tic; [consistModel,consistMatricesSUX,BlockedRxns] = prepareFastGapFill(model);
tpre=toc;
% Add statistics to the table
Stats{cnt,col+1} = num2str(length(BlockedRxns.allRxns));cnt = cnt+1;
Stats{cnt,col+1} = num2str(length(BlockedRxns.solvableRxns));cnt = cnt+1;
[a,b] = size(consistModel.S);
Stats{cnt,col+1} = strcat(num2str(a),'x',num2str(b));cnt = cnt+1;
[a,b] = size(consistMatricesSUX.S);
Stats{cnt,col+1} = strcat(num2str(a),'x',num2str(b));cnt = cnt+1;

%% Gapfill
epsilon = 1e-4;
tic; [AddedRxns] = fastGapFill(consistMatricesSUX,epsilon, weights);
tgap=toc;
% Add statistics to the table
Stats{cnt,col+1} = num2str(length(AddedRxns.rxns));cnt = cnt+1;

%% Post-Gapfill
IdentifyPW = 0;
[AddedRxnsExtended] = postProcessGapFillSolutions(AddedRxns,model,BlockedRxns);
clear AddedRxns;
% Add statistics to the table
Stats{cnt,col+1} = num2str(AddedRxnsExtended.Stats.metabolicSol);cnt = cnt+1;
Stats{cnt,col+1} = num2str(AddedRxnsExtended.Stats.transportSol);cnt = cnt+1;
Stats{cnt,col+1} = num2str(AddedRxnsExtended.Stats.exchangeSol);cnt = cnt+1;
Stats{cnt,col+1} = num2str(tpre);cnt = cnt+1;
Stats{cnt,col+1} = num2str(tgap);cnt = cnt+1;
clear a b
% Generate reaction list
col = col + 1;
RxnList{1,col}=filename;RxnList(2:length(AddedRxnsExtended.rxns)+1,col) = AddedRxnsExtended.rxns; col = col + 1;
RxnList{1,col}=filename;RxnList(2:length(AddedRxnsExtended.rxns)+1,col) = AddedRxnsExtended.rxnFormula; col = col + 1;
RxnList{1,col}=filename;RxnList(2:length(AddedRxnsExtended.rxns)+1,col) = AddedRxnsExtended.subSystem; col = col + 1;
clear AddedRxnsExtended tgap tpre j BlockedRxns i cnt consistM*

%% Analysis
RxnList
% Print deadend metabolites
metid = detectDeadEnds(model,false);
horzcat(model.mets(metid),model.metNames(metid))


