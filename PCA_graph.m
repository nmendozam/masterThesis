% This script plots the PCA scatter plot for the different metabolic
% scenarios

%% Load abundance and expression data
abundanceTable_FileName = 'out/abundanceTable.mat';
expressionTable_FileName = 'out/expressionTable.mat';
load(abundanceTable_FileName, 'abundanceTable');
load(expressionTable_FileName, 'expressionTable');

abundanceTable_pa = abundanceTable(:, 1:6); % Only take palmitate
expressionTable_pa = expressionTable(:, 1:6); % Only take palmitate

% Make a single table
omicDataTable = [abundanceTable_pa, expressionTable_pa];
omicData = table2array(omicDataTable);

% Set -1 values to NaN (This is the meaning in cobratoolbox)
omicData(omicData == -1) = NaN;
NaN_mask = any(~isnan(omicData), 2);

% Normalize omic data
omicDataNormalized = normalize(omicData(NaN_mask, :));
writematrix(omicDataNormalized, "out/omicDataNormalized.xlsx")

% Get the principal components
[coeff,pca_scores,latent,tsquared,explained,mu] = pca(omicDataNormalized);

scatter(pca_scores(:,1),pca_scores(:,2),"filled")
axis equal
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
hold on

%% ##########################
%% Graph tibolone
%% ##########################
abundanceTable_tib = abundanceTable(:, 7:12); % Only take tibolone
expressionTable_tib = expressionTable(:, 1:6); % Only take palmitate

% Make a single table
omicDataTable = [abundanceTable_tib, expressionTable_tib];
omicData = table2array(omicDataTable);

% Set -1 values to NaN (This is the meaning in cobratoolbox)
omicData(omicData == -1) = NaN;
NaN_mask = any(~isnan(omicData), 2);

% Normalize omic data
omicDataNormalized = normalize(omicData(NaN_mask, :));

% Get the principal components
[coeff,pca_scores,latent,tsquared,explained,mu] = pca(omicDataNormalized);

scatter(pca_scores(:,1),pca_scores(:,2),"filled", "MarkerFaceColor", 'r')

%% ##########################
%% Graph controls
%% ##########################

abundanceTable_ctl = abundanceTable(:, 13:); % Only take control
expressionTable_ctl = expressionTable(:, 1:6); % Only take palmitate

% Make a single table
omicDataTable = [abundanceTable_ctl, expressionTable_ctl];
omicData = table2array(omicDataTable);

% Set -1 values to NaN (This is the meaning in cobratoolbox)
omicData(omicData == -1) = NaN;
NaN_mask = any(~isnan(omicData), 2);

% Normalize omic data
omicDataNormalized = normalize(omicData(NaN_mask, :));

% Get the principal components
[coeff,pca_scores,latent,tsquared,explained,mu] = pca(omicDataNormalized);

scatter(pca_scores(:,1),pca_scores(:,2), "x", "MarkerFaceColor", 'g')
legend('Palmitate', 'Tibolone', 'Control')