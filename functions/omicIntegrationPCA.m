function result = omicIntegrationPCA(abundanceTable, expressionTable)

% Make a single table
omicDataTable = [abundanceTable, expressionTable];
omicData = table2array(omicDataTable);

% Set -1 values to NaN (This is the meaning in cobratoolbox)
omicData(omicData == -1) = NaN;
NaN_mask = any(~isnan(omicData), 2);

% Normalize omic data
omicDataNormalized = normalize(omicData(NaN_mask, :));

% Get the principal components
[coeff,pca_scores,latent,tsquared,explained,mu] = pca(omicDataNormalized);
reducedDimension = coeff(:,1:2); % Take the first two dimentions

% Get the sthe integrated data
omicDataNormalized(isnan(omicDataNormalized)) = 0;
reducedFeatureMatrix = (omicDataNormalized * reducedDimension) * latent(1:2) / sum(latent(1:2));
% Return NaN to -1 again
result = ones(size(omicData,1), 1) * -1;
result(NaN_mask) = reducedFeatureMatrix;

end