function expresionTable = EnsemblToEntrez(expresionTable, EntrezFileName)

% EntrezFileName is the path to the file downloaded from  https://www.ensembl.org/biomart/martview

f = fopen(EntrezFileName);
C = textscan(f, '%s %s', 'HeaderLines', 1);
fclose(f);
EnsemblIDs = regexprep(expresionTable.(1),"\.[0-9]*","");
EnrezIDs = cell(numel(EnsemblIDs), 1);
for i=1:numel(EnsemblIDs)
    match = strcmp(C{1},EnsemblIDs(i));
    if sum(match) > 0
        EnrezIDs{i} = C{2}{match};
    end
end

expresionTable = addvars(expresionTable, EnrezIDs, 'After', 1);

end