function model = harmonizeReconECnumbers(model)

%% model EC numbers Harmonization
model.rxnECNumbers = regexprep(model.rxnECNumbers,";$","");
model.rxnECNumbers = strrep(model.rxnECNumbers,";"," & ");
model.rxnECNumbers = strrep(model.rxnECNumbers,","," & ");
model.rxnECNumbers = strrep(model.rxnECNumbers,"or","|");
model.rxnECNumbers = strrep(model.rxnECNumbers,"EC:","");
model.rxnECNumbers = strrep(model.rxnECNumbers,"TCDB:","");


end