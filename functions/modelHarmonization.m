function model = standarizeReconECnumbers(model)

% model reactions names Harmonization
model.rxns(ismember(model.rxns,'ATPM'))={'DM_atp_c_'};
model.rxns(ismember(model.rxns,'ATPhyd'))={'DM_atp_c_'};
model.rxns(ismember(model.rxns,'DM_atp(c)'))={'DM_atp_c_'};
model.rxns(ismember(model.rxns,'EX_biomass_reaction'))={'biomass_reaction'};
model.rxns(ismember(model.rxns,'EX_biomass_maintenance'))={'biomass_maintenance'};
model.rxns(ismember(model.rxns,'EX_biomass_maintenance_noTrTr'))={'biomass_maintenance_noTrTr'};
% Use of brackets
model.rxns = regexprep(model.rxns,'\(','\[');
model.rxns = regexprep(model.rxns,'\)','\]');
model.rxns = regexprep(model.rxns,'Ex_','EX_');
model.rxns = regexprep(model.rxns,'Sink_','sink_');
model.rxns = regexprep(model.rxns,'-','_'); 

end