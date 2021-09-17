function [modelExchanges, selExc] = findExchangeReactions(model)

modelexchanges1 = strmatch('Ex_',model.rxns);
modelexchanges4 = strmatch('EX_',model.rxns);
modelexchanges2 = strmatch('DM_',model.rxns);
modelexchanges3 = strmatch('sink_',model.rxns);
selExc = (find( full((sum(abs(model.S)==1,1) ==1) & (sum(model.S~=0) == 1))))';
modelExchanges = unique([modelexchanges1;modelexchanges2;modelexchanges3;modelexchanges4;selExc]);

end