 function [epsiup,epsidown,dataup,datadown] = get_upcast_epsi(data,CTDtime,up,down)

% CTDtime  is actually not time yet it is a epsi sample stamp

 
%epsidown=cellfun( @(x) find(data.EPSItime>=CTDtime(x(1)) & data.EPSItime<=CTDtime(x(end))),down,'un',0);
%epsiup=cellfun( @(x) find(data.EPSItime>=CTDtime(x(1)) & data.EPSItime<=CTDtime(x(end))),up,'un',0);
epsidown=cellfun( @(x) find(data.nbsample>=CTDtime(x(1)) & data.nbsample<=CTDtime(x(end))),down,'un',0);
epsiup=cellfun( @(x) find(data.nbsample>=CTDtime(x(1)) & data.nbsample<=CTDtime(x(end))),up,'un',0);

dataup=cellfun(@(x) structfun(@(y) y(x),data,'un',0),epsiup,'un',0);
datadown=cellfun(@(x) structfun(@(y) y(x),data,'un',0),epsidown,'un',0);

end
