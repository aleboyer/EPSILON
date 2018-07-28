 function [epsiup,epsidown,dataup,datadown] = EPSI_getcastepsi(data,CTDtime,up,down)

%  input:
%  output:
%
%  Created by Arnaud Le Boyer on 7/28/18.
%  Copyright © 2018 Arnaud Le Boyer. All rights reserved.

data=rmfield(data,'index');
data=rmfield(data,'sderror');
data=rmfield(data,'chsum1');
data=rmfield(data,'alti');
data=rmfield(data,'chsumepsi');

 
epsidown=cellfun( @(x) find(data.epsitime>=CTDtime(x(1)) & data.epsitime<=CTDtime(x(end))),down,'un',0);
epsiup=cellfun( @(x) find(data.epsitime>=CTDtime(x(1)) & data.epsitime<=CTDtime(x(end))),up,'un',0);

dataup=cellfun(@(x) structfun(@(y) y(x),data,'un',0),epsiup,'un',0);
datadown=cellfun(@(x) structfun(@(y) y(x),data,'un',0),epsidown,'un',0);

end
