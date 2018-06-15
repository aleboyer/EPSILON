 function cast = get_cast_epsiWW(Epsi,CTDProfiles)
indcast=cellfun( @(x) find(Epsi.EPSItime>=x.time(1) & Epsi.EPSItime<=x.time(end)),CTDProfiles,'un',0);
cast=cellfun(@(x) structfun(@(y) y(x),Epsi,'un',0),indcast,'un',0);
end
