root_data='/Users/aleboyer/ARNAUD/SCRIPPS/';
root_script='/Users/aleboyer/ARNAUD/SCRIPPS/EPSILON/';


Cruise_name='NISKINE'; % 
WW_name='epsifish2'; % 
deployement='d6';
epsifile = sprintf('epsi_%s.mat',deployement);
ctdfile  = sprintf('ctd_%s.mat',deployement);

% need seawater to use sw_bfrq
addpath /Users/aleboyer/ARNAUD/SCRIPPS/WireWalker/scripts/mixing_library/mixing_library/private1/seawater
addpath Toolbox/

ctdpath=sprintf('%s/%s/%s/%s/ctd/',root_data,Cruise_name,WW_name,deployement);
epsipath=sprintf('%s/%s/%s/%s/epsi/',root_data,Cruise_name,WW_name,deployement);
WWpath=sprintf('%s/%s/%s/%s/L1/',root_data,Cruise_name,WW_name,deployement);

name_ctd=[WW_name '_ctd_' deployement];

Epsi = load([epsipath epsifile]);
CTD  = load([ctdpath ctdfile]);
%%
[CTDProfile.up,CTDProfile.down,CTDProfile.dataup,CTDProfile.datadown] = ...
                                               get_upcast_sbe(CTD,5);
                                           
ind_ok=find(cellfun(@length,CTDProfile.up)>0);                                           
CTDProfile.up=CTDProfile.up(ind_ok);
CTDProfile.down=CTDProfile.down(ind_ok);
CTDProfile.datadown=CTDProfile.datadown(ind_ok);
CTDProfile.dataup=CTDProfile.dataup(ind_ok);

[EpsiProfile.up,EpsiProfile.down,EpsiProfile.dataup,EpsiProfile.datadown] =...
                                                     get_upcast_epsi(Epsi,CTD.time,CTDProfile.up,CTDProfile.down);


save([WWpath 'Profiles_' name_ctd],'CTDProfile','EpsiProfile','-v7.3');



