root_data='/Users/aleboyer/ARNAUD/SCRIPPS/PLUMEX18/';
root_script='/Users/aleboyer/ARNAUD/SCRIPPS/EPSILON/';
Cruise_name='Plumex_Feb2018'; % 
WW_name='EPSIFISH'; % 
deployement='mar3';
epsifile = 'Plumex18_matthewtest_epsi.mat';
ctdfile  = 'Plumex18_matthewtest_ctd.mat';


% need seawater to use sw_bfrq
addpath /Users/aleboyer/ARNAUD/SCRIPPS/WireWalker/scripts/mixing_library/mixing_library/private1/seawater
addpath Toolbox/

ctdpath=sprintf('%s/%s/%s/%s/ctd/',root_data,Cruise_name,WW_name,deployement);
epsipath=sprintf('%s/%s/%s/%s/epsi/',root_data,Cruise_name,WW_name,deployement);
WWpath=sprintf('%s/%s/%s/%s/L1/',root_data,Cruise_name,WW_name,deployement);

name_ctd=[WW_name '_ctd_' deployement];

Epsi = load([epsipath epsifile]);
CTD  = load([ctdpath ctdfile]);

[CTDProfile.up,CTDProfile.down,CTDProfile.dataup,CTDProfile.datadown] = ...
                                               get_upcast_sbe(CTD);
Epsi.Sensor5=Epsi.Sensor1*nan;
[EpsiProfile.up,EpsiProfile.down,EpsiProfile.dataup,EpsiProfile.datadown] =...
                                                     get_upcast_epsi(Epsi,CTD.ctdtime,CTDProfile.up,CTDProfile.down);


save([WWpath 'Profiles_' name_ctd],'CTDProfile','EpsiProfile','-v7.3');



