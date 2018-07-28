
root_data='../';
root_script='/Users/aleboyer/ARNAUD/SCRIPPS/EPSILON/';
epsifile = 'epsi_d3.mat';
ctdfile  = 'ctd_d3.mat';


Cruise_name='pregranite'; % 
WW_name='epsifish'; % 
deployement='d3';



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
                                               get_upcast_sbe(CTD,3);
Epsi.Sensor5=Epsi.Sensor1*nan;
indok=cellfun(@length,CTDProfile.up);
indok=find(indok>0);

CTDProfile.up=CTDProfile.up(indok);
CTDProfile.down=CTDProfile.down(indok);
CTDProfile.dataup=CTDProfile.dataup(indok);
CTDProfile.datadown=CTDProfile.datadown(indok);

[EpsiProfile.up,EpsiProfile.down,EpsiProfile.dataup,EpsiProfile.datadown] =...
                                                     get_upcast_epsi(Epsi,CTD.time,CTDProfile.up,CTDProfile.down);

save([WWpath 'Profiles_' name_ctd],'CTDProfile','EpsiProfile','-v7.3');



