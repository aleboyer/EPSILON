root_data='../';
root_script='/Users/aleboyer/ARNAUD/SCRIPPS/EPSILON/';
epsifile = 'epsi_rbrsync_EPSI.mat';
ctdfile  = 'Profiles_WW_rbr_d2.mat';

Cruise_name='NISKINE'; % 
vehicle_name='WW'; % 
deployement='d2';

% need seawater to use sw_bfrq
addpath /Users/aleboyer/ARNAUD/SCRIPPS/WireWalker/scripts/mixing_library/mixing_library/private1/seawater
addpath Toolbox/

ctdpath=sprintf('%s/%s/%s/%s/rbr/',root_data,Cruise_name,vehicle_name,deployement);
epsipath=sprintf('%s/%s/%s/%s/epsi/',root_data,Cruise_name,vehicle_name,deployement);
WWpath=sprintf('%s/%s/%s/%s/L1/',root_data,Cruise_name,vehicle_name,deployement);

name_ctd=[vehicle_name '_ctd_' deployement];


Epsi = load([epsipath epsifile]);
CTD  = load([ctdpath ctdfile]);

CTDProfiles=CTD.RBRprofiles;
Epsi.Sensor5 = Epsi.Sensor1*nan;
EpsiProfiles  = get_cast_epsiWW(Epsi,CTDProfiles);

save([WWpath 'Profiles_' name_ctd],'CTDProfile','EpsiProfile','-v7.3');



