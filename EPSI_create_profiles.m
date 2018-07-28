function EPSI_create_profiles(Meta_Data)

%  input:
%  output:
%
%  Created by Arnaud Le Boyer on 7/28/18.
%  Copyright © 2018 Arnaud Le Boyer. All rights reserved.


CTDpath=Meta_Data.CTDpath;
Epsipath=Meta_Data.Epsipath;
L1path=Meta_Data.L1path;

CTD=load([CTDpath 'ctd_' Meta_Data.deployement '.mat'],'ctdtime','T','P','S','sig');
EPSI=load([Epsipath 'epsi_' Meta_Data.deployement '.mat']);


[CTDProfile.up,CTDProfile.down,CTDProfile.dataup,CTDProfile.datadown] = ...
                                               EPSI_getcastctd(CTD,20);
                                           
[EpsiProfile.up,EpsiProfile.down,EpsiProfile.dataup,EpsiProfile.datadown] =...
                                                     EPSI_getcastepsi(EPSI,CTD.ctdtime,CTDProfile.up,CTDProfile.down);

                                                 
fprintf('Saving data in %sProfiles_%s.mat\n',L1path,Meta_Data.deployement)

save([L1path 'Profiles_' Meta_Data.deployement '.mat'],'CTDProfile','EpsiProfile','-v7.3');



