root_data='/Volumes/Ahua/data_archive/WaveChasers-DataArchive/EPSI_SPROUL/'
Cruise_name='EPSI_SPROUL'; % 
WW_name='EPSI'; % 
deployement='d3';

% need seawater to use sw_bfrq
addpath Toolbox/

ctdpath=sprintf('%s/%s/%s/%s/ctd/',root_data,Cruise_name,WW_name,deployement);
WWpath=sprintf('%s/%s/%s/%s/L1/',root_data,Cruise_name,WW_name,deployement);
epsipath=sprintf('%s/%s/%s/%s/epsi/',root_data,Cruise_name,WW_name,deployement);
name_ctd=[WW_name '_ctd_' deployement];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0=clock;
close all

% 	get data 
load([WWpath 'Profiles_' name_ctd '.mat'],'EpsiProfile')
load([WWpath 'Profiles_' name_ctd '.mat'],'CTDProfile')
CTDprofiles=CTDProfile.datadown;
Profiles=EpsiProfile.datadown;


% add pressure from ctd to the epsi profile. This should be temporary until
% the addition of the pressure sensor on Epsi
nb_profiles=length(Profiles);
for i=1:nb_profiles
    Profiles{i}.P=interp1(CTDprofiles{i}.time,CTDprofiles{i}.P,Profiles{i}.time);
    Profiles{i}.T=interp1(CTDprofiles{i}.time,CTDprofiles{i}.T,Profiles{i}.time);
    Profiles{i}.time=datetime(Profiles{i}.time,'ConvertFrom','posixtime');
end

%% Parameters fixed by data structure
tscan     =  15;
%FS        = round(1./nanmean(diff(Profiles{1}.time))/86400);               % sample rate channels
FS        = 320;               % sample rate channels

FPO7_noise(Profiles{1},tscan,FS)
fig=gcf;
fig.PaperPosition = [0 0 15 8];
print('-dpng2','../FIGURE/EPSILON/EPSISPROUL_d2_P1_FPO7_noise.png')

Shear_noise(Profiles{1},tscan,FS)
fig=gcf;
fig.PaperPosition = [0 0 15 8];
print('-dpng2','../FIGURE/EPSILON/EPSISPROUL_d2_P1_Shear_noise.png')

Accel_noise(Profiles{1},tscan,FS)
fig=gcf;
fig.PaperPosition = [0 0 15 8];
print('-dpng2','../FIGURE/EPSILON/EPSISPROUL_d2_P1_Accel_noise.png')

% Profile2=structfun(@(x) x(5000:end),Profiles{2},'un',0);
% FPO7_noise(Profile2,tscan,FS)
% Shear_noise(Profile2,tscan,FS)
% 
% Profile3=structfun(@(x) x(5000:end),Profiles{3},'un',0);
% FPO7_noise(Profile3,tscan,FS)
% Shear_noise(Profile3,tscan,FS)
% 
% Profile4=structfun(@(x) x(5000:end),Profiles{4},'un',0);
% FPO7_noise(Profile4,tscan,FS)
% Shear_noise(Profile4,tscan,FS)

