clear
%% 
if exist('/Volumes/aleboyer/ARNAUD/SCRIPPS/PLUMEX18/','dir')
    root_data='/Volumes/aleboyer/ARNAUD/SCRIPPS/PLUMEX18/';
    root_script='/Volumes/aleboyer/ARNAUD/SCRIPPS/EPSILON/';
else
    root_data='/Volumes/DataDrive/NISKINE/';
    root_script='/Users/aleboyer/ARNAUD/SCRIPPS/EPSILON/';
end


Cruise_name='NISKINE'; % 
WW_name='epsifish1'; % 
deployement='d5';

%% add the needed toobox 
%addpath /Users/aleboyer/ARNAUD/SCRIPPS/WireWalker/scripts/mixing_library/mixing_library/private1/seawater
addpath toolbox/
addpath toolbox/CTD
addpath toolbox/FILTER
addpath toolbox/process
addpath toolbox/seawater2/
addpath toolbox/PLOTS/


%% define path
WWpath=sprintf('%s/%s/%s/%s/L1/',root_data,Cruise_name,WW_name,deployement);
epsipath=sprintf('%s/%s/%s/%s/epsi/',root_data,Cruise_name,WW_name,deployement);
ctdpath=sprintf('%s/%s/%s/%s/ctd/',root_data,Cruise_name,WW_name,deployement);
name_ctd=[WW_name '_ctd_' deployement];
name_epsi=[WW_name '_epsi_' deployement];

%% 	get data 
load([ WWpath '/Profiles_' name_ctd],'EpsiProfile','CTDProfile')
CTD_Profiles=CTDProfile.datadown;
EPSI_Profiles=EpsiProfile.datadown;

%% Parameters fixed by data structure 
tscan     =  3;                                                            % length of 1 scan in second
%FS        = round(1./nanmean(diff(Profiles{1}.time)));                     % sample rate channels
FS        = 325;                     % sample rate channels
df        = 1/tscan;                                                       % number of samples per scan (1s) in channels
% shear 103,104 = s2,s1
Sv        = [42.14,43.36]; % epsifish1
%Sv        = [53.00,48.16]; % epsifish2
q=3.2;

f=(df:df:FS/2)'; % frequency vector for spectra
Epsilon = struct([]);
chi     = struct([]);

% add pressure from ctd to the epsi profile. This should be ter mporary until
% the addition of the pressure sensor on Epsi

% epsi2 d2 23rd cast is bad 
% epsi2 d3 15th cast is bad 
indok=cellfun(@(x) length(x.time),EPSI_Profiles);
indok=find(indok>10000);
EPSI_Profiles=EPSI_Profiles(indok);
CTD_Profiles=CTD_Profiles(indok);
timeaxis=zeros(1,length(EPSI_Profiles));

%for i=1:length(EPSI_Profiles)
for i=17:length(EPSI_Profiles)
    fprintf('Profile %i over %i \n',i,length(EPSI_Profiles))
    %TODO correct the double(Sensor5) earlier in the process
    [ctdtime,IA,IB]=unique(CTD_Profiles{i}.time);
    P=CTD_Profiles{i}.P(IA);
    S=real(CTD_Profiles{i}.S(IA));
    T=CTD_Profiles{i}.T(IA);
    EPSI_Profiles{i}.Sensor1=EPSI_Profiles{i}.Sensor2;
    EPSI_Profiles{i}.Sensor5=EPSI_Profiles{i}.Sensor1*0+1;
    EPSI_Profiles{i}.Sensor7=EPSI_Profiles{i}.Sensor1*0+1;
    EPSI_Profiles{i}.P=interp1(ctdtime,P,EPSI_Profiles{i}.time);
    EPSI_Profiles{i}.T=interp1(ctdtime,T,EPSI_Profiles{i}.time);
    EPSI_Profiles{i}.S=interp1(ctdtime,S,EPSI_Profiles{i}.time);
    MS{i}=calc_turbulence_epsi_FastCTD(EPSI_Profiles{i},tscan,f,Sv);
    %Epsilon{i}.timeaxis=nanmean(Profiles{i}.rbrtime);
end
save([WWpath 'Turbulence_Profiles.mat'],'MS','EPSI_Profiles')

%load([WWpath 'Turbulence_Profiles.mat'],'MS','EPSI_Profiles')
indok=cellfun(@isempty,MS);
Epsilon_class=calc_binned_epsi(MS(~indok));
Chi_class=calc_binned_chi(MS(~indok));

[F1,F2]=plot_binned_epsilon(Epsilon_class,name_epsi);
print(F1,[WWpath deployement '_binned_epsilon1_t3s.png'],'-dpng2')
print(F2,[WWpath deployement '_binned_epsilon2_t3s.png'],'-dpng2')

[F1,F2]=plot_binned_chi(Chi_class,name_epsi,1:40);
print(F1,[WWpath deployement '_binned_chi22_c_t3s.png'],'-dpng2')
print(F2,[WWpath deployement '_binned_chi21_c_t3s.png'],'-dpng2')


            
MSempty=cellfun(@isempty,MS);
Map_pr=cellfun(@(x) (x.pr),MS(~MSempty),'un',0);
zaxis=min([Map_pr{:}]):.5:max([Map_pr{:}]);
Map_epsilon=cellfun(@(x) interp1(x.pr,x.epsilon(:,1),zaxis),MS(~MSempty),'un',0);
Map_time=cell2mat(cellfun(@(x) mean(x.time),MS(~MSempty),'un',0));

Map_epsilon=cell2mat(Map_epsilon.');

name_ctd=['ctd_' deployement];
if strcmp(ctdpath,'/Volumes/DataDrive/NISKINE//NISKINE/epsifish2/d5/ctd/')
    load([ctdpath 'ctd_sd_d5.mat'])
else
    load([ctdpath name_ctd])
end
[tt,pp]=meshgrid(Map_time,zaxis);
F=scatteredInterpolant(time.',P.',sig.');
sig2D=F(tt,pp);
sig2D=sig2D-1000;
level=min(sig2D(:)):.01:max(sig2D(:));
L=length(level);
L1=round(L-L/4);
close all
figure;
colormap('jet')
pcolor(Map_time,zaxis,log10(real(Map_epsilon.')));shading flat;axis ij
colorbar
hold on
contour(Map_time,zaxis,sig2D,[level(1:5:L-L1) level(L-L1:L)],'k')
hold off
caxis([-9,-4])
set(gca,'XTickLabelRotation',45)
set(gca,'XTick',Map_time(1:5:end))
datetick
cax=colorbar;
xlabel(['Start date :' datestr(Map_time(1),'mm-dd-yyyy')],'fontsize',15)
set(gca,'fontsize',15)
ylabel(cax,'log_{10}(\epsilon)','fontsize',20)
ylabel('Depth (m)','fontsize',20)
%xlim([RBRgrid.time(1) RBRgrid.time(end)])
ylim([0 700])

fig=gcf;
fig.PaperPosition = [0 0 15 10];
fig.PaperOrientation='Portrait';
print([WWpath deployement '_EpsiMap1.png'],'-dpng2')


            
MSempty=cellfun(@isempty,MS);
Map_pr=cellfun(@(x) (x.pr),MS(~MSempty),'un',0);
zaxis=min([Map_pr{:}]):.5:max([Map_pr{:}]);
Map_epsilon=cellfun(@(x) interp1(x.pr,x.epsilon(:,2),zaxis),MS(~MSempty),'un',0);
Map_time=cell2mat(cellfun(@(x) mean(x.time),MS(~MSempty),'un',0));

Map_epsilon=cell2mat(Map_epsilon.');



figure;
colormap('jet')
pcolor(Map_time,zaxis,log10(real(Map_epsilon.')));shading flat;axis ij
colorbar
hold on
contour(Map_time,zaxis,sig2D,[level(1:5:L-L1) level(L-L1:L)],'k')
hold off
caxis([-9,-4])
set(gca,'XTickLabelRotation',45)
set(gca,'XTick',Map_time(1:5:end))
datetick
cax=colorbar;
xlabel(['Start date :' datestr(Map_time(1),'mm-dd-yyyy')],'fontsize',15)
set(gca,'fontsize',15)
ylabel(cax,'log_{10}(\epsilon)','fontsize',20)
ylabel('Depth (m)','fontsize',20)
%xlim([RBRgrid.time(1) RBRgrid.time(end)])
ylim([0 700])

fig=gcf;
fig.PaperPosition = [0 0 15 10];
fig.PaperOrientation='Portrait';
print([WWpath deployement '_EpsiMap2.png'],'-dpng2')

