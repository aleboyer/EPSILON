%% 
root_data='../EPSIWW/';
root_script='/Users/aleboyer/ARNAUD/SCRIPPS/EPSILON/';
epsifile = 'Profiles_EPSI_rbr_d1.mat';
ctdfile  = 'Profiles_EPSI_rbr_d1.mat';
Cruise_name='WW'; % 
vehicle_name='EPSI'; % 
deployement='d1';


%% add the needed toobox 
%addpath /Users/aleboyer/ARNAUD/SCRIPPS/WireWalker/scripts/mixing_library/mixing_library/private1/seawater
addpath toolbox/
addpath toolbox/CTD
addpath toolbox/FILTER
addpath toolbox/process
addpath toolbox/PLOTS


%% define path
WWpath=sprintf('%s/%s/%s/%s/L1/',root_data,Cruise_name,vehicle_name,deployement);
epsipath=sprintf('%s/%s/%s/%s/epsi/',root_data,Cruise_name,vehicle_name,deployement);
name_ctd=[vehicle_name '_ctd_' deployement];

%% 	get data 
%load('/Users/aleboyer/ARNAUD/SCRIPPS/EPSIWW/WW/EPSI/d1/rbr/Profiles_EPSI_rbr_d1.mat')
load('/Users/aleboyer/ARNAUD/SCRIPPS/EPSIWW/WW/EPSI/d1/epsi/Profiles_EPSI_rbr_d1.mat')
%load([ WWpath '/Profiles_' name_ctd],'EpsiProfiles','RBRProfiles')
%load([WWpath 'WW_grid.mat'],'RBRgrid')
EpsiProfiles=EpsiProfile.dataup;
CTDProfiles=RBRProfile.dataup;
LEpsiProfile=cellfun(@(x) length(x.rbrtime),EpsiProfiles);
indok=find(LEpsiProfile>10);
EpsiProfiles=EpsiProfiles(indok);
CTDProfiles=CTDProfiles(indok);

%% Parameters fixed by data structure
tscan     =  3;                                                            % length of 1 scan in second
FS        = round(1./nanmean(diff(EpsiProfiles{1}.rbrtime))/86400);               % sample rate channels
df        = 1/tscan;                                                       % number of samples per scan (1s) in channels
Sv        = [70.63 26.43]; % SV for sh1 and sh2 from Meta_SP1810_d1.dat

f=(df:df:FS/2)'; % frequency vector for spectra
Epsilon = struct([]);
chi     = struct([]);
timeaxis=zeros(1,length(EpsiProfiles));

q=3.2;


% add pressure from ctd to the epsi profile. This should be ter mporary until
% the addition of the pressure sensor on Epsi
for i=1:length(EpsiProfiles)
    %TODO correct the double(Sensor5) earlier in the process
    EpsiProfiles{i}.P=interp1(CTDProfiles{i}.time,CTDProfiles{i}.P,EpsiProfiles{i}.rbrtime);
    EpsiProfiles{i}.T=interp1(CTDProfiles{i}.time,CTDProfiles{i}.T,EpsiProfiles{i}.rbrtime);
    S=sw_salt(CTDProfiles{i}.C*10./sw_c3515,CTDProfiles{i}.T,CTDProfiles{i}.P);
    EpsiProfiles{i}.S=interp1(CTDProfiles{i}.time,S,EpsiProfiles{i}.rbrtime);
    MS{i}=calc_turbulence_epsiWW(EpsiProfiles{i},tscan,f,Sv);
    %Epsilon{i}.timeaxis=nanmean(Profiles{i}.rbrtime);
end
save([WWpath 'Turbulence_Profiles.mat'],'MS','EpsiProfiles','-v7.3')

load([WWpath 'Turbulence_Profiles.mat'],'MS','EPSI_Profiles')
Epsilon_class=calc_binned_epsi(MS);
Chi_class=calc_binned_chi(MS);

[F1,F2]=plot_binned_epsilon(Epsilon_class,'NISKINE WW d2 ');
print(F1,[WWpath deployement '_binned_epsilon1_t3s.png'],'-dpng2')
print(F2,[WWpath deployement '_binned_epsilon2_t3s.png'],'-dpng2')

[F1,F2]=plot_binned_chi(Chi_class,'NISKINE WW d2',20:40);
print(F1,[WWpath deployement '_binned_chi22_c_t5s.png'],'-dpng2')
print(F2,[WWpath deployement '_binned_chi21_c_t5s.png'],'-dpng2')

indok=cellfun(@isempty,MS);
MS=MS(~indok);
load([WWpath 'WW_grid.mat'],'RBRgrid')

Map_pr=cellfun(@(x) (x.pr),MS,'un',0);
zaxis=min([Map_pr{:}]):.5:max([Map_pr{:}]);
Map_epsilon1=cellfun(@(x) interp1(x.pr,x.epsilon(:,1),zaxis),MS,'un',0);
Map_epsilon1=cell2mat(Map_epsilon1.');
RBRgrid.Epsilon1=interp1(zaxis,Map_epsilon1.',RBRgrid.z);
Map_epsilon2=cellfun(@(x) interp1(x.pr,x.epsilon(:,2),zaxis),MS,'un',0);
Map_epsilon2=cell2mat(Map_epsilon2.');
RBRgrid.Epsilon2=interp1(zaxis,Map_epsilon2.',RBRgrid.z);

save([WWpath vehicle_name '_grid.mat'],'RBRgrid')

load([WWpath vehicle_name '_grid.mat'],'RBRgrid')
level=nanmin(RBRgrid.rho(:)):.01:nanmax(RBRgrid.rho(:));
L=length(level);

t1=RBRgrid.time(indok(1));
t2=t1+.5;
D=150;
close all
figure;
subplot('Position',[.1 .1 .6 .8])
%colormap('parula')
colormap('redbluecmap')
pcolor(RBRgrid.time,RBRgrid.z,log10(1e-2.*RBRgrid.Epsilon1));shading flat;axis ij
hold on
contour(RBRgrid.time,RBRgrid.z,RBRgrid.rho,[level(1:10:L-D) level(L-D:2:L)],'k')
hold off
caxis([-9,-6])
set(gca,'XTickLabelRotation',45)
set(gca,'XTick',t1:.1:t2)
set(gca,'XTickLabel',datestr(t1:.1:t2,'HH'))
%datetick
cax=colorbar;
set(gca,'fontsize',15)
ylabel(cax,'log_{10}(\epsilon) / W.kg^{-1}','fontsize',20)
ylabel('Depth (m)','fontsize',20)
xlabel(['Start:' datestr(t1)],'fontsize',20)
title('Epsilon - shear 2 - WW - La Jolla','fontsize',20)
ylim([2 50])
xlim([t1 t2])

eventlog=find(RBRgrid.time>=datenum('21-Jun-2017 07:00:00'),1,'first');
eventlog1=find(RBRgrid.time>=datenum('21-Jun-2017 08:00:00'),1,'first');
eventlog2=find(RBRgrid.time>=datenum('21-Jun-2017 09:00:00'),1,'first');
subplot('Position',[.82 .1 .15 .8])
semilogx(1e-2.*RBRgrid.Epsilon1(:,eventlog),RBRgrid.z,'linewidth',2)
hold on
semilogx(1e-2.*RBRgrid.Epsilon1(:,eventlog1),RBRgrid.z,'linewidth',2)
semilogx(1e-2.*RBRgrid.Epsilon1(:,eventlog2),RBRgrid.z,'linewidth',2)
legend(datestr(RBRgrid.time(eventlog),'HH:MM:SS'),...
       datestr(RBRgrid.time(eventlog1),'HH:MM:SS'),...
       datestr(RBRgrid.time(eventlog2),'HH:MM:SS'))
set(gca,'XTick',[1e-10 1e-8 1e-6])   
set(gca,'fontsize',15)
axis ij
ylim([2 50])
xlim([1e-10 1e-4])
title('Profiles','fontsize',15)
xlabel('\epsilon / W.kg^{-1}','fontsize',20)


fig=gcf;
fig.PaperPosition = [0 0 15 10];
fig.PaperOrientation='Portrait';
print(sprintf('%s/%s_EpsiMap1.png',WWpath,name_ctd),'-dpng2')

close all
figure;
subplot('Position',[.1 .1 .6 .8])
%colormap('parula')
colormap('redbluecmap')
pcolor(RBRgrid.time(indok),RBRgrid.z,log10(RBRgrid.Epsilon2));shading interp;axis ij
hold on
contour(RBRgrid.time(indok),RBRgrid.z,RBRgrid.rho(:,indok),[level(1:2:L-D) level(L-D:5:L)],'k')
hold off
caxis([-8,-4.5])
set(gca,'XTickLabelRotation',45)
set(gca,'XTick',t1:t2-t1:RBRgrid.time(end))
set(gca,'XTickLabel',datestr(t1:t2-t1:RBRgrid.time(end),'HH'))
datetick
cax=colorbar;
set(gca,'fontsize',15)
ylabel(cax,'log_{10}(\epsilon) / W.kg^{-1}','fontsize',20)
ylabel('Depth (m)','fontsize',20)
xlabel(['Start:' datestr(t1)],'fontsize',20)
title('Epsilon - shear 2 - WW - SP1810','fontsize',20)
xlim([t1 t2])
ylim([2 300])

subplot('Position',[.82 .1 .15 .8])
semilogx(RBRgrid.Epsilon2(:,1),RBRgrid.z,'linewidth',2)
hold on
semilogx(RBRgrid.Epsilon2(:,2),RBRgrid.z,'linewidth',2)
semilogx(RBRgrid.Epsilon2(:,3),RBRgrid.z,'linewidth',2)
legend(datestr(RBRgrid.time(indok(1)),'HH:MM:SS'),...
       datestr(RBRgrid.time(indok(2)),'HH:MM:SS'),...
       datestr(RBRgrid.time(indok(3)),'HH:MM:SS'))
set(gca,'XTick',[1e-10 1e-8 1e-6])   
set(gca,'fontsize',15)
axis ij
ylim([2 300])
xlim([1e-8 1e-4])
title('Successive Profiles','fontsize',15)
xlabel('\epsilon / W.kg^{-1}','fontsize',20)

fig=gcf;
fig.PaperPosition = [0 0 15 10];
fig.PaperOrientation='Portrait';
print(sprintf('%s/%s_EpsiMap2.png',WWpath,name_ctd),'-dpng2')



%% plot chi
Map_pr=cellfun(@(x) (x.pr),chi,'un',0);
zaxis=min([Map_pr{:}]):.5:max([Map_pr{:}]);
Map_chi=cellfun(@(x) interp1(x.pr,x.chi,zaxis),chi,'un',0);
Map_chi=cell2mat(Map_chi.');
RBRgrid.chi=interp1(zaxis,Map_chi.',RBRgrid.z);
%save([WWpath vehicle_name '_grid.mat'],'RBRgrid')


level=min(RBRgrid.rho(:)):.01:max(RBRgrid.rho(:));
L=length(level);

close all
figure;
colormap('jet')
pcolor(RBRgrid.time,RBRgrid.z,log10(-RBRgrid.chi));shading interp;axis ij
hold on
contour(RBRgrid.time,RBRgrid.z,RBRgrid.rho,[level(1:5:L-200) level(L-200:L)],'k')
hold off
caxis([-18,-8])
set(gca,'XTickLabelRotation',45)
set(gca,'XTick',RBRgrid.time(1:5:end))
datetick
%xlim([RBRgrid.time(1) RBRgrid.time(110)])
cax=colorbar;
xlabel(['Start date :' datestr(RBRgrid.time(1),'mm-dd-yyyy')],'fontsize',15)
set(gca,'fontsize',15)
ylabel(cax,'log_{10}(-\chi)','fontsize',20)
ylabel('Depth (m)','fontsize',20)
%xlim([RBRgrid.time(1) RBRgrid.time(end)])
ylim([2 48])

fig=gcf;
fig.PaperPosition = [0 0 15 10];
fig.PaperOrientation='Portrait';
print(sprintf('../FIGURE/CHI/%s_chiMap1.png',name_rbr),'-dpng2')


close all
figure;
colormap('jet')
pcolor(RBRgrid.time,RBRgrid.z,log10(-RBRgrid.chi));shading interp;axis ij
hold on
contour(RBRgrid.time,RBRgrid.z,RBRgrid.rho,[level(1:5:L-200) level(L-200:L)],'k')
hold off
caxis([-18,-6])
set(gca,'XTickLabelRotation',45)
set(gca,'XTick',RBRgrid.time(1:5:end))
datetick
xlim([RBRgrid.time(280) RBRgrid.time(370)])
cax=colorbar;
xlabel(['Start date :' datestr(RBRgrid.time(280),'mm-dd-yyyy')],'fontsize',15)
set(gca,'fontsize',15)
ylabel(cax,'log_{10}(\chi)','fontsize',20)
ylabel('Depth (m)','fontsize',20)
%xlim([RBRgrid.time(1) RBRgrid.time(end)])
ylim([2 48])

fig=gcf;
fig.PaperPosition = [0 0 15 10];
fig.PaperOrientation='Portrait';
print(sprintf('../FIGURE/CHI/%s_EpsiMap2.png',name_rbr),'-dpng2')








