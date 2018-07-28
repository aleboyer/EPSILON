function EPSI_batchprocess_epsifish(Meta_Data)

L1path=Meta_Data.L1path;

%% add the needed toobox  move that to create Meta file
%addpath /Users/aleboyer/ARNAUD/SCRIPPS/WireWalker/scripts/mixing_library/mixing_library/private1/seawater
addpath(genpath('toolbox/'))


if ~exist([L1path 'Turbulence_Profiles.mat'],'file')
    %% 	get data
    load([L1path 'Profiles_' Meta_Data.deployement '.mat'],'CTDProfile','EpsiProfile');
    CTD_Profiles=CTDProfile.datadown;
    EPSI_Profiles=EpsiProfile.datadown;
    
    %% Parameters fixed by data structure
    tscan     =  3;                                                            % length of 1 scan in second
    FS        = str2double(Meta_Data.Firmware.sampling_frequency(1:3));                     % sample rate channels
    df        = 1/tscan;                                                       % number of samples per scan (1s) in channels
    
    f=(df:df:FS/2)'; % frequency vector for spectra
    MS = struct([]);
    
    % add pressure from ctd to the epsi profile. This should be ter mporary until
    % the addition of the pressure sensor on Epsi
    for i=1:length(EPSI_Profiles)
        EPSI_Profiles{i}.P=interp1(CTD_Profiles{i}.ctdtime,CTD_Profiles{i}.P,EPSI_Profiles{i}.epsitime);
        EPSI_Profiles{i}.T=interp1(CTD_Profiles{i}.ctdtime,CTD_Profiles{i}.T,EPSI_Profiles{i}.epsitime);
        EPSI_Profiles{i}.S=interp1(CTD_Profiles{i}.ctdtime,CTD_Profiles{i}.S,EPSI_Profiles{i}.epsitime);
        MS{i}=calc_turbulence_epsi_FastCTD(EPSI_Profiles{i},tscan,f,45,Meta_Data);
    end
    save([L1path 'Turbulence_Profiles.mat'],'MS')
else
    load([L1path 'Turbulence_Profiles.mat'],'MS')
end

% compite binned epsilon for all profiles
Epsilon_class=calc_binned_epsi(MS);
Chi_class=calc_binned_chi(MS);

% plot binned epsilon for all profiles
[F1,F2]=plot_binned_epsilon(Epsilon_class,[Meta_Data.mission '-' Meta_Data.deployement]);
print(F1,[L1path Meta_Data.deployement '_binned_epsilon1_t3s.png'],'-dpng2')
print(F2,[L1path Meta_Data.deployement '_binned_epsilon2_t3s.png'],'-dpng2')

%[F1,F2]=plot_binned_chi(Chi_class,Meta_Data,[Meta_Data.mission '-' Meta_Data.deployement]);
%print(F1,[L1path Meta_Data.deployement '_binned_chi22_c_t3s.png'],'-dpng2')
%print(F2,[L1path Meta_Data.deployement '_binned_chi21_c_t3s.png'],'-dpng2')


MSempty=cellfun(@isempty,MS);
Map_pr=cellfun(@(x) (x.pr),MS(~MSempty),'un',0);
zaxis=min([Map_pr{:}]):.5:max([Map_pr{:}]);
Map_epsilon2=cellfun(@(x) interp1(x.pr,x.epsilon(:,2),zaxis),MS(~MSempty),'un',0);
Map_epsilon1=cellfun(@(x) interp1(x.pr,x.epsilon(:,1),zaxis),MS(~MSempty),'un',0);
Map_chi1=cellfun(@(x) interp1(x.pr,x.chi(:,1),zaxis),MS(~MSempty),'un',0);
Map_chi2=cellfun(@(x) interp1(x.pr,x.chi(:,2),zaxis),MS(~MSempty),'un',0);


Map_time=cell2mat(cellfun(@(x) mean(x.time),MS(~MSempty),'un',0));

Map_epsilon1=cell2mat(Map_epsilon1.');
Map_epsilon2=cell2mat(Map_epsilon2.');
Map_chi1=cell2mat(Map_chi1.');
Map_chi2=cell2mat(Map_chi2.');
save([L1path  'Turbulence_grid.mat'],'Map_epsilon1','Map_epsilon2','Map_chi1','Map_chi2','Map_time','zaxis')

close all

% epsilon 1 
figure;
colormap('jet')
pcolor(Map_time,zaxis,log10(real(Map_epsilon1.')));shading flat;axis ij
colorbar
caxis([-9,-6])
set(gca,'XTickLabelRotation',45)
datetick
cax=colorbar;
xlabel(['Start date :' datestr(Map_time(1),'mm-dd-yyyy')],'fontsize',15)
set(gca,'fontsize',15)
ylabel(cax,'log_{10}(\epsilon)','fontsize',20)
ylabel('Depth (m)','fontsize',20)

fig=gcf;
fig.PaperPosition = [0 0 15 10];
fig.PaperOrientation='Portrait';
print(sprintf('%sEpsiMap1.png',L1path),'-dpng2')

% epsilon 2
figure;
colormap('jet')
pcolor(Map_time,zaxis,log10(real(Map_epsilon2.')));shading flat;axis ij
colorbar
caxis([-9,-6])
set(gca,'XTickLabelRotation',45)
datetick
cax=colorbar;
xlabel(['Start date :' datestr(Map_time(1),'mm-dd-yyyy')],'fontsize',15)
set(gca,'fontsize',15)
ylabel(cax,'log_{10}(\epsilon)','fontsize',20)
ylabel('Depth (m)','fontsize',20)

fig=gcf;
fig.PaperPosition = [0 0 15 10];
fig.PaperOrientation='Portrait';
print(sprintf('%sEpsiMap2.png',L1path),'-dpng2')

% chi 1 
figure;
colormap('jet')
pcolor(Map_time,zaxis,log10(real(Map_chi1.')));shading flat;axis ij
colorbar
caxis([-9,-6])
set(gca,'XTickLabelRotation',45)
datetick
cax=colorbar;
xlabel(['Start date :' datestr(Map_time(1),'mm-dd-yyyy')],'fontsize',15)
set(gca,'fontsize',15)
ylabel(cax,'log_{10}(\chi)','fontsize',20)
ylabel('Depth (m)','fontsize',20)

fig=gcf;
fig.PaperPosition = [0 0 15 10];
fig.PaperOrientation='Portrait';
print(sprintf('%sEpsichi1.png',L1path),'-dpng2')

% epsilon 1 
figure;
colormap('jet')
pcolor(Map_time,zaxis,log10(real(Map_chi2.')));shading flat;axis ij
colorbar
caxis([-9,-6])
set(gca,'XTickLabelRotation',45)
datetick
cax=colorbar;
xlabel(['Start date :' datestr(Map_time(1),'mm-dd-yyyy')],'fontsize',15)
set(gca,'fontsize',15)
ylabel(cax,'log_{10}(\chi)','fontsize',20)
ylabel('Depth (m)','fontsize',20)

fig=gcf;
fig.PaperPosition = [0 0 15 10];
fig.PaperOrientation='Portrait';
print(sprintf('%sEpsichi2.png',L1path),'-dpng2')





