%% 
root_data='/Users/aleboyer/ARNAUD/SCRIPPS/SHEAR_PROBE/';
root_script='/Users/aleboyer/ARNAUD/SCRIPPS/WireWalker/WW_process_profile/';
Cruise_name='EPSIWW'; % 
WW_name='EPSI'; % 
deployement='d1';

%% add the needed toobox 
%addpath /Users/aleboyer/ARNAUD/SCRIPPS/WireWalker/scripts/mixing_library/mixing_library/private1/seawater
addpath toolbox/
addpath toolbox/CTD
addpath toolbox/FILTER
addpath toolbox/process



%% define path
WWpath=sprintf('%s/%s/WW/%s/%s/L1/',root_data,Cruise_name,WW_name,deployement);
epsipath=sprintf('%s/%s/WW/%s/%s/epsi/',root_data,Cruise_name,WW_name,deployement);
name_rbr=[WW_name '_rbr_' deployement];

%% 	get data 
load([ epsipath '/Profiles_' name_rbr],'EpsiProfile','RBRProfile')
load([WWpath 'EPSI_grid.mat'],'RBRgrid')
RBRprofiles=RBRProfile.dataup;
Profiles=EpsiProfile.dataup;

%% Parameters fixed by data structure
tscan     =  3;                                                            % length of 1 scan in second
FS        = round(1./nanmean(diff(Profiles{1}.time))/86400);               % sample rate channels
df        = 1/tscan;                                                       % number of samples per scan (1s) in channels
Sv        = 63.29285714;
Cs        = 683;

f=(df:df:FS/2)'; % frequency vector for spectra
Epsilon = struct([]);
chi     = struct([]);
timeaxis=zeros(1,length(Profiles));

q=3.2;

% add pressure from ctd to the epsi profile. This should be temporary until
% the addition of the pressure sensor on Epsi
%for i=1:length(Profiles)
for i=1:10
    %TODO correct the double(Sensor5) earlier in the process
    Profiles{i}.Sensor5=double(Profiles{i}.Sensor5); % TODO Sensor used for dev phase, it will be replace by conductivity
    Profiles{i}.P=interp1(RBRprofiles{i}.time,RBRprofiles{i}.P,Profiles{i}.rbrtime);
    Profiles{i}.T=interp1(RBRprofiles{i}.time,RBRprofiles{i}.T,Profiles{i}.rbrtime);
    Profiles{i}.S=interp1(RBRprofiles{i}.time,RBRprofiles{i}.S,Profiles{i}.rbrtime);
    [Epsilon{i}]=calc_eps_epsi(Profiles{i},tscan,f,Cs,Sv,i);
    [chi{i}]=calc_chi_epsi(Profiles{i},tscan,f,i);
    
    %Epsilon{i}.timeaxis=nanmean(Profiles{i}.rbrtime);
end

i=1;j=10;
for j=10:50
    figure
    [ktest,Ptest]=batchelor(Epsilon{i}.epsilon(j),chi{i}.chi(j),...
                            Epsilon{i}.kvis{j},chi{i}.ktemp{j},...
        q);
    
    loglog(Epsilon{i}.k{j},Epsilon{i}.Pshear{j},'b')
    hold on
    loglog(chi{i}.k{j},chi{i}.Ptgradk{j},'r')
    loglog(Epsilon{i}.kpan{j},Epsilon{i}.Ppan{j},'m')
    loglog(ktest,Ptest,'k')
    hold off
    legend('u_z','\phi_{TG}','Panchev','Batchelor','location','northwest')
    title(['Profile ' int2str(i) ' Scan=' int2str(j) ', pr=' ...
        num2str(Epsilon{i}.pr(j)) 'm']);
    text(3e2,1e-7,['\epsilon' sprintf('=%1.3e',Epsilon{i}.epsilon(j))]);
    text(3e2,1e-6,[ '\chi' sprintf('=%1.3e',chi{i}.chi(j))]);
    xlabel('k (cpm)','fontsize',15)
    ylabel('\phi(k)','fontsize',15)
    set(gca,'fontsize',15)
    fig=gcf;
    fig.PaperPosition = [0 0 8 8];
    print('-dpng',sprintf('../FIGURE/EPS_CHI/Uz_PHItg_Panch_Batch_%i.png',j))
    close all

end
            
%save([epsipath 'EpsiProfile.mat'],'Epsilon','chi','Profiles')

Map_pr=cellfun(@(x) (x.pr),Epsilon,'un',0);
zaxis=min([Map_pr{:}]):.5:max([Map_pr{:}]);
Map_epsilon=cellfun(@(x) interp1(x.pr,x.epsilon,zaxis),Epsilon,'un',0);
Map_epsilon=cell2mat(Map_epsilon.');
RBRgrid.Epsilon=interp1(zaxis,Map_epsilon.',RBRgrid.z);
save([WWpath WW_name '_grid.mat'],'RBRgrid')

load([WWpath WW_name '_grid.mat'],'RBRgrid')
level=min(RBRgrid.rho(:)):.01:max(RBRgrid.rho(:));
L=length(level);

t1=datenum('21-Jun-2017 00:00:00');
t2=datenum('21-Jun-2017 01:00:00');

close all
figure;
subplot('Position',[.1 .1 .6 .8])
%colormap('parula')
colormap('redbluecmap')
pcolor(RBRgrid.time,RBRgrid.z,log10(RBRgrid.Epsilon));shading interp;axis ij
hold on
contour(RBRgrid.time,RBRgrid.z,RBRgrid.rho,[level(1:10:L-200) level(L-200:L)],'k')
hold off
caxis([-9,-6])
set(gca,'XTickLabelRotation',45)
set(gca,'XTick',t1:t2-t1:RBRgrid.time(end))
set(gca,'XTickLabel',datestr(t1:t2-t1:RBRgrid.time(end),'HH'))
%datetick
%xlim([RBRgrid.time(1) RBRgrid.time(110)])
cax=colorbar;
xlabel(['Start date :' datestr(RBRgrid.time(1),'mm-dd-yyyy HH:MM:SS')],'fontsize',15)
set(gca,'fontsize',15)
ylabel(cax,'log_{10}(\epsilon) / W.kg^{-1}','fontsize',20)
ylabel('Depth (m)','fontsize',20)
title('EPSI on WireWalker','fontsize',20)
xlim([RBRgrid.time(1) RBRgrid.time(110)])
ylim([2 48])

subplot('Position',[.82 .1 .15 .8])
semilogx(RBRgrid.Epsilon(:,1),RBRgrid.z,'linewidth',2)
hold on
semilogx(RBRgrid.Epsilon(:,2),RBRgrid.z,'linewidth',2)
semilogx(RBRgrid.Epsilon(:,3),RBRgrid.z,'linewidth',2)
legend(datestr(RBRgrid.time(1),'HH:MM:SS'),...
       datestr(RBRgrid.time(2),'HH:MM:SS'),...
       datestr(RBRgrid.time(3),'HH:MM:SS'))
set(gca,'XTick',[1e-10 1e-8 1e-6])   
set(gca,'fontsize',15)
axis ij
ylim([2 48])
xlim([1e-9 1e-6])
title('Successive Profiles','fontsize',15)
xlabel('\epsilon / W.kg^{-1}','fontsize',20)


fig=gcf;
fig.PaperPosition = [0 0 15 10];
fig.PaperOrientation='Portrait';
print(sprintf('../FIGURE/EPSILON/%s_EpsiMap1.png',name_rbr),'-dpng2')


close all
figure;
colormap('jet')
pcolor(RBRgrid.time,RBRgrid.z,log10(RBRgrid.Epsilon));shading interp;axis ij
hold on
contour(RBRgrid.time,RBRgrid.z,RBRgrid.rho,[level(1:5:L-200) level(L-200:L)],'k')
hold off
caxis([-9,-6])
set(gca,'XTickLabelRotation',45)
set(gca,'XTick',RBRgrid.time(1:5:end))
datetick
xlim([RBRgrid.time(280) RBRgrid.time(370)])
cax=colorbar;
xlabel(['Start date :' datestr(RBRgrid.time(280),'mm-dd-yyyy')],'fontsize',15)
set(gca,'fontsize',15)
ylabel(cax,'log_{10}(\epsilon)','fontsize',20)
ylabel('Depth (m)','fontsize',20)
%xlim([RBRgrid.time(1) RBRgrid.time(end)])
ylim([2 48])

fig=gcf;
fig.PaperPosition = [0 0 15 10];
fig.PaperOrientation='Portrait';
print(sprintf('../FIGURE/EPSILON/%s_EpsiMap2.png',name_rbr),'-dpng2')



%% plot chi
Map_pr=cellfun(@(x) (x.pr),chi,'un',0);
zaxis=min([Map_pr{:}]):.5:max([Map_pr{:}]);
Map_chi=cellfun(@(x) interp1(x.pr,x.chi,zaxis),chi,'un',0);
Map_chi=cell2mat(Map_chi.');
RBRgrid.chi=interp1(zaxis,Map_chi.',RBRgrid.z);
%save([WWpath WW_name '_grid.mat'],'RBRgrid')


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








