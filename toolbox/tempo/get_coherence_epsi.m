%root_data = '/Volumes/Ahua/data_archive/WaveChasers-DataArchive/EPSI_SPROUL/';%EPSI_SPROUL/EPSI_SPROUL/EPSI/d3/';
root_data = '/Users/aleboyer/ARNAUD/SCRIPPS/SHEAR_PROBE/';


root_script='/Users/aleboyer/ARNAUD/SCRIPPS/WireWalker/WW_process_profile/';
Cruise_name='EPSI_SPROUL'; % 
WW_name='EPSI'; % 
deployement='d2';

% need seawater to use sw_bfrq
addpath Toolbox/

ctdpath=sprintf('%s/%s/%s/%s/ctd/',root_data,Cruise_name,WW_name,deployement);
WWpath=sprintf('%s/%s/%s/%s/L1/',root_data,Cruise_name,WW_name,deployement);
name_ctd=[WW_name '_ctd_' deployement];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0=clock;
close all

% 	get data 
load([WWpath 'Profiles_' name_ctd '.mat'],'EpsiProfile','CTDProfile')
CTDprofiles=CTDProfile.datadown{1};
Profiles=EpsiProfile.datadown;


% add pressure from ctd to the epsi profile. 
nb_profiles=length(Profiles);
%for i=1:nb_profiles
for i=1
    
    %TODO correct the double(Sensor5) earlier in the process
    Profiles{i}.Sensor5=double(Profiles{i}.Sensor5);
   % Profiles{i}.P=interp1(CTDprofiles{i}.time,CTDprofiles{i}.P,Profiles{i}.time);
    Profiles{i}.P=interp1(CTDprofiles.time,CTDprofiles.P,Profiles{i}.time);
end

%% Parameters fixed by data structure
WWvel     = .71;                                                            %50 cm/s
tscan     =  3;
scan      = WWvel*tscan;                                                   % depth of a scan 3 second
FS        = round(1./nanmean(diff(Profiles{1}.time)));                     % sample rate channels
df        = 1/tscan;                                                       % number of samples per scan (1s) in channels


%% Parameters for processing of microstructure data (epsilon and chi)
f=(df:df:FS/2)'; % frequency vector for spectra
Pco12=zeros(nb_profiles,length(f));
Pco13=Pco12;
Pco14=Pco12;
timeaxis=Pco12;

fig=figure('Position',[100,100,1600,2000]);
v = VideoWriter('Coherence_EPSISPROUL.avi');
v.FrameRate=10;
open(v)

for j=1:nb_profiles
%for j=1:10
    timeaxis(j)=nanmean(Profiles{j}.time);
    [Pco12(j,:),Pco13(j,:),Pco14(j,:)]=calc_coherence_epsi(Profiles{j},tscan,f,j);
    frame=getframe(fig);
    writeVideo(v,frame)
    clf;
end
close(v)

close all
pcolor(timeaxis(:,1),f,Pco13.');shading flat
datetick
xlabel([datestr(timeaxis(1,1),'mm-dd-yyyy') ' start'],'fontsize',20)
ylabel('Hz','fontsize',20)
title('Coherence Shear-Ax','fontsize',20)
xlim(timeaxis([1 end],1))
caxis([0 0.5])
cax=colorbar;
ylabel(cax,'Coherence','fontsize',20)
set(gca,'fontsize',20)

print('../FIGURE/COHERENCE/Pco13.png','-dpng2')


