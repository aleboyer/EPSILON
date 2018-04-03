root_data = '/Volumes/Ahua/data_archive/WaveChasers-DataArchive/EPSI_SPROUL/';%EPSI_SPROUL/EPSI_SPROUL/EPSI/d3/';

Cruise_name='EPSI_SPROUL'; % 
WW_name='EPSI'; % 
deployement='d3';

%% add the needed toobox 
%addpath /Users/aleboyer/ARNAUD/SCRIPPS/WireWalker/scripts/mixing_library/mixing_library/private1/seawater
addpath toolbox/
addpath toolbox/CTD
addpath toolbox/FILTER
addpath toolbox/process

%% define path
WWpath=sprintf('%s/%s/%s/%s/L1/',root_data,Cruise_name,WW_name,deployement);
epsipath=sprintf('%s/%s/%s/%s/epsi/',root_data,Cruise_name,WW_name,deployement);
name_ctd=[WW_name '_ctd_' deployement];


load([WWpath 'EpsiProfile.mat'],'Epsilon','Epsilon_mmp','Profiles')
Profile=Profiles{1};


Lscan     = length(Epsilon{1}.indscan{1});
Lscan_mmp = length(Epsilon_mmp{1}.indscan{1});
nbscan=length(Epsilon{1}.indscan);
nbscan_mmp=length(Epsilon_mmp{1}.indscan);
P1=zeros(nbscan,Lscan);
P2=zeros(nbscan,Lscan);
P3=zeros(nbscan_mmp,Lscan_mmp);
P4=zeros(nbscan_mmp,Lscan_mmp);
P11=zeros(nbscan,Lscan);
P22=zeros(nbscan,Lscan);
P33=zeros(nbscan_mmp,Lscan_mmp);
P44=zeros(nbscan_mmp,Lscan_mmp);
%% Parameters fixed by data structure
FS  = round(1./nanmean(diff(Profile.time)));                     % sample rate channels
FS_mmp  = round(1./nanmean(diff(Profile.time_mmpshear)));                     % sample rate channels
df  = Profile.time(Epsilon{1}.indscan{1}(end));
df  = df -Profile.time(Epsilon{1}.indscan{1}(1));
df=1./ceil(df);
f     = (df:df:FS/2)'; % frequency vector for spectra
f_mmp = (df:df:FS_mmp/2)'; % frequency vector for spectra
G=9.81;

Sv     = [48.66,50.54];
Sv_mmp =[1 1]*34.2;
H_MMP=load('toolbox/MMP/MMP4_shear_Htransfer.mat'); 
h_freqmmp=interp1(H_MMP.f,H_MMP.h_freq,f_mmp);
h_freq=get_filters_MADRE(f);

for j=1:nbscan
    speed=Epsilon{1}.w(j)/100;
    Pco12=zeros(1,length(f));
    Pco13=Pco12;
    Pco14=Pco12;
    timeaxis=Pco12;
    htotal = (Sv.*speed/(2*G)).^2 .* h_freq.shear .* haf_oakey(f,speed);        % should add epsi filter
    htotalmmp = (Sv_mmp.*speed/(2*G)).^2 .* h_freqmmp .* haf_oakey(f_mmp,speed);        % should add epsi filter
    
    window = hanning(Lscan);
    windowmmp = hanning(Lscan_mmp);
    wc2=1/mean(window.^2);% window correction factor
    wc2mmp=1/mean(windowmmp.^2);% window correction factor
    data1 = Profile.Sensor3(Epsilon{1}.indscan{j});
    data2 = Profile.Sensor4(Epsilon{1}.indscan{j});
    s1    = Profile.s1(Epsilon_mmp{1}.indscan{j});
    s2    = Profile.s2(Epsilon_mmp{1}.indscan{j});

    A1=window.*data1(:);
    A2=window.*data2(:);
    A3=windowmmp.*s1(:);
    A4=windowmmp.*s2(:);
    
    P1(j,:)=fft(A1);P2(j,:)=fft(A2);P3(j,:)=fft(A3);P4(j,:)=fft(A4);

    P11(j,:)=abs(fft(A1)).^2;
    P22(j,:)=abs(fft(A2)).^2;
    P33(j,:)=abs(fft(A3)).^2;
    P44(j,:)=abs(fft(A4)).^2;
    P11(j,:)=P11(j,:)./Lscan^2/df*wc2;
    P22(j,:)=P22(j,:)./Lscan^2/df*wc2;
    P33(j,:)=P33(j,:)./Lscan_mmp^2/df*wc2mmp;
    P44(j,:)=P44(j,:)./Lscan_mmp^2/df*wc2mmp;
    
end

Pi3=interp1([flipud(-f_mmp) ;f_mmp],P3.',[flipud(-f); f]).';
Pi4=interp1([flipud(-f_mmp); f_mmp],P4.',[flipud(-f); f]).';
Pi33=interp1([flipud(-f_mmp); f_mmp],P33.',[flipud(-f); f]).';
Pi44=interp1([flipud(-f_mmp); f_mmp],P44.',[flipud(-f); f]).';

P12=conj(P1).*P2./(Lscan)^2/df*wc2;
P13=conj(P1).*Pi3./(Lscan)^2/df*wc2;
P14=conj(P1).*Pi4./(Lscan)^2/df*wc2;

CoP12=P12(100:end-100,:).^2./(P11(100:end-100,:).*P22(100:end-100,:));
CoP13=P13(100:end-100,:).^2./(P11(100:end-100,:).*Pi33(100:end-100,:));
CoP14=P14(100:end-100,:).^2./(P11(100:end-100,:).*Pi44(100:end-100,:));

CoP12=abs(mean(CoP12,1));
CoP13=abs(mean(CoP13,1));
CoP14=abs(mean(CoP14,1));


Pa1=mean(P11(100:end-20,:),1);
Pa2=mean(P22(20:end-20,:),1);
Pa3=mean(P33(20:end-20,:),1);
Pa4=mean(P44(20:end-20,:),1);
% Pa1=2*Pa1(1:Lscan/2)./htotal(:,1).';
% Pa2=2*Pa2(1:Lscan/2)./htotal(:,2).';
% Pa3=2*Pa3(1:Lscan_mmp/2)./htotalmmp.';
% Pa4=2*Pa4(1:Lscan_mmp/2)./htotalmmp.';

Pa1=2*Pa1(1:Lscan/2)./htotal(:,1).';
Pa2=2*Pa2(1:Lscan/2)./htotal(:,2).';
Pa3=2*Pa3(1:Lscan_mmp/2)./htotalmmp.';
Pa4=2*Pa4(1:Lscan_mmp/2)./htotalmmp.';


Pco12=CoP12(1:Lscan/2);
Pco13=CoP13(1:Lscan/2);
Pco14=CoP14(1:Lscan/2);

Fn    = .5*FS;  % Nyquist frequency
FR    = 2.5;      % Full range
Nb_by = 24;     % number of bytes
%[(full range = 3V)/2^24 ]^2 / Fn
noise24= (FR/2^Nb_by)^2 /Fn;
Nb_by = 16;     % number of bytes
%[(full range = 3V)/2^24 ]^2 / Fn
noise16= (FR/2^Nb_by)^2 /Fn;



ax(1)=subplot(211);
loglog(f*1.19,Pa1/100,'k')
hold on
loglog(f_mmp*1.19,Pa3/100,'r')
loglog(f_mmp,f_mmp*0+noise16,'k--')
hold off

ax(2)=subplot(212);
semilogx(f*1.19^2,Pco13)
ylim([0,1])

linkaxes(ax,'x')
