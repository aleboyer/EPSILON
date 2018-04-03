root_data='/Users/aleboyer/ARNAUD/SCRIPPS/SHEAR_PROBE/';
root_script='/Users/aleboyer/ARNAUD/SCRIPPS/WireWalker/WW_process_profile/';
Cruise_name='EPSIWW'; % 
WW_name='EPSI'; % 
deployement='d1';

% need seawater to use sw_bfrq
addpath Toolbox/

rbrpath=sprintf('%s/%s/WW/%s/%s/rbr/',root_data,Cruise_name,WW_name,deployement);
WWpath=sprintf('%s/%s/WW/%s/%s/L1/',root_data,Cruise_name,WW_name,deployement);
name_rbr=[WW_name '_rbr_' deployement];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0=clock;
close all

% 	get data 
load(['../../EPSIWW/WW/EPSI/d1/epsi/' 'Profiles_' name_rbr '.mat'],'EpsiProfile')
load(['../../EPSIWW/WW/EPSI/d1/epsi/' 'Profiles_' name_rbr '.mat'],'RBRProfile')
RBRprofiles=RBRProfile.dataup;
Profiles=EpsiProfile.dataup;


% add pressure from ctd to the epsi profile. This should be temporary until
% the addition of the pressure sensor on Epsi
nb_profiles=length(Profiles);
for i=1:nb_profiles
    %TODO correct the double(Sensor5) earlier in the process
    Profiles{i}.Sensor5=double(Profiles{i}.Sensor5);
    Profiles{i}.P=interp1(RBRprofiles{i}.time,RBRprofiles{i}.P,Profiles{i}.rbrtime);
end

%% Parameters fixed by data structure
tscan     =  15;
%FS        = round(1./nanmean(diff(Profiles{1}.time))/86400);               % sample rate channels
FS        = 325;               % sample rate channels
FSrbr     = round(1./nanmean(diff(RBRprofiles{1}.time))/86400);               % sample rate channels
df        = 1/tscan;                                                       % number of samples per scan (1s) in channels
Sv        = 63.29285714;
Cs        = 683;


%% Parameters for processing of microstructure data (epsilon and chi)
f    = (df:df:FS/2)'; % frequency vector for spectra
frbr = (df:df:FSrbr/2)'; % frequency vector for spectra
%Pco12=zeros(nb_profiles,length(f));


%bit noise
Fn    = .5*FS;  % Nyquist frequency
Fnrbr = .5*FSrbr;  % Nyquist frequency
FR    = 3.3;      % Full range
Nb_by = 24;     % number of bytes
%[(full range = 3V)/2^24 ]^2 / f_N where f_N =200 =#
noise24= (FR/2^Nb_by)^2 /Fn;
Nb_by = 16;     % number of bytes
%[(full range = 3V)/2^24 ]^2 / f_N where f_N =200 =#
noise16= (FR/2^Nb_by)^2 /Fn;

%Johnson noise
k  = 1.38E-23;
T0  = 273.15;
R  = 664e3;
vJ0 = sqrt(4*k*T0*R*FS);
T20  = 273.15+20;
vJ20 = sqrt(4*k*T20*R*FS);
T50  = 273.15+50;
vJ50 = sqrt(4*k*T50*R*FS);



Ptot=0;
Ptotrbr=0;
for j=1:length(Profiles)
    nb_profile=j;
    Profile=Profiles{j};
    RBRProfile=RBRprofiles{j};
    
    dataP=max(Profile.Sensor1)-Profile.Sensor1;
    dataP=dataP./max(dataP);
    datarbrP=RBRProfile.T-min(RBRProfile.T);
    datarbrP=datarbrP./max(datarbrP);
    
    timeepsi=RBRProfile.time(1):1/325/86400:RBRProfile.time(1)+(length(dataP)/325-1/325)/86400;
    % for j=1:nb_profiles
    % %for j=1:10
    %     timeaxis(j)=nanmean(Profiles{j}.rbrtime);
    %     [Pco12(j,:),Pco13(j,:),Pco14(j,:)]=calc_coherence_epsi(Profiles{j},tscan,f,j);
    % end
    
    %% Length of the Profile
    T       = length(Profile.time);
    Trbr    = length(RBRProfile.time);
    df      = f(1);
    %% define number of scan in the profile
    Lscan     = tscan*2*f(end);
    Lscanrbr  = tscan*2*frbr(end);
    nbscan    = floor(T/Lscan);
    nbscanrbr = floor(Trbr/Lscanrbr);
    
    %% we compute spectra on scan with 50% overlap
    nbscan    = 2*nbscan-1;
    nbscanrbr = 2*nbscanrbr-1;
    
    %% define the index in the profile for each scan
    indscan = arrayfun(@(x) (1+floor(Lscan/2)*(x-1):1+floor(Lscan/2)*(x-1)+Lscan-1),1:nbscan,'un',0);
    indscanrbr = arrayfun(@(x) ...
                           (1+floor(Lscanrbr/2)*(x-1):1+floor(Lscanrbr/2)*(x-1)+Lscanrbr-1),...
                           1:nbscanrbr,'un',0);
    
    %% define fmax the cut off frequency after which the vibration of the
    window = hanning(Lscan);
    windowrbr = hanning(Lscanrbr);

    P=zeros(nbscan,Lscan);
    Prbr=zeros(nbscanrbr,Lscanrbr);
    
    gain_adc=1;
    sinc4_freq=f/2/f(end);
    hfilt = gain_adc.*(sinc(sinc4_freq)).^4;

    
    for i=1:nbscan
%        data = Profile.Sensor1(indscan{i}(1:Lscan));
%        datarbr = RBRProfile.T(indscanrbr{i}(1:Lscanrbr));
        data = dataP(indscan{i}(1:Lscan));
        datarbr = datarbrP(indscanrbr{i}(1:Lscanrbr));
        A=window.*data(:);
        B=windowrbr.*datarbr(:);
        P(i,:)=abs(fft(A)).^2;
        Prbr(i,:)=abs(fft(B)).^2;
    end
    nu=2*floor(Lscan/nbscan);
    err_low = nu/chi2inv(.05/2,nu);
    err_high = nu/chi2inv(1-.05/2,nu);
    
    nurbr=2*floor(Lscanrbr/nbscanrbr);
    err_lowrbr = nurbr/chi2inv(.05/2,nurbr);
    err_highrbr = nurbr/chi2inv(1-.05/2,nurbr);

    
    Pa1=mean(P,1)./Lscan^2/f(1);
    Pa1=2*Pa1(1:floor(Lscan/2))./hfilt.';

    Pa1rbr=mean(Prbr,1)./Lscanrbr^2/f(1);
    Pa1rbr=2*Pa1rbr(1:floor(Lscanrbr/2)).';
    
    Ptot=Ptot+Pa1;
    Ptotrbr=Ptotrbr+Pa1rbr;
    %% plot coherence
    ax(1)=subplot(311);
    ax(2)=subplot(312);
    ax(3)=subplot(313);
    time1=86400*(RBRProfile.time-RBRProfile.time(1));
    
    
    l1=plot(ax(1),time1,RBRProfile.P,'linewidth',2);
    axis(ax(1),'ij')
    xlim(ax(1),time1([1 end]))
    ylim(ax(1),[min(Profile.P), max(Profile.P)])
    ylabel(ax(1),'Presurre (db)','fontsize',15)
    set(ax(1),'fontsize',15)
    title(ax(1),['Profile ' num2str(nb_profile) '-' datestr(nanmean(Profile.rbrtime))])
    
    %channels=detrend(Profile.Sensor1.','constant').';
    channels=Profile.Sensor1;
    
    % if = 
    
    l2rbr=plot(ax(2),86400*(RBRProfile.time-RBRProfile.time(1)),datarbrP,'r','linewidth',2);
    hold(ax(2),'on')
%    l2=plot(ax(2),1.11*(Profile.time-Profile.time(1))*86400,dataP,'b','linewidth',2);
    l2=plot(ax(2),1.142*(timeepsi-timeepsi(1))*86400,dataP,'b','linewidth',2);
    hold(ax(2),'off')
    xlim(ax(2),time1([1 end]))
    ylim(ax(2),[min(dataP),max(dataP)])

    ylabel(ax(2),'V','fontsize',15)
    xlabel(ax(2),'time (second) ','fontsize',15)
    set(ax(2),'fontsize',15)
    legend(ax(2),[l2,l2rbr],{'FPO7','RBR'},'location','northwest')
    
    
    
    l3=loglog(ax(3),f,Pa1,'b','linewidth',2);
    hold(ax(3),'on')
    l3rbr=loglog(ax(3),frbr,Pa1rbr,'r','linewidth',2);
    loglog(f,[err_low;err_high]*Pa1,'b--')
    n24=loglog(ax(3),f,Pa1*0+noise24,'g--');
    n16=loglog(ax(3),f,Pa1*0+noise16,'m--');
    J0=loglog(ax(3),f,Pa1*0+vJ0^2,'c--');
    xlim(ax(3),f([1 end]))
    ylim(ax(3),[.5*noise24, 5])
    legend([l3,l3rbr,n24,n16,J0],{sprintf('mean(FPO7 Spectrum)- N=%i',nbscan),'mean RBR spectrum','24 bit','16 bit','Johnson 0C'})
    ylabel(ax(3),'V^2 (T^2_{RBR}) / Hz','fontsize',15)
    xlabel('Hz','fontsize',15)
    set(ax(3),'fontsize',15)
    ax(3).YTick=[1e-12 1e-9 1e-6 1e-3 1];
    
    fig=gcf;
    fig.PaperPosition = [0 0 10 15];

    print('-dpng2',sprintf('../FIGURE/FPO7/FPO7_P%i.png',nb_profile))
    close all
    
end

close all
ax(3)=subplot('Position',[.1 .1 .8 .8]);
l3=loglog(ax(3),f,Ptot./length(Profiles),'b','linewidth',2);
hold(ax(3),'on')
n24=loglog(ax(3),f,Pa1*0+noise24,'g--');
n16=loglog(ax(3),f,Pa1*0+noise16,'r--');
J0=loglog(ax(3),f,Pa1*0+vJ0^2,'c--');
xlim(ax(3),f([1 end]))
legend([l3,n24,n16,J0],{sprintf('mean(Spectrum)- N=%s','all Profiles'),'24 bit','16 bit','Johnson 0C'})
ylabel(ax(3),'V^2 / Hz','fontsize',15)
xlabel('Hz','fontsize',15)
set(ax(3),'fontsize',15)
ax(3).YTick=[1e-9 1e-6 1e-3 1 1e3 1e6];

print('-dpng2','../FIGURE/FPO7/FPO7_tot.png')





