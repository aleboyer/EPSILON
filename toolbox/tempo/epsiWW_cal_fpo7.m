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

Ptot=0;
Ptotrbr=0;
j=200

nb_profile=j;
Profile=Profiles{j};
RBRProfile=RBRprofiles{j};

dataP=Profile.Sensor1;
datarbrP=RBRProfile.T;

timeepsi=RBRProfile.time(1):1/325/86400:RBRProfile.time(1)+(length(dataP)/325-1/325)/86400;

%% plot data
timecaleps=1.142*(timeepsi-timeepsi(1))*86400;
timecalrbr=86400*(RBRProfile.time-RBRProfile.time(1));

dataP=interp1(timecaleps,dataP,timecalrbr);
datarbrP=datarbrP(~isnan(dataP));
timecalrbr=timecalrbr(~isnan(dataP));
dataP=dataP(~isnan(dataP));

CALFPO7=polyfit(dataP,datarbrP,3);

ax(1)=subplot(411);
ax(2)=subplot(412);
ax(3)=subplot(413);
ax(4)=subplot(414);

l2rbr=plot(ax(1),timecalrbr,datarbrP,'r','linewidth',2);
l2=plot(ax(2),timecalrbr,dataP,'b','linewidth',2);

xlim(ax(1),time1([1 end]))
ylim(ax(1),[min(datarbrP),max(datarbrP)])
xlim(ax(2),time1([1 end]))
ylim(ax(2),[min(dataP),max(dataP)])

ylabel(ax(1),'Temperature (^{\circ}C)','fontsize',15)
ylabel(ax(2),'Volt','fontsize',15)
xlabel(ax(2),'time (second) ','fontsize',15)
axis(ax(2),'ij')
legend(ax(1),l2rbr,'RBR','location','northwest')
legend(ax(2),l2,'FPO7','location','northwest')

plot(ax(3),dataP,datarbrP,'linewidth',3)
hold(ax(3),'on')
plot(ax(3),dataP,polyval(CALFPO7,dataP),'r','linewidth',2)
text(1.6,datarbrP(end-10),'CALIBRATION POLYFIT n=3','fontsize',20,'Parent',ax(3))
hold(ax(3),'off')
ylabel(ax(3),'Temperature (^{\circ}C)','fontsize',15)
xlabel(ax(3),'Volt','fontsize',15)

l2rbr=plot(ax(4),timecalrbr,datarbrP,'r','linewidth',2);
hold(ax(4),'on')
l2=plot(ax(4),timecaleps,polyval(CALFPO7,Profile.Sensor1),'b','linewidth',2);
hold(ax(4),'off')
xlim(ax(4),time1([1 end]))
ylim(ax(4),[min(datarbrP),max(datarbrP)])
ylabel(ax(4),'Temperature (^{\circ}C)','fontsize',15)
xlabel(ax(4),'time (second) ','fontsize',15)
legend(ax(4),[l2rbr,l2],{'RBR','FPO7'},'location','northwest')

set(ax(1),'fontsize',15)
set(ax(2),'fontsize',15)
set(ax(3),'fontsize',15)

fig=gcf;
fig.PaperPosition = [0 0 10 10];
print('-dpng','../FIGURE/calibration_FPO7.png')


Profile = compute_fallrate(Profile);


%% Length of the Profile

TFPO7= polyval(CALFPO7,Profile.Sensor1);
TRBR = datarbrP(timecalrbr<=timecaleps(end));
timerbr=timecalrbr(timecalrbr<=timecaleps(end));


%% Parameters fixed by data structure
tscan     =  5;
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





T       = length(Profile.time);
Trbr    = length(timerbr);
df      = f(1);
%% define number of scan in the profile
Lscan     = tscan*2*f(end);
nbscan    = floor(T/Lscan);

%% we compute spectra on scan with 50% overlap
nbscan    = 2*nbscan-1;

%% define the index in the profile for each scan
indscan = arrayfun(@(x) (1+floor(Lscan/2)*(x-1):1+floor(Lscan/2)*(x-1)+Lscan-1),1:nbscan,'un',0);

for i=1:length(indscan)
    indscanrbr{i} = find(timerbr>=timecaleps(indscan{i}(1)) & ...
                         timerbr<=timecaleps(indscan{i}(end)) );
    if length(indscanrbr{i})==34 
        % this is bad
        %I force the length of indscan rbr to 35 
        % because I know that it will oscillate between 34 and 35
        % I simply want a a structure of index array that will match indscan 
        % with the same size array.
        indscanrbr{i}=[indscanrbr{i};indscanrbr{i}(end)+1];
    end
end




%% define the index in the profile for each scan
total_w       = cellfun(@(x) nanmean(Profile.w(x)),indscan); 
ind_upcast = find((-total_w)>20);
Epsilon.indscan = indscan(ind_upcast);
Epsilon.w       = cellfun(@(x) nanmean(Profile.w(x)),indscan(ind_upcast)); 
Epsilon.pr      = cellfun(@(x) nanmean(Profile.P(x)),indscan(ind_upcast)); % needed to compute Kvis 



%% define fmax the cut off frequency after which the vibration of the

P=zeros(nbscan,length(indscan{1})/2);
Prbr=zeros(nbscan,floor(length(indscanrbr{1})/2));

for i=1:nbscan
    data = TFPO7(indscan{i});
    data = data-mean(data);
    datarbr = TRBR(indscanrbr{i});
    datarbr = datarbr-mean(datarbr);
    [P(i,:),f1] = pwelch(data,length(indscan{i}),0,f,2*f(end));

    df     = 1./timerbr(indscanrbr{i}(end)-indscanrbr{i}(1));               % sample rate channels
    frbr = (df:df:FSrbr/2+df)'; % frequency vector for spectra
    [Prbr(i,:),f1rbr] = pwelch(datarbr,length(indscanrbr{i}),0,frbr,2*frbr(end));

end
i=ind_upcast(30)
data = TFPO7(indscan{i});
data = data-mean(data);
datarbr = TRBR(indscanrbr{i});
datarbr = datarbr-mean(datarbr);

figure
subplot(311)
plot(timecaleps(indscan{i}),data,'b','linewidth',2)
hold on
plot(timerbr(indscanrbr{i}),datarbr,'r','linewidth',2)
legend('FPO7','RBR','location','Northwest')
xlabel('second','fontsize',20)
ylabel('^{\circ}C','fontsize',20)
set(gca,'fontsize',15)

gain_adc=1;
sinc4_freq=f/2/f(end);
hfilt = gain_adc.*(sinc(sinc4_freq)).^4;
k=f1/speed;

speed=-Epsilon.w(i)/100; % convert to m/s
tau=0.005 * speed^(-0.32); % thermistor time constant
magsq=1 ./ (1+((2*pi*tau).*f).^2).^2; % magnitude-squared
phase=-2*atan( 2*pi*f*tau);


h_th_total=hfilt .* magsq(:);

subplot(312)
loglog(f1,magsq,'b','linewidth',2)
hold on
loglog(f1,hfilt,'r','linewidth',2)
loglog(f1,h_th_total,'k','linewidth',2)
legend('H_{FPO7}','sinc4','total','location','Northwest')
xlabel('Hz','fontsize',20)
ylabel('H','fontsize',20)
set(gca,'fontsize',20)


subplot(313)
loglog(f1,P(i,:),'b','linewidth',2)
hold on
plot(f1rbr,Prbr(i,:),'r','linewidth',2)
ylim([min(P(i,:)),max(P(i,:))])
set(gca,'Ytick',[1e-10 1e-8 1e-6 1e-4])

lTF=plot(f1,h_th_total./max(h_th_total)*max(P(i,:)),'k','linewidth',2);
legend(lTF,'Transfer function')
xlabel('Hz','fontsize',20)
ylabel('^{\circ}C^2 /Hz','fontsize',20)
set(gca,'fontsize',15)

fig=gcf;
fig.PaperPosition = [0 0 8 8];
print('-dpng','../FIGURE/Temp_spectra_cal.png')



Ptempk = Emp_Corr_fac * (P(i,:)*speed)./h_th_total';
Ptgradk=(2*pi.*k).^2 .* Ptempk(:) ;

% TODO implement 
Epsilon.indscan = total_indscan(ind_upcast);
Epsilon.w       = cellfun(@(x) nanmean(Profile.w(x)),indscan(ind_upcast)); 
Epsilon.t       = cellfun(@(x) nanmean(Profile.T(x)),indscan(ind_upcast)); % needed to compute Kvis 
Epsilon.s       = cellfun(@(x) nanmean(Profile.S(x)),indscan(ind_upcast)); % needed to compute Kvis 
Epsilon.pr      = cellfun(@(x) nanmean(Profile.P(x)),indscan(ind_upcast)); % needed to compute Kvis 

i=30

%% get cut off frequency
str=['fc_index=' algorithm '_' electronicsid '_mmp' ...
    '(fth,Pth,displ_chi_spec,pr_chi(j),j,speed);'];
displ_chi_spec='no';
fc_index=FPO7_cutoff(f,P(i,:),displ_chi_spec,Epsilon.pr(i),i,speed);
% scale with wavenumber
kcth=f(fc_index)/speed;
dk=df/speed;


ktemp=kt(Epsilon.s(i),Epsilon.t(i),Epsilon.pr(i));

chi=6*ktemp*dk.*sum(Ptgradk(1:fc_index));
eps_chi(j,i) = abs(n2(idn2)*chi(j,i)/(2*0.2*dthetadz(idn2)^2));




figure
loglog(k,Ptempk,'b','linewidth',2)
hold on
loglog(k,Ptgradk,'r','linewidth',2)
loglog([kcth kcth],[min(Ptempk),max(Ptgradk)],'k--')
text(kcth,max(Ptgradk),'Cut off wavenumber','fontsize',20)
text(3e-1,6e-15,['\chi' sprintf('=%1.3e K^2 s^{-1}',chi)],'fontsize',20)
xlabel('cpm','fontsize',20)
%ylabel('\color{blue}\,^{\circ}C^2  \,^{\circ}C^2m^{-2}_{red} /cpm','interpreter','tex','fontsize',20)
ylabel('\color{blue} \circ C^2 \color{red}\circ C^2m^{-2} \color{black}/ cpm','interpreter','tex','fontsize',20)
%title(['\fontsize{16}black {\color{magenta}magenta ','\color[rgb]{0 .5 .5}teal \color{red}red} black again'],'interpreter','tex')
set(gca,'fontsize',20)
legend('\phi_T','\phi_{TG}','location','northwest')
print('-dpng','../FIGURE/Tempgrad_WNB_spectra_cal.png')





nu=2*floor(Lscan/nbscan);
err_low = nu/chi2inv(.05/2,nu);
err_high = nu/chi2inv(1-.05/2,nu);

nurbr=2*floor(Lscanrbr/nbscanrbr);
err_lowrbr = nurbr/chi2inv(.05/2,nurbr);
err_highrbr = nurbr/chi2inv(1-.05/2,nurbr);





