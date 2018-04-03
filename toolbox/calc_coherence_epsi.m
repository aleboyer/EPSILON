function [Pco12,Pco13,Pco14,fig]=calc_coherence_epsi(Profile,tscan,f,nb_profile)

%% Length of the Profile
T       = length(Profile.time);
df      = f(1);
%% define number of scan in the profile
Lscan   = tscan*2*f(end);
nbscan  = floor(T/Lscan);

%% we compute spectra on scan with 50% overlap
nbscan=2*nbscan-1;

%% define the fall rate of the Profile. 
%  Add Profile.w with w the vertical vel. 
%  We are using the pressure from other sensors (CTD);
Profile = compute_fallrate(Profile);

%% define the index in the profile for each scan
Epsilon.indscan = arrayfun(@(x) (1+floor(Lscan/2)*(x-1):1+floor(Lscan/2)*(x-1)+Lscan-1),1:nbscan,'un',0);
Epsilon.w       = cellfun(@(x) nanmean(Profile.w(x)),Epsilon.indscan,'un',0); 
Epsilon.pr      = cellfun(@(x) nanmean(Profile.P(x)),Epsilon.indscan,'un',0);

%% define fmax the cut off frequency after which the vibration of the 
window = hanning(Lscan);
wc2    = 1/mean(window.^2);% window correction factor
P1=zeros(nbscan,Lscan);
P2=zeros(nbscan,Lscan);
P3=zeros(nbscan,Lscan);
P4=zeros(nbscan,Lscan);
testP1=zeros(nbscan,Lscan);
P11=zeros(nbscan,Lscan);
P22=zeros(nbscan,Lscan);
P33=zeros(nbscan,Lscan);
P44=zeros(nbscan,Lscan);
window=hanning(Lscan);
%for j=1:nbscan
for j=50:nbscan-50
    data = Profile.Sensor3(Epsilon.indscan{j}(1:Lscan));
    Az   = Profile.Sensor6(Epsilon.indscan{j}(1:Lscan));
    Ax   = Profile.Sensor7(Epsilon.indscan{j}(1:Lscan));
    Ay   = Profile.Sensor8(Epsilon.indscan{j}(1:Lscan));
%    A1=window.*data(:);A2=window.*Az(:);A3=window.*Ax(:);A4=window.*Ay(:);
    A1=data(:);A2=Az(:);A3=Ax(:);A4=Ay(:);
    P1(j,:)=fft(A1);P2(j,:)=fft(A2);P3(j,:)=fft(A3);P4(j,:)=fft(A4);
    %P11(j,:)=abs(fft(A1)).^2/Lscan.^2/df*wc2;
    P11(j,:)=abs(fft(A1)).^2/Lscan;
    
    P22(j,:)=abs(fft(A2)).^2/Lscan;
    P33(j,:)=abs(fft(A3)).^2/Lscan;
    P44(j,:)=abs(fft(A4)).^2/Lscan;
    
end
P12=conj(P1).*P2./(Lscan)^2/df*wc2;
P13=conj(P1).*P3./(Lscan)^2/df*wc2;
P14=conj(P1).*P4./(Lscan)^2/df*wc2;
CoP12=P12.^2./(P11.*P22);
CoP13=P13.^2./(P11.*P33);
CoP14=P14.^2./(P11.*P44);

CoP12(CoP12==0)=nan;
CoP13(CoP13==0)=nan;
CoP14(CoP14==0)=nan;

CoP12=abs(nanmean(CoP12,1));
CoP13=abs(nanmean(CoP13,1));
CoP14=abs(nanmean(CoP14,1));

P11(P11==0)=nan;
P22(P22==0)=nan;
P33(P33==0)=nan;
P44(P44==0)=nan;

Pa1=nanmean(P11,1);
Pa2=nanmean(P22,1);
Pa3=nanmean(P33,1);
Pa4=nanmean(P44,1);
Pa1=2*Pa1(1:Lscan/2);
Pa2=2*Pa2(1:Lscan/2);
Pa3=2*Pa3(1:Lscan/2);
Pa4=2*Pa4(1:Lscan/2);

Pco12=CoP12(1:Lscan/2);
Pco13=CoP13(1:Lscan/2);
Pco14=CoP14(1:Lscan/2);

%% plot coherence 
ax(1)=subplot(411);
ax(2)=subplot(412);
ax(3)=subplot(413);
ax(4)=subplot(414);

%l1=plot(ax(1),Profile.rbrtime,Profile.P,'linewidth',2);
l1=plot(ax(1),Profile.time,Profile.P,'linewidth',2);
datetick(ax(1))
axis(ax(1),'ij')
%xlim(ax(1),Profile.rbrtime([1 end]))
xlim(ax(1),Profile.time([1 end]))
ylim(ax(1),[min(Profile.P), max(Profile.P)])
ylabel(ax(1),'Presurre (db)','fontsize',15)
set(ax(1),'fontsize',15)
%title(ax(1),['Profile ' num2str(nb_profile) '-' datestr(nanmean(Profile.rbrtime))])
title(ax(1),['Profile ' num2str(nb_profile) ])

channels=detrend([Profile.Sensor6;Profile.Sensor7;Profile.Sensor8;Profile.Sensor3;Profile.Sensor4;].','constant').';
%l2=plot(ax(2),Profile.rbrtime,channels,'linewidth',2);
l2=plot(ax(2),Profile.time,channels,'linewidth',2);
datetick(ax(2))
%xlim(ax(2),Profile.rbrtime([1 end]))
xlim(ax(2),Profile.time([1 end]))
ylim(ax(2),[min(min(channels(:,Lscan:end-Lscan))),max(max(channels(:,Lscan:end-Lscan)))])
ylabel(ax(2),'V and g','fontsize',15)
xlabel(ax(2),'time ','fontsize',15)
set(ax(2),'fontsize',15)
legend(ax(2),l2,{'Az','Ax','Ay','s1','s2'},'location','northwest')


Fn    = .5*f(end);  % Nyquist frequency
FR    = 2.5;        % Full range
Nb_by = 24;         % number of bytes
%[(full range = 3V)/2^24 ]^2 / Fn
noise24= (FR/2^Nb_by)^2 /Fn;
Nb_by = 16;     % number of bytes
%[(full range = 3V)/2^24 ]^2 / Fn
noise16= (FR/2^Nb_by)^2 /Fn;    
h_freq=get_filters_MADRE(f);
speed=.65;
Sv=48.66;G=9.81
htotal = h_freq.shear .* haf_oakey(f,speed);        % should add epsi filter

l3=loglog(ax(3),f,[Pa1;Pa2;Pa3;Pa4],'linewidth',2);
hold(ax(3),'on')
%loglog(ax(3),f,f*0+noise16,'linewidth',2);
hold(ax(3),'off')

xlim(ax(3),f([1 end]))
ylim(ax(3),[min(min([Pa1;Pa2;Pa3;Pa4])), max(max([Pa1;Pa2;Pa3;Pa4]))])
ylabel(ax(3),'V^2 and g^2 / Hz','fontsize',15)
set(ax(3),'Xscale','log','Yscale','log')
legend('Shear','Ax','','')
set(ax(3),'fontsize',15)

l4=semilogx(ax(4),f,[Pco12;Pco13;Pco14],'linewidth',2);
xlim(ax(4),f([1 end]))
ylim(ax(4),[min(min([Pco12;Pco13;Pco14])), max(max([Pco12;Pco13;Pco14]))])
legend(ax(4),l4,{'Shear-Az','Shear-Ax','Shear-Ay'})
xlabel(ax(4),'Hz','fontsize',15)
ylabel(ax(4),'Coherence','fontsize',15)
set(ax(4),'fontsize',15)
print('-dpng2',sprintf('../FIGURE/COHERENCE/Coherence_P%i.png',nb_profile))

