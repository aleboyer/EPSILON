function FPO7_noise(Profile,tscan,FS)
figure

if isfield(Profile,'T')==0
    disp('Add CTD data to Profile')
else
    
    %% plot frequency spectrum of the FPO7 and compare with bit and Johnson Noise
    
    %% Parameters fixed by data structure
    df        = 1/tscan;                                                       % number of samples per scan (1s) in channels
    %% Parameters for processing of microstructure data (epsilon and chi)
    f    = (df:df:FS/2)'; % frequency vector for spectra
    %bit noise
    Fn    = .5*FS;  % Nyquist frequency
    FR    = 2.5;      % Full range
    Nb_by = 24;     % number of bytes
    %[(full range = 3V)/2^24 ]^2 / Fn
    noise24= (FR/2^Nb_by)^2 /Fn;
    Nb_by = 16;     % number of bytes
    %[(full range = 3V)/2^24 ]^2 / Fn
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
    
    
    CALFPO7_2=polyfit(Profile.Sensor2,Profile.T,3);
    CALFPO7_1=polyfit(Profile.Sensor1,Profile.T,3);
    FPO7_T1 = polyval(CALFPO7_1,Profile.Sensor1);
    FPO7_T2 = polyval(CALFPO7_2,Profile.Sensor2);
    
    %% Length of the Profile
    T       = length(Profile.time);
    %% define number of scan in the profile
    Lscan     = tscan*2*f(end);
    nbscan    = floor(T/Lscan);
    
    %% we compute spectra on scan with 50% overlap
    nbscan    = 2*nbscan-1;
    
    %% define the index in the profile for each scan
    indscan = arrayfun(@(x) (1+floor(Lscan/2)*(x-1):1+floor(Lscan/2)*(x-1)+Lscan-1),1:nbscan,'un',0);
    
    %% define fmax the cut off frequency after which the vibration of the
    window = hanning(Lscan);
    P1=zeros(nbscan,Lscan);
    P2=zeros(nbscan,Lscan);
    gain_adc=1;
    sinc4_freq=f/2/f(end);
    hfilt = gain_adc.*(sinc(sinc4_freq)).^4;
    
    for i=1:nbscan
        %        data = Profile.Sensor1(indscan{i}(1:Lscan));
        %        datarbr = RBRProfile.T(indscanrbr{i}(1:Lscanrbr));
        data1 = Profile.Sensor1(indscan{i}(1:Lscan));
        data2 = Profile.Sensor2(indscan{i}(1:Lscan));
        A1=window.*data1(:);
        A2=window.*data2(:);
        P1(i,:)=abs(fft(A1)).^2;
        P2(i,:)=abs(fft(A2)).^2;
    end
    nu=2*floor(Lscan/nbscan);
    err_low = nu/chi2inv(.05/2,nu);
    err_high = nu/chi2inv(1-.05/2,nu);
    
    Pa1=mean(P1,1)./Lscan^2/f(1);
    Pa1=2*Pa1(1:floor(Lscan/2))./hfilt.';
    Pa2=mean(P2,1)./Lscan^2/f(1);
    Pa2=2*Pa2(1:floor(Lscan/2))./hfilt.';
    
    %% plot coherence
    ax(1)=subplot(311);
    ax(2)=subplot(312);
    ax(3)=subplot(313);
    
    l1=plot(ax(1),Profile.time,Profile.P,'linewidth',2);
    axis(ax(1),'ij')
    xlim(ax(1),Profile.time([1 end]))
    ylim(ax(1),[min(Profile.P), max(Profile.P)])
    ylabel(ax(1),'Presure (db)','fontsize',15)
    set(ax(1),'fontsize',15)
    title(ax(1),['FPO7 Noise ' datestr(nanmean(Profile.time))])
    
    
    l2ctd=plot(ax(2),Profile.time-Profile.time(1),Profile.T,'g','linewidth',2);
    hold(ax(2),'on')
    l2=plot(ax(2),Profile.time-Profile.time(1),FPO7_T1,'b','linewidth',2);
    l3=plot(ax(2),Profile.time-Profile.time(1),FPO7_T2,'r','linewidth',2);
    hold(ax(2),'off')
    ylim(ax(2),[min(FPO7_T2),max(FPO7_T2)])
    
    ylabel(ax(2),'Celcius','fontsize',15)
    xlabel(ax(2),'time (second) ','fontsize',15)
    set(ax(2),'fontsize',15)
    legend(ax(2),[l2,l3,l2ctd],{'FPO7 1','FPO7 2','CTD'},'location','northeast')
    
    l31=loglog(ax(3),f,Pa1,'b','linewidth',2);
    hold(ax(3),'on')
    l32=loglog(ax(3),f,Pa2,'r','linewidth',2);
    loglog(f,[err_low;err_high]*Pa1,'b--','linewidth',2)
    n24=loglog(ax(3),f,Pa1*0+noise24,'g--','linewidth',2);
    n16=loglog(ax(3),f,Pa1*0+noise16,'m--','linewidth',2);
    J0=loglog(ax(3),f,Pa1*0+vJ0^2,'c--','linewidth',2);
    xlim(ax(3),f([1 end]))
    ylim(ax(3),[.5*noise24, 5])
    legend([l31,l32,n24,n16,J0],{sprintf('FPO7 1 - N=%i',nbscan),sprintf('FPO7 2 - N=%i',nbscan),'24 bit','16 bit','Johnson 0 \circ C'},'location','west')
    ylabel(ax(3),'Volt^2 / Hz','fontsize',15)
    xlabel('Hz','fontsize',15)
    set(ax(3),'fontsize',15)
    ax(3).YTick=[1e-12 1e-9 1e-6 1e-3 1];
end
    
    
    
    
