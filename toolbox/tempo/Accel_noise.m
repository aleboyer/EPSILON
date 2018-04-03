function Accel_noise(Profile,tscan,FS)
figure
if isfield(Profile,'P')==0
    disp('Add CTD data to Profile')
else
    
    %% plot frequency spectrum of the Accel and compare with bit and Spec noise floor
    
    %% Parameters fixed by data structure
    df        = 1/tscan;                                                       % number of samples per scan (1s) in channels
    %% Parameters for processing of microstructure data (epsilon and chi)
    f    = (df:df:FS/2)'; % frequency vector for spectra
    %bit noise
    Fn    = .5*FS;  % Nyquist frequency
    FR    = 3.3;      % Full range
    Nb_by = 24;     % number of bytes
    %[(full range = 3V)/2^24 ]^2 / Fn
    noise24= (FR/2^Nb_by)^2 /Fn;
    Nb_by = 16;     % number of bytes
    %[(full range = 3V)/2^24 ]^2 / Fn
    noise16= (FR/2^Nb_by)^2 /Fn;
    
    
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
    P3=zeros(nbscan,Lscan);
    gain_adc=1;
    sinc4_freq=f/2/f(end);
    hfilt = gain_adc.*(sinc(sinc4_freq)).^4;
    
    for i=1:nbscan
        %        data = Profile.Sensor1(indscan{i}(1:Lscan));
        %        datarbr = RBRProfile.T(indscanrbr{i}(1:Lscanrbr));
        data1 = Profile.Sensor6(indscan{i}(1:Lscan));
        data2 = Profile.Sensor7(indscan{i}(1:Lscan));
        data3 = Profile.Sensor8(indscan{i}(1:Lscan));
        A1=window.*data1(:);
        A2=window.*data2(:);
        A3=window.*data3(:);
        P1(i,:)=abs(fft(A1)).^2;
        P2(i,:)=abs(fft(A2)).^2;
        P3(i,:)=abs(fft(A3)).^2;
    end
    nu=2*floor(Lscan/nbscan);
    err_low = nu/chi2inv(.05/2,nu);
    err_high = nu/chi2inv(1-.05/2,nu);
    
    Pa1=mean(P1,1)./Lscan^2/f(1);
    Pa1=2*Pa1(1:floor(Lscan/2))./hfilt.';
    Pa2=mean(P2,1)./Lscan^2/f(1);
    Pa2=2*Pa2(1:floor(Lscan/2))./hfilt.';
    Pa3=mean(P3,1)./Lscan^2/f(1);
    Pa3=2*Pa3(1:floor(Lscan/2))./hfilt.';
    
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
    title(ax(1),['Acceleration noise ' datestr(nanmean(Profile.time))])
    
    
    hold(ax(2),'on')
    l2=plot(ax(2),Profile.time-Profile.time(1),Profile.Sensor6,'b','linewidth',2);
    l3=plot(ax(2),Profile.time-Profile.time(1),Profile.Sensor7,'r','linewidth',2);
    l4=plot(ax(2),Profile.time-Profile.time(1),Profile.Sensor8,'k','linewidth',2);
    hold(ax(2),'off')
    ylim(ax(2),[min([Profile.Sensor6(:);...
                     Profile.Sensor7(:);...
                     Profile.Sensor8(:)]),...
                max([Profile.Sensor6(:);...
                     Profile.Sensor7(:);...
                     Profile.Sensor8(:)])])

                     
                     
    ylabel(ax(2),'g','fontsize',15)
    xlabel(ax(2),'time (second) ','fontsize',15)
    set(ax(2),'fontsize',15)
    legend(ax(2),[l2,l3,l4],{'Ax','Ay','Az'},'location','northeast')
    
    l31=loglog(ax(3),f,Pa1,'b','linewidth',2);
    hold(ax(3),'on')
    l32=loglog(ax(3),f,Pa2,'r','linewidth',2);
    l33=loglog(ax(3),f,Pa3,'k','linewidth',2);
    loglog(f,[err_low;err_high]*Pa1,'b--','linewidth',2)
    Partnoise=loglog(ax(3),f,45e-6^2+0*f,'c--','linewidth',2);
    
    xlim(ax(3),f([1 end]))
    ylim(ax(3),[.5*noise24, 5])
    legend([l31,l32,l33,Partnoise],{sprintf('Ax - N=%i',nbscan),...
                                  sprintf('Ay - N=%i',nbscan),...
                                  sprintf('Az - N=%i',nbscan),...
                                  'Part Noise'},'location','west')
    
    ylabel(ax(3),'g^2 / Hz','fontsize',15)
    xlabel('Hz','fontsize',15)
    set(ax(3),'fontsize',15)
    ax(3).YTick=[1e-12 1e-9 1e-6 1e-3 1];
end
    
    
    
    
