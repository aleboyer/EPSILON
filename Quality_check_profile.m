question=input('redownload? (1:oui)');

addpath toolbox/
addpath toolbox/CTD
addpath toolbox/FILTER
addpath toolbox/process
addpath toolbox/seawater2/
addpath toolbox/PLOTS/

if question == 1
    clear;
%    root_dir='/Users/aleboyer/ARNAUD/SCRIPPS/SP1810/WW/d1/L1/';
    %root_dir='/Users/aleboyer/ARNAUD/SCRIPPS/EPSIWW/WW/EPSI/d1/L1/';
    %root_dir='/Users/aleboyer/ARNAUD/SCRIPPS/SP1810/EPSIfish/d9/L1/';
    %root_dir='/Users/aleboyer/ARNAUD/SCRIPPS/pregranite/epsifish/d3/L1/';
    %root_dir='/Users/aleboyer/ARNAUD/SCRIPPS/NISKINE/epsifish2/d3/L1/';
    root_dir='/Users/aleboyer/ARNAUD/SCRIPPS/NISKINE/WW/d2/L1/';
%    load('/Users/aleboyer/ARNAUD/SCRIPPS/pregranite/epsifish/d3/L1/Profiles_epsifish_ctd_d3.mat') 
%    load('/Users/aleboyer/ARNAUD/SCRIPPS/NISKINE/epsifish2/d3/L1/Profiles_epsifish2_ctd_d3.mat')
    load([ root_dir 'Profiles_WW_ctd_d2.mat'])
%     load([root_dir 'Profiles_epsifish2_ctd_d3.mat'])
%     load([root_dir 'Profiles_EPSIfish_ctd_d9.mat'])
%     load([root_dir 'Profiles_epsifish_ctd_d3.mat'])
%      EPSIWW february 2017
%     load('/Users/aleboyer/ARNAUD/SCRIPPS/EPSIWW/WW/EPSI/d1/epsi/Profiles_EPSI_rbr_d1.mat')
%     load('/Users/aleboyer/ARNAUD/SCRIPPS/EPSIWW/WW/EPSI/d1/rbr/Profiles_EPSI_rbr_d1.mat')
    load([root_dir 'Turbulence_Profiles.mat'])
    
else
    clear MS1
    close all
    clear data
end
tscan=15;
%EpsiProfiles=EpsiProfile.datadown;
%CTDProfiles=CTDProfile.datadown;

% Trick for NISKINE d2 WW
%CTDProfiles=CTDProfiles(10:33);

for id_profile=1:length(EpsiProfiles)
%for id_profile=4
    try
        clear data
        %id_profile=180;
        %Sv=[53.00,57.57];
        Sv=[42.14,56.71];
        %Sv=[53.00 48.16];
        
        %CTDProfiles=CTDProfiles.datadown;
        flag_vehicle=-1; % -1 if wirewalker 1 if epsifish
        
        %Profile=EPSI_Profiles{id_profile};
        %CTDProfiles=RBRProfile.dataup;
        %EpsiProfiles=EpsiProfile.dataup;
        %EpsiProfiles{id_profile}.time=EpsiProfiles{id_profile}.rbrtime; % epsiWW february
        
        % Trick for SP1810 d1 WW
        %EpsiProfiles{id_profile}.time=EpsiProfiles{id_profile}.EPSItime;
        % Trick for NISKINE d2 WW
        EpsiProfiles{id_profile}.time=EpsiProfiles{id_profile}.EPSItime;

        
        Profile=EpsiProfiles{id_profile};
        % for NISKINE espifish2
%         [time,IA]=unique(CTDProfiles{id_profile}.time);
%         Profile.T=interp1(time, ...
%             CTDProfiles{id_profile}.T(IA), ...
%             EpsiProfiles{id_profile}.time);
%         Profile.S=interp1(time, ...
%             CTDProfiles{id_profile}.S(IA), ...
%             EpsiProfiles{id_profile}.time);
%         Profile.P=interp1(time, ...
%             CTDProfiles{id_profile}.P(IA), ...
%             EpsiProfiles{id_profile}.time);
        
        % for WW
        Profile.T=interp1(CTDProfiles{id_profile}.time, ...
            CTDProfiles{id_profile}.T, ...
            EpsiProfiles{id_profile}.time);
        Profile.S=interp1(CTDProfiles{id_profile}.time, ...
            CTDProfiles{id_profile}.S, ...
            EpsiProfiles{id_profile}.time);
        Profile.P=interp1(CTDProfiles{id_profile}.time, ...
            CTDProfiles{id_profile}.P, ...
            EpsiProfiles{id_profile}.time);
        f=1/15:1/15:320/2;
        
        %% Gravity  ... of the situation :)
        G       = 9.81;
        
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
        Profile = compute_fallrate_downcast(Profile);
        Profile.w=flag_vehicle*Profile.w;
        if nanmedian(Profile.w)>10
            Profile.w=Profile.w/100;
        end
        
        % on line calibration of the FPO7s
        ind_nonan1=find(~isnan(Profile.Sensor1));
        ind_nonan2=find(~isnan(Profile.Sensor2));
        CALFPO7_1=polyfit(Profile.Sensor1(ind_nonan1),Profile.T(ind_nonan1),3);
        CALFPO7_2=polyfit(Profile.Sensor2(ind_nonan2),Profile.T(ind_nonan2),3);
        %CALFPO7_2=polyfit(Profile.Sensor2,Profile.T,1);
        
        %TODO probably check nan at previous step
        Profile.Sensor1=fillmissing(Profile.Sensor1,'linear');
        Profile.Sensor2=fillmissing(Profile.Sensor2,'linear');
        Profile.Sensor3=fillmissing(Profile.Sensor3,'linear');
        Profile.Sensor4=fillmissing(Profile.Sensor4,'linear');
        Profile.Sensor5=fillmissing(Profile.Sensor5,'linear');
        Profile.Sensor6=fillmissing(Profile.Sensor6,'linear');
        Profile.Sensor7=fillmissing(Profile.Sensor7,'linear');
        Profile.Sensor8=fillmissing(Profile.Sensor8,'linear');
        
        %EPSIWW
        Profile.Sensor5=double(Profile.Sensor5);
        
        Profile.Sensor1=filloutliers(Profile.Sensor1,'linear');
        Profile.Sensor2=filloutliers(Profile.Sensor2,'linear');
        Profile.Sensor3=filloutliers(Profile.Sensor3,'linear');
        Profile.Sensor4=filloutliers(Profile.Sensor4,'linear');
        Profile.Sensor5=filloutliers(Profile.Sensor5,'linear');
        Profile.Sensor6=filloutliers(Profile.Sensor6,'linear');
        Profile.Sensor7=filloutliers(Profile.Sensor7,'linear');
        Profile.Sensor8=filloutliers(Profile.Sensor8,'linear');
        
        
        %% define the index in the profile for each scan
        total_indscan = arrayfun(@(x) (1+floor(Lscan/2)*(x-1):1+floor(Lscan/2)*(x-1)+Lscan-1),1:nbscan,'un',0);
        total_w       = cellfun(@(x) nanmean(Profile.w(x)),total_indscan);
        
        ind_downcast = find((total_w)>.20);
        nbscan=length(ind_downcast);
        
        MS1.indscan    = total_indscan(ind_downcast);
        MS1.nbscan    = nbscan;
        MS1.fmax      = 45; % arbitrary cut off frequency usually extract from coherence spectra shear/accel
        MS1.nbchannel = nbscan;
        
        MS1.w       = cellfun(@(x) nanmean(Profile.w(x)),total_indscan(ind_downcast));
        MS1.t       = cellfun(@(x) nanmean(Profile.T(x)),total_indscan(ind_downcast)); % needed to compute Kvis
        MS1.s       = cellfun(@(x) nanmean(Profile.S(x)),total_indscan(ind_downcast)); % needed to compute Kvis
        MS1.pr      = cellfun(@(x) nanmean(Profile.P(x)),total_indscan(ind_downcast)); % needed to compute Kvis
        MS1.time    = cellfun(@(x) nanmean(Profile.time(x)),total_indscan(ind_downcast)); % needed to compute Kvis
        
        % loop to estimate eps
        % get and split the data
        data(1,:,:) = cell2mat(cellfun(@(x) polyval(CALFPO7_1,Profile.Sensor1(x)),MS1.indscan,'un',0).');
        data(2,:,:) = cell2mat(cellfun(@(x) polyval(CALFPO7_2,Profile.Sensor2(x)),MS1.indscan,'un',0).');
        data(3,:,:) = cell2mat(cellfun(@(x) Profile.Sensor3(x),MS1.indscan,'un',0).');
        data(4,:,:) = cell2mat(cellfun(@(x) Profile.Sensor4(x),MS1.indscan,'un',0).');
        % acceleration divided by speed to get the same units (s-1)
        data(5,:,:) = cell2mat(cellfun(@(x,y) Profile.Sensor6(x)./y,MS1.indscan,num2cell(MS1.w),'un',0).');
        data(6,:,:) = cell2mat(cellfun(@(x,y) Profile.Sensor7(x)./y,MS1.indscan,num2cell(MS1.w),'un',0).');
        data(7,:,:) = cell2mat(cellfun(@(x,y) Profile.Sensor8(x)./y,MS1.indscan,num2cell(MS1.w),'un',0).');
        
        
        % Profile Power and Co spectrum and Coherence. (Coherence still needs to be averaged over few scans afterwork)
        [f1,P1,P11,Co12]=get_profile_spectrum(data,f);
        %TODO comment on the Co12 sturcutre and think about reducing the size of
        %the Coherence spectra (doublon)
        
        indf1=find(f1>=0);
        indf1=indf1(1:end-1);
        f1=f1(indf1);
        Lf1=length(indf1);
        
        P11=2*P11(:,:,indf1);
        %% get MADRE filters
        h_freq=get_filters_MADRE('MADRE2.1',f1);
        
        %% correct transfert functions for accel spectra
        P11(5:7,:,:)=P11(5:7,:,:)./...
            shiftdim(repmat(ones(nbscan,1)*h_freq.electAccel,[1,1,3]),2).^2;
        
        %% correct transfert functions for shear spectra
        
        TF1 =@(x) (Sv.'.*x/(2*G)).^2 .* h_freq.shear .* haf_oakey(f1,x);     % should add epsi filter
        TFshear=cell2mat(cellfun(@(x) TF1(x),num2cell(MS1.w),'un',0).');
        TFshear=reshape(TFshear,[2,nbscan,Lf1]);
        P11(3:4,:,:) = P11(3:4,:,:) ./ TFshear;      % vel frequency spectra m^2/s^-2 Hz^-1
        
        %% correct transfert functions for temperature spectra
        Emp_Corr_fac=1;
        TFtemp=cell2mat(cellfun(@(x) h_freq.FPO7(x),num2cell(MS1.w),'un',0).');
        TFtemp=shiftdim(repmat(TFtemp,[1,1,2]),2);
        P11(1:2,:,:) = Emp_Corr_fac * P11(1:2,:,:)./TFtemp;
        
        %f=1/3:1/3:320/2;
        %Epsilon_corrected=calc_eps_epsi_accel_corrected(Profile,3,f,Sv);
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        ax(1)=subplot('Position',[.48 .62 .5 .3]);
        loglog(f(1:end-1),smoothdata(nanmedian(squeeze(P11(3,:,:)))),'linewidth',2)
        hold on
        loglog(f(1:end-1),smoothdata(nanmedian(squeeze(P11(4,:,:)))),'linewidth',2)
        loglog(f(1:end-1),smoothdata(nanmedian(squeeze(P11(5,:,:)))),'linewidth',2)
        loglog(f(1:end-1),smoothdata(nanmedian(squeeze(P11(6,:,:)))),'linewidth',2)
        loglog(f(1:end-1),smoothdata(nanmedian(squeeze(P11(7,:,:)))),'linewidth',2)
        
        minY=min(smoothdata(nanmedian(squeeze(P11(3,:,:)))));
        maxY=max(smoothdata(nanmedian(squeeze(P11(3,:,:)))));
        minY=min([minY min(smoothdata(nanmedian(squeeze(P11(4,:,:)))))]);
        maxY=max([maxY max(smoothdata(nanmedian(squeeze(P11(4,:,:)))))]);
        minY=min([minY min(smoothdata(nanmedian(squeeze(P11(5,:,:)))))]);
        maxY=max([maxY max(smoothdata(nanmedian(squeeze(P11(5,:,:)))))]);
        minY=min([minY min(smoothdata(nanmedian(squeeze(P11(6,:,:)))))]);
        maxY=max([maxY max(smoothdata(nanmedian(squeeze(P11(6,:,:)))))]);
        
        hold off
        legend('Shear1','Shear2','Accel X/speed','Accel Y/speed','Accel Z/speed','location','eastoutside')
%        legend('Shear1','Shear2','Accel X/speed','Accel Y/speed','location','eastoutside')
        ylabel('s^{-2}/Hz','fontsize',20)
        xlabel('Hz','fontsize',20)
        grid on
        xlim(f([1 end]))
        ylim([minY maxY])
        set(ax(1),'fontsize',15)
        
        ax(2)=subplot('Position',[.48 .3 .5 .2]);
        semilogx(f(1:end-1),smoothdata(abs(mean(squeeze(Co12(3,4,:,indf1))))),'linewidth',2)
        hold on
        semilogx(f(1:end-1),smoothdata(abs(mean(squeeze(Co12(3,5,:,indf1))))),'linewidth',2)
        semilogx(f(1:end-1),smoothdata(abs(mean(squeeze(Co12(3,6,:,indf1))))),'linewidth',2)
        hold off
        legend('Sh1-AccelX','Sh1-AccelY','Sh1-AccelZ','location','eastoutside')
%        legend('Sh1-AccelX','Sh1-AccelY','location','eastoutside')
        title('Coherence','fontsize',20)
        xlim(f([1 end]))
        set(gca,'XtickLabel',[])
        set(ax(2),'fontsize',15)
        grid on
        
        ax(3)=subplot('Position',[.48 .08 .5 .2]);
        semilogx(f(1:end-1),smoothdata(abs(mean(squeeze(Co12(4,4,:,indf1))))),'linewidth',2)
        hold on
        semilogx(f(1:end-1),smoothdata(abs(mean(squeeze(Co12(4,5,:,indf1))))),'linewidth',2)
        semilogx(f(1:end-1),smoothdata(abs(mean(squeeze(Co12(4,6,:,indf1))))),'linewidth',2)
        hold off
%        legend('Sh2-AccelX','Sh2-AccelY','location','eastoutside')
        legend('Sh2-AccelX','Sh2-AccelY','Sh2-AccelZ','location','eastoutside')
        xlim(f([1 end]))
        set(ax(3),'fontsize',15)
        xlabel('Hz','fontsize',20)
        grid on
        
        
        ax(4)=subplot('Position',[.1 .57 .15 .35]);
        semilogx(smoothdata(MS{id_profile}.epsilon(:,1)),MS{id_profile}.pr,'linewidth',2)
        hold(ax(4),'on')
        semilogx(smoothdata(MS{id_profile}.epsilon(:,2)),MS{id_profile}.pr,'linewidth',2)
        hold(ax(4),'off')
        legend('epsilon1','epsilon2','location','southeast')
        axis ij
        ylim(sort(MS{id_profile}.pr([1 end])))
        xlim([min(MS{id_profile}.epsilon(:)) ...
            max(MS{id_profile}.epsilon(:))])
        xlabel('\epsilon W.kg^{-1}','fontsize',15)
        ylabel('Depth (m)','fontsize',15)
        grid on
        
        ax(5)=subplot('Position',[.3 .57 .1 .35]);
        semilogx(MS{id_profile}.w,MS{id_profile}.pr,'linewidth',2)
        axis ij
        ylim(sort(MS{id_profile}.pr([1 end])))
        xlim([min(MS{id_profile}.w(:)) max(MS{id_profile}.w(:))])
        xlabel('speed (m s^{-1})','fontsize',15)
        ylabel('Depth (m)','fontsize',15)
        grid on
        
        Epsilon_class=calc_binned_epsi(MS(id_profile));
        
        ax(6)=subplot('Position',[.08 .3 .37 .2]);
        ax(7)=subplot('Position',[.08 .1 .37 .2]);
        
        plot_binned_epsilon_sanity_profile(Epsilon_class,' ',ax(6),ax(7))
        
        fig=gcf;
        fig.PaperPosition=[0 0 15 10];
        print(sprintf('%sSanity_Profile_%i',root_dir,id_profile),'-dpng2')
        pause(.1)
        cla(ax(1));
        cla(ax(2));
        cla(ax(3));
        cla(ax(4));
        cla(ax(5));
        cla(ax(6));
        cla(ax(7));
    catch
        fprintf('issue with profile %i\n',id_profile);
    end
end