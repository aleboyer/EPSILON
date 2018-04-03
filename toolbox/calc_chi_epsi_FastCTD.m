function [chi]=calc_chi_epsi_FastCTD(Profile,tscan,f,nb_profile)
% calc_chi3G_mmp.m
%	Called from batchprocess_epsi_WW to estimate chis, from th1/th2 spectra
%	integrated to cutoff based on voltage noise floor.
%	Intervals thinner than 0.5m or with w_eps<0.2m/s are skipped (chi1,2=NaN).
%	Does NOT re-compute if drop_flag=2 and chi<drop>.mat already exists.
%	REQUIRES the following, set up in batchprocess4G_mmp.m, setup_epschi3_mmp.m:
%		drop,procdata,cruise,mmpid, hfperscan,dt_hf,eps_step,FS_hf, f,df;
%		nchi,cntr_scan,pr_chi,w_eps, t,s,ktemp;  save_chi_spec,displ_chi_spec;
%	CREATES file chi<drop>.mat, saving:  chi,kcth [nchi-by-n_tlch matrices],
%		pr_chi,w_eps,t,s [vectors(1:nchi)], tlch [1 and/or 2], eps_step
%		Optionally, Ptg1<drop>.mat,Ptg2<drop>.mat save the temp.grad spectra.
%	Returned for later plotting are:  chi, pr_chi, tlch, n_tlch;
% REVISED  Sept-2001 by Dave Winkel,
%		from 02sep96 M.Gregg version calc_chi2_mmp.m
%
% REVISED  June-2017 by Arnaud Le Boyer

%% Length of the Profile
T       = length(Profile.nbsample);
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
if sum(abs(Profile.w)<10)==0
    Profile.w=Profile.w/100;
end

%% define the index in the profile for each scan
total_indscan = arrayfun(@(x) (1+floor(Lscan/2)*(x-1):1+floor(Lscan/2)*(x-1)+Lscan-1),1:nbscan,'un',0);
total_w       = cellfun(@(x) nanmean(Profile.w(x)),total_indscan); 

ind_downcast = find((total_w)>.20);

chi.indscan = total_indscan(ind_downcast);
chi.w       = cellfun(@(x) nanmean(Profile.w(x)),total_indscan(ind_downcast)); 
chi.t       = cellfun(@(x) nanmean(Profile.T(x)),total_indscan(ind_downcast)); % needed to compute Kvis 
chi.s       = cellfun(@(x) nanmean(Profile.S(x)),total_indscan(ind_downcast)); % needed to compute Kvis 
chi.pr      = cellfun(@(x) nanmean(Profile.P(x)),total_indscan(ind_downcast)); % needed to compute Kvis 

% on line calibration of the FPO7s
CALFPO7_1=polyfit(Profile.Sensor1,Profile.T,3);
CALFPO7_2=polyfit(Profile.Sensor2,Profile.T,3);


%% read filters
h_freq=get_filters_MADRE('MADRE2.1',f);
%
Emp_Corr_fac = 1; % empirical correction factor, until figured out


if mod(nb_profile,1)==0
    fig=figure('Position',[100,100,1600,2000]);
    v = VideoWriter(sprintf('../PLUMEX18/CHI_EPSIWW_P%i.avi',nb_profile));
    v.FrameRate=10;
    open(v)
end

% loop to estimate eps
for j=1:length(ind_downcast)
    speed=chi.w(j);%/100; % convert to m/s
    if speed > .2
        

        % select data per scan
        data1 = polyval(CALFPO7_1,Profile.Sensor1(chi.indscan{j}(1:Lscan)));
        data2 = polyval(CALFPO7_2,Profile.Sensor2(chi.indscan{j}(1:Lscan)));
        
        ktemp=kt(chi.s(j),chi.t(j),chi.pr(j));
        % compute kinematic viscosity
        % remove the mean to compute the fft
        data1 = data1-mean(data1);
        data2 = data2-mean(data2);
        [P1,f1] = pwelch(data1,Lscan,0,f,2*f(end));
        [P2,~]  = pwelch(data2,Lscan,0,f,2*f(end));
        
        if length(f1)~=length(f)
            warning('watch for the length of the scan %i',j);
        end
        
        % vertical wavenumber spectrum of the temperature on one scan
        chi.Ptempk1{j} = Emp_Corr_fac * (P1*speed)./h_freq.FPO7(speed);
        chi.Ptempk2{j} = Emp_Corr_fac * (P2*speed)./h_freq.FPO7(speed);
        
        
        % convert temperature gradient
        k  = f/speed;
        dk = df/speed;
        chi.Ptgradk1{j}= (2*pi.*k).^2 .* chi.Ptempk1{j}(:) ;
        chi.Ptgradk2{j}= (2*pi.*k).^2 .* chi.Ptempk2{j}(:) ;
        chi.k{j}      = k;                                             % k
        chi.ktemp{j}  = ktemp;
        % compute cut off frequency
        displ_chi_spec='no';
        chi.fc_index1(j)=FPO7_cutoff(f,P1,displ_chi_spec,chi.pr(j),j,speed);
        chi.fc_index2(j)=FPO7_cutoff(f,P2,displ_chi_spec,chi.pr(j),j,speed);
        chi.chi1(j)=6*ktemp*dk.*sum(chi.Ptgradk1{j}(1:chi.fc_index1(j)));
        chi.chi2(j)=6*ktemp*dk.*sum(chi.Ptgradk2{j}(1:chi.fc_index2(j)));
        
        %TODO: compute dn2 and dthetadz
        %eps_chi(j,i) = abs(n2(idn2)*chi(j,i)/(2*0.2*dthetadz(idn2)^2));

        %
        if mod(nb_profile,1)==0
            % display each spectrum if displ_shear_spec='yes'
            
            ref1(1,1)=max(chi.Ptgradk1{j}); ref1(2,1)=min(chi.Ptgradk1{j});
            ref2(1,1)=max(chi.Ptgradk2{j}); ref2(2,1)=min(chi.Ptgradk2{j});
            kcref1(1,1)=chi.fc_index1(j); kcref1(2,1)=chi.fc_index1(j);
            kcref2(1,1)=chi.fc_index2(j); kcref2(2,1)=chi.fc_index2(j);
            loglog(chi.k{j},chi.Ptgradk1{j},'b');
            hold on
            loglog(chi.k{j},chi.Ptgradk2{j},'k');
            loglog(kcref1,ref1,'b.-'); %loglog(kcref,ref,'x')
            loglog(kcref2,ref2,'k.-'); %loglog(kcref,ref,'x')
            title(['Profile ' int2str(nb_profile) 'Scan=' int2str(j) ', pr=' ...
                num2str(chi.pr(j)) ', chi1=' num2str(chi.chi1(j)) ', chi2=' num2str(chi.chi2(j)) ]);
            xlabel('k (cpm)','fontsize',15)
            ylabel('\phi_{TG}(k)','fontsize',15)
            set(gca,'fontsize',15)
            ylim([min(chi.Ptgradk2{j}) max(chi.Ptgradk2{j})])
            xlim([1e-3 1e3])
            hold off
            frame=getframe(fig);
            writeVideo(v,frame)
            clf;
        end
    end
end% end nbscan, end profile processing
if mod(nb_profile,1)==0
    close(v)
    close all
end



