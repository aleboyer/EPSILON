function [chi]=calc_chi_epsi(Profile,tscan,f,nb_profile)
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
total_indscan = arrayfun(@(x) (1+floor(Lscan/2)*(x-1):1+floor(Lscan/2)*(x-1)+Lscan-1),1:nbscan,'un',0);
total_w       = cellfun(@(x) nanmean(Profile.w(x)),total_indscan); 

ind_upcast = find((-total_w)>20);

chi.indscan = total_indscan(ind_upcast);
chi.w       = cellfun(@(x) nanmean(Profile.w(x)),total_indscan(ind_upcast)); 
chi.t       = cellfun(@(x) nanmean(Profile.T(x)),total_indscan(ind_upcast)); % needed to compute Kvis 
chi.s       = cellfun(@(x) nanmean(Profile.S(x)),total_indscan(ind_upcast)); % needed to compute Kvis 
chi.pr      = cellfun(@(x) nanmean(Profile.P(x)),total_indscan(ind_upcast)); % needed to compute Kvis 


%% read filters
h_freq=get_filters_MADRE(f);
%
Emp_Corr_fac = 1/30; % empirical correction factor, until figured out
% temporay calibration of the FPO7 probe computed in epsiWW_cal_fpo7.m
% maybe an option would be to recalibrate fpo7 for each profile ?? 
load('tempo_CALFPO7.mat','CALFPO7')

if mod(nb_profile,10)==0
    fig=figure('Position',[100,100,1600,2000]);
    v = VideoWriter(sprintf('../FIGURE/CHI/CHI_EPSIWW_P%i.avi',nb_profile));
    v.FrameRate=10;
    open(v)
end

% loop to estimate eps
for j=1:length(ind_upcast)
    speed=-chi.w(j)/100; % convert to m/s
    if speed > .2
        
        % select data per scan
        data = polyval(CALFPO7,Profile.Sensor1(chi.indscan{j}(1:Lscan)));
        
        ktemp=kt(chi.s(j),chi.t(j),chi.pr(j));
        % compute kinematic viscosity
        % remove the mean to compute the fft
        data = data-mean(data);
        [P,f1] = pwelch(data,Lscan,0,f,2*f(end));
        
        if length(f1)~=length(f)
            warning('watch for the length of the scan %i',j);
        end
        
        % vertical wavenumber spectrum of the temperature on one scan
        chi.Ptempk{j} = Emp_Corr_fac * (P*speed)./h_freq.FPO7(speed);
        
        
        % convert temperature gradient
        k=f/speed;
        dk=df/speed;
        chi.Ptgradk{j}= (2*pi.*k).^2 .* chi.Ptempk{j}(:) ;
        chi.k{j}      = k;                                             % k
        chi.ktemp{j}  = ktemp;
        % compute cut off frequency
        displ_chi_spec='no';
        chi.fc_index(j)=FPO7_cutoff(f,P,displ_chi_spec,chi.pr(j),j,speed);
        chi.chi(j)=6*ktemp*dk.*sum(chi.Ptgradk{j}(1:chi.fc_index(j)));
        
        %TODO: compute dn2 and dthetadz
        %eps_chi(j,i) = abs(n2(idn2)*chi(j,i)/(2*0.2*dthetadz(idn2)^2));

        %
        if mod(nb_profile,10)==0
            % display each spectrum if displ_shear_spec='yes'
            
            ref(1,1)=max(chi.Ptgradk{j}); ref(2,1)=min(chi.Ptgradk{j});
            kcref(1,1)=chi.fc_index(j); kcref(2,1)=chi.fc_index(j);
            loglog(chi.k{j},chi.Ptgradk{j});
            hold on
            loglog(kcref,ref,'r'); %loglog(kcref,ref,'x')
            title(['Profile ' int2str(nb_profile) 'Scan=' int2str(j) ', pr=' ...
                num2str(chi.pr(j)) ', chi=' num2str(chi.chi(j)) ]);
            xlabel('k (cpm)','fontsize',15)
            ylabel('\phi_{TG}(k)','fontsize',15)
            set(gca,'fontsize',15)
            ylim([min(chi.Ptgradk{j}) max(chi.Ptgradk{j})])
            hold off
            frame=getframe(fig);
            writeVideo(v,frame)
            clf;
        end
    end
end% end nbscan, end profile processing
if mod(nb_profile,10)==0
    close(v)
    close all
end



