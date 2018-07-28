function [Epsilon]=calc_eps3G_epsi(Profile,tscan,f,Cs,Sv,nb_profile)
% calc_eps3G_mmp.m
%	Called from batchprocess4G_mmp to estimate epsilons, from v1/v2 spectra
%	integrated to cutoff chosen with help from Panchev spectra.
%	Intervals thinner than 0.5m or with w_eps<0.2m/s are skipped (eps1,2=NaN).
%	Does NOT re-compute if drop_flag=2 and eps<drop>.mat already exists.
%	REQUIRES the following, set up in batchprocess4G_mmp.m, setup_epschi3_mmp.m:
%		drop,procdata,cruise,mmpid, hfperscan,dt_hf,eps_step,FS_hf, f,df;
%		neps,cntr_scan,pr_eps,w_eps, t,s,kvis;  save_eps_spec,displ_shear_spec;
%	CREATES file eps<drop>.mat, saving:  epsilon,kc [neps-by-n_epsch matrices],
%		pr_eps,w_eps,t,s [vectors(1:neps)], epsch [1 and/or 2], eps_step
%		(Stopped saving eps1,eps2,kmax with this version);
%		Optionally, Psh1<drop>.mat,Psh2<drop>.mat save the shear spectra.
%	Returned for later plotting are:  epsilon, pr_eps, epsch, n_epsch;
%	Left for chi processing are:  eps_step,w_eps,cntr_scan, and many REQUIREDS
% 
%   June 2017: convert the previous routine to a function
%               input: - Profile structure- 
%                        Profile.shear, Profile.temp, Profile.sal Profile.Pressure, 
%                        ....,
%                      - tscan: length(second) of one scan. 1 scan give 1
%                      value of epsilon
%                      - f frequency array.
%                      - Cs (ohms)   impedance of the probe.
%                      - Sv (Volt/m) calibration parameter .
%
%
% REVISED  June-2017 by Arnaud Le Boyer
% REVISED  June-2000 by Dave Winkel,
%		from 02sep96 M.Gregg version calc_eps3_mmp.m

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
Profile = compute_fallrate(Profile);

%% define the index in the profile for each scan
total_indscan = arrayfun(@(x) (1+floor(Lscan/2)*(x-1):1+floor(Lscan/2)*(x-1)+Lscan-1),1:nbscan,'un',0);
total_w       = cellfun(@(x) nanmean(Profile.w(x)),total_indscan); 

ind_upcast = find((-total_w)>20);

Epsilon.indscan = total_indscan(ind_upcast);
Epsilon.w       = cellfun(@(x) nanmean(Profile.w(x)),total_indscan(ind_upcast)); 
Epsilon.t       = cellfun(@(x) nanmean(Profile.T(x)),total_indscan(ind_upcast)); % needed to compute Kvis 
Epsilon.s       = cellfun(@(x) nanmean(Profile.S(x)),total_indscan(ind_upcast)); % needed to compute Kvis 
Epsilon.pr      = cellfun(@(x) nanmean(Profile.P(x)),total_indscan(ind_upcast)); % needed to compute Kvis 


%% read filters
h_freq=get_filters_MADRE(f);
%%
fmax=30; %defined by get_coherence_epsi.m -> the WW is vibrating around 30 hz 

if mod(nb_profile,10)==0
    fig=figure('Position',[100,100,1600,2000]);
    v = VideoWriter(sprintf('../FIGURE/EPSILON/Epsilon_EPSIWW_P%i.avi',nb_profile));
    v.FrameRate=10;
    open(v)
end

% loop to estimate eps
for j=1:length(ind_upcast)
    speed=-Epsilon.w(j)/100; % convert to m/s
    if speed > .2
        
        % select data per scan
        data = Profile.Sensor3(Epsilon.indscan{j}(1:Lscan));
        
        % compute kinematic viscosity
        kvis=nu(Epsilon.s(j),Epsilon.t(j),Epsilon.pr(j));
        % remove the mean to compute the fft
        data = data-mean(data);
        [P,f1] = pwelch(data,Lscan,0,f,2*f(end));
        
        if length(f1)~=length(f)
            warning('watch for the length of the scan %i',j);
        end
        % convert frequency to wavenumber
        k=f/speed;
        dk=df/speed;
        % evaluate and apply transfer function
        htotal = (Sv*speed/(2*G))^2 .* h_freq.* haf_oakey(f,speed);        % should add epsi filter
        Pvelk  = (P*speed) ./ htotal;                                      % vel spec as function of k
        Epsilon.Pshear{j} = (2*pi*k).^2 .* Pvelk;                          % shear spec  as function of k
        Epsilon.k{j}      = k;                                             % k
        
        % calc epsilon by integrating to k with 90% variance of Panchev spec
        % unless spectrum is noisy at lower k.
        %
        % Set kmax for integration to highest bin below pump spike,
        % which is between 49 and 52 Hz in a 1024-pt spectrum
        Epsilon.kmax(j)=fmax/speed; % Lowest estimate below pump spike in 1024-pt record
        % Check that data window > 0.5 m, as needed for initial estimate
        if tscan*speed>0.5
            [Epsilon.epsilon(j),Epsilon.kc(j)]=eps1_mmp(k,Epsilon.Pshear{j},kvis,speed,dk,Epsilon.kmax(j));
            [Epsilon.kpan{j},Epsilon.Ppan{j}] = panchev(Epsilon.epsilon(j),kvis);
        end
        %
        if mod(nb_profile,10)==0
            % display each spectrum if displ_shear_spec='yes'
            
            ref(1,1)=max(Epsilon.Pshear{j}); ref(2,1)=min(Epsilon.Pshear{j});
            kcref(1,1)=Epsilon.kc(j); kcref(2,1)=Epsilon.kc(j);
            kmaxref(1,1)=Epsilon.kmax(j); kmaxref(2,1)=Epsilon.kmax(j);
            loglog(Epsilon.k{j},Epsilon.Pshear{j});
            hold on
            loglog(kcref,ref,'r'); %loglog(kcref,ref,'x')
            loglog(kmaxref,ref,'g'); %loglog(kmaxref,ref,'+')
            title(['Profile ' int2str(nb_profile) 'Scan=' int2str(j) ', pr=' ...
                num2str(Epsilon.pr(j)) ', epsilon=' num2str(Epsilon.epsilon(j)) ]);
            xlabel('k (cpm)','fontsize',15)
            ylabel('\phi(k)','fontsize',15)
            set(gca,'fontsize',15)
            loglog(Epsilon.kpan{j},Epsilon.Ppan{j},'bo');
            ylim([0.1*min(Epsilon.Ppan{j}) 1e3*max(Epsilon.Ppan{j})])
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



