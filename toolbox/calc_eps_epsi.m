function Epsilon=calc_eps_epsi(Profile,tscan,f,Sv)
% calc_eps_mmp.m
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
%                        we are using Profile.shear, Profile.w in (m/s), 
%                        Profile.T,Profile.S (to compute nu)
%                      - tscan: length(second) of one scan. 1 scan give 1
%                      value of epsilon
%                      - f frequency array.
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

%% define the index in the profile for each scan
total_indscan = arrayfun(@(x) (1+floor(Lscan/2)*(x-1):1+floor(Lscan/2)*(x-1)+Lscan-1),1:nbscan,'un',0);
total_w       = cellfun(@(x) nanmean(Profile.w(x)),total_indscan); 

ind_cast = find((total_w)>.2);

Epsilon.indscan = total_indscan(ind_cast);
Epsilon.w       = cellfun(@(x) nanmean(Profile.w(x)),total_indscan(ind_cast)); 
Epsilon.t       = cellfun(@(x) nanmean(Profile.T(x)),total_indscan(ind_cast)); % needed to compute Kvis 
Epsilon.s       = cellfun(@(x) nanmean(Profile.S(x)),total_indscan(ind_cast)); % needed to compute Kvis 
Epsilon.pr      = cellfun(@(x) nanmean(Profile.P(x)),total_indscan(ind_cast)); % needed to compute Kvis 


%% read filters
h_freq=get_filters_MADRE('MADRE2.1',f);
%%
fmax=30; %defined by get_coherence_epsi.m -> the WW is vibrating around 30 hz 

% loop to estimate eps
for j=1:length(ind_cast)
    speed=Epsilon.w(j); % convert to m/s
    if speed > .2
        % select data per scan
        data1 = Profile.Sensor3(Epsilon.indscan{j}(1:Lscan));
        data2 = Profile.Sensor4(Epsilon.indscan{j}(1:Lscan));
        
        % compute kinematic viscosity
        kvis=nu(Epsilon.s(j),Epsilon.t(j),Epsilon.pr(j));
        % remove the mean to compute the fft
        data1 = data1-mean(data1);
        data2 = data2-mean(data2);
        [P1,f1] = pwelch(data1,Lscan,0,f,2*f(end));
        [P2,~]  = pwelch(data2,Lscan,0,f,2*f(end));
        
        if length(f1)~=length(f)
            warning('watch for the length of the scan %i',j);
        end
        % convert frequency to wavenumber
        k=f/speed;
        dk=df/speed;
        % evaluate and apply transfer function
        htotal1 = (Sv(1).*speed/(2*G)).^2 .* h_freq.shear .* haf_oakey(f,speed);        % should add epsi filter
        htotal2 = (Sv(2).*speed/(2*G)).^2 .* h_freq.shear .* haf_oakey(f,speed);        % should add epsi filter
        Pvelk1  = (P1*speed) ./ htotal1;                                      % vel spec as function of k
        Pvelk2  = (P2*speed) ./ htotal2;                                      % vel spec as function of k
        Epsilon.Pshear1{j} = (2*pi*k).^2 .* Pvelk1;                          % shear spec  as function of k
        Epsilon.Pshear2{j} = (2*pi*k).^2 .* Pvelk2;                          % shear spec  as function of k
        Epsilon.k{j}       = k;                                             % k
        Epsilon.kvis{j}    = kvis;
        % calc epsilon by integrating to k with 90% variance of Panchev spec
        % unless spectrum is noisy at lower k.
        %
        % Set kmax for integration to highest bin below pump spike,
        % which is between 49 and 52 Hz in a 1024-pt spectrum
        Epsilon.kmax(j)=fmax/speed; % Lowest estimate below pump spike in 1024-pt record
        % Check that data window > 0.5 m, as needed for initial estimate
        if tscan*speed>0.5
            [Epsilon.epsilon1(j),Epsilon.kc1(j)]=eps1_mmp(k,Epsilon.Pshear1{j},kvis,speed,dk,Epsilon.kmax(j));
            [Epsilon.epsilon2(j),Epsilon.kc2(j)]=eps1_mmp(k,Epsilon.Pshear2{j},kvis,speed,dk,Epsilon.kmax(j));
            [Epsilon.kpan1{j},Epsilon.Ppa1{j}] = panchev(Epsilon.epsilon1(j),kvis);
            [Epsilon.kpan2{j},Epsilon.Ppa2{j}] = panchev(Epsilon.epsilon2(j),kvis);
        end
        %
    end
end% end nbscan, end profile processing



