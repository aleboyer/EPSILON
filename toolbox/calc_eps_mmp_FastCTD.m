function [Epsilon_mmp]=calc_eps_mmp_FastCTD(Profile,Epsilon,nb_profile,tscan,f_mmp,Sv)
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
T       = length(Profile.time_mmpshear);
df      = f_mmp(1);
%% define number of scan in the profile
Lscan   = tscan*2*f_mmp(end);
nbscan  = floor(T/Lscan);

%% we compute spectra on scan with 50% overlap
nbscan=2*nbscan-1;


%% define the index in the profile for each scan
total_indscan = arrayfun(@(x) (1+floor(Lscan/2)*(x-1):1+floor(Lscan/2)*(x-1)+Lscan-1),1:nbscan,'un',0);
total_w       = Epsilon.w; 
ind_downcast = find((total_w)>20);

Epsilon_mmp.indscan = total_indscan(ind_downcast);
Epsilon_mmp.w       = Epsilon.w; 
Epsilon_mmp.t       = Epsilon.t; % needed to compute Kvis 
Epsilon_mmp.s       = Epsilon.s; % needed to compute Kvis 
Epsilon_mmp.pr      = Epsilon.pr; % needed to compute Kvis 



%% read filters
H_MMP=load('toolbox/MMP/MMP4_shear_Htransfer.mat'); 
if nargin<6
    Sv=H_MMP.Sv;
    Sv=[1 1]*34.2;
end
h_freq=interp1(H_MMP.f,H_MMP.h_freq,f_mmp);

%h_freq=get_filters_MADRE(f);
%%
%fmax=30; %defined by get_coherence_epsi.m -> the WW is vibrating around 30 hz 
fmax=45; %defined by get_coherence_epsi.m -> the WW is vibrating around 30 hz 

if mod(nb_profile,10)==0
    fig=figure('Position',[100,100,1600,2000]);
%    v = VideoWriter(sprintf('../FIGURE/EPSILON/Epsilon_EPSIWW_P%i.avi',nb_profile));
    v = VideoWriter(sprintf('../FIGURE/EPSILON/MMPEpsilon_EPSISPROUL_P%i.avi',nb_profile));
    v.FrameRate=10;
    open(v)
end

% loop to estimate eps
for j=1:length(ind_downcast)
    speed=Epsilon_mmp.w(j)/100; % convert to m/s
    if speed > .2
        
        % select data per scan
        data1 = Profile.s1(Epsilon_mmp.indscan{j}(1:Lscan));
        data2 = Profile.s2(Epsilon_mmp.indscan{j}(1:Lscan));
        
        % compute kinematic viscosity
        kvis=nu(Epsilon.s(j),Epsilon_mmp.t(j),Epsilon.pr(j));
        % remove the mean to compute the fft
        data1 = data1-mean(data1);
        data2 = data2-mean(data2);
        [P1,f1] = pwelch(data1,Lscan,0,f_mmp,2*f_mmp(end));
        [P2,~]  = pwelch(data2,Lscan,0,f_mmp,2*f_mmp(end));
        
        if length(f1)~=length(f_mmp)
            warning('watch for the length of the scan %i',j);
        end
        % convert frequency to wavenumber
        k=f_mmp/speed;
        dk=df/speed;
        % evaluate and apply transfer function
        htotal = (Sv.*speed/(2*G)).^2 .* h_freq .* haf_oakey(f_mmp,speed);        % should add epsi filter
        Pvelk1  = (P1*speed) ./ htotal(:,1);                                      % vel spec as function of k
        Pvelk2  = (P2*speed) ./ htotal(:,2);                                      % vel spec as function of k
        Epsilon_mmp.Pshear1{j} = (2*pi*k).^2 .* Pvelk1;                          % shear spec  as function of k
        Epsilon_mmp.Pshear2{j} = (2*pi*k).^2 .* Pvelk2;                          % shear spec  as function of k
        Epsilon_mmp.k{j}      = k;                                             % k
        Epsilon_mmp.kvis{j}   =kvis;
        % calc epsilon by integrating to k with 90% variance of Panchev spec
        % unless spectrum is noisy at lower k.
        %
        % Set kmax for integration to highest bin below pump spike,
        % which is between 49 and 52 Hz in a 1024-pt spectrum
        Epsilon_mmp.kmax(j)=fmax/speed; % Lowest estimate below pump spike in 1024-pt record
        % Check that data window > 0.5 m, as needed for initial estimate
        if tscan*speed>0.5
            [Epsilon_mmp.epsilon1(j),Epsilon_mmp.kc1(j)]=eps1_mmp(k,Epsilon_mmp.Pshear1{j},kvis,speed,dk,Epsilon_mmp.kmax(j));
            [Epsilon_mmp.epsilon2(j),Epsilon_mmp.kc2(j)]=eps1_mmp(k,Epsilon_mmp.Pshear2{j},kvis,speed,dk,Epsilon_mmp.kmax(j));
            [Epsilon_mmp.kpan1{j},Epsilon_mmp.Ppan1{j}] = panchev(Epsilon_mmp.epsilon1(j),kvis);
            [Epsilon_mmp.kpan2{j},Epsilon_mmp.Ppan2{j}] = panchev(Epsilon_mmp.epsilon2(j),kvis);
            %
            if mod(nb_profile,10)==0
                % display each spectrum if displ_shear_spec='yes'
                
                ref1(1,1)=max(Epsilon_mmp.Pshear1{j}); ref1(2,1)=min(Epsilon_mmp.Pshear1{j});
                ref2(1,1)=max(Epsilon_mmp.Pshear2{j}); ref2(2,1)=min(Epsilon_mmp.Pshear2{j});
                kcref1(1,1)=Epsilon_mmp.kc1(j); kcref1(2,1)=Epsilon_mmp.kc1(j);
                kcref2(1,1)=Epsilon_mmp.kc2(j); kcref2(2,1)=Epsilon_mmp.kc2(j);
                kmaxref(1,1)=Epsilon_mmp.kmax(j); kmaxref(2,1)=Epsilon_mmp.kmax(j);
                loglog(Epsilon_mmp.k{j},Epsilon_mmp.Pshear1{j},'b');
                hold on
                loglog(Epsilon_mmp.k{j},Epsilon_mmp.Pshear2{j},'k');
                loglog(kcref1,ref1,'b'); %loglog(kcref,ref,'x')
                loglog(kcref2,ref2,'k'); %loglog(kcref,ref,'x')
                loglog(kmaxref,ref1,'g'); %loglog(kmaxref,ref,'+')
                title(['Profile ' int2str(nb_profile) 'Scan=' int2str(j) ', pr=' ...
                    num2str(Epsilon_mmp.pr(j)) ', epsilon1=' num2str(Epsilon_mmp.epsilon1(j)) ...
                    ', epsilon2=' num2str(Epsilon_mmp.epsilon2(j))]);
                xlabel('k (cpm)','fontsize',15)
                ylabel('\phi(k)','fontsize',15)
                set(gca,'fontsize',15)
                loglog(Epsilon_mmp.kpan1{j},Epsilon_mmp.Ppan1{j},'bo');
                loglog(Epsilon_mmp.kpan2{j},Epsilon_mmp.Ppan2{j},'ko');
                %ylim([0.1*min(Epsilon.Ppan1{j}) 1e3*max(Epsilon.Ppan1{j})])
                xlim([1e-1 1e3])
                ylim([1e-9 1e-1])
                hold off
                frame=getframe(fig);
                writeVideo(v,frame)
                clf;
            end
        end
    end
end% end nbscan, end profile processing
if mod(nb_profile,10)==0
    close(v)
    close all
end



