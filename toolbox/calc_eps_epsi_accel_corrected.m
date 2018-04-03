function Epsilon=calc_eps_epsi_accel_corrected(Profile,tscan,f,Sv)
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
Lscan   = floor(tscan*2*f(end));
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
if Epsilon.s>1
    Epsilon.s=Epsilon.s./1000;
end


%% read filters
h_freq=get_filters_MADRE('MADRE2.1',f);
%%
fmax=20; %defined by get_coherence_epsi.m -> the WW is vibrating around 30 hz 

% loop to estimate eps
P11=zeros(length(ind_cast),Lscan);
P22=zeros(length(ind_cast),Lscan);
P33=zeros(length(ind_cast),Lscan);
P1a=zeros(length(ind_cast),Lscan);
P2a=zeros(length(ind_cast),Lscan);
P3a=zeros(length(ind_cast),Lscan);
P13a=zeros(length(ind_cast),Lscan);
P23a=zeros(length(ind_cast),Lscan);

Sh1=zeros(length(ind_cast),length(f));
Sh2=zeros(length(ind_cast),length(f));
P1=zeros(length(ind_cast),length(f));
P2=zeros(length(ind_cast),length(f));
P3=zeros(length(ind_cast),length(f));
H1=zeros(length(ind_cast),length(f));
H2=zeros(length(ind_cast),length(f));

for j=1:length(ind_cast)
    speed=Epsilon.w(j); % convert to m/s
    if speed > .2
        % select data per scan
        data1 = Profile.Sensor3(Epsilon.indscan{j}(1:Lscan));
        data2 = Profile.Sensor4(Epsilon.indscan{j}(1:Lscan));
        a1 = Profile.Sensor6(Epsilon.indscan{j}(1:Lscan));
        
        % compute kinematic viscosity
        kvis=nu(Epsilon.s(j),Epsilon.t(j),Epsilon.pr(j));
        % remove the mean to compute the fft
        data1 = data1-mean(data1);
        data2 = data2-mean(data2);
        a1p  = a1./speed;
        acc1  = (a1p-mean(a1p));
%        [P1,f1] = pwelch(data1,floor(Lscan),0,f,2*f(end));
        window = hanning(Lscan);
        wc2=1/mean(window.^2);            % window correction factor
        
        Am1 = window.'.*data1;
        Am2 = window.'.*data2;
        Aa1 = window.'.*acc1;
        
        P1a(j,:)  = fft(Am1);
        P2a(j,:)  = fft(Am2);
        P3a(j,:)  = fft(Aa1);
        P11(j,:) = abs(fft(Am1)).^2;
        P22(j,:) = abs(fft(Am2)).^2;
        P33(j,:) = abs(fft(Aa1)).^2;
        
        P11(j,:) = P11(j,:)./Lscan^2/df*wc2;
        P22(j,:) = P22(j,:)./Lscan^2/df*wc2;
        P33(j,:) = P33(j,:)./Lscan^2/df*wc2;
        
        coef = Lscan/f(end);               % duree d'un segment
        f1 = [0:floor(Lscan/2)-1 0 -floor(Lscan/2)+1:-1]' / coef;
        indf1=find(f1>=0);
        indf1=indf1(1:end-1);

        P1=2*P11(j,indf1).';
        P2=2*P22(j,indf1).';
        P3(j,:)=2*P33(j,indf1).'./h_freq.shear;

        if length(indf1)~=length(f)
            warning('watch for the length of the scan %i',j);
        end
        % convert frequency to wavenumber
        k=f/speed;
        k1=f1.'/speed;
        dk=df/speed;
        % evaluate and apply transfer function
        htotal = (Sv.*speed/(2*G)).^2 .* h_freq.shear .* haf_oakey(f,speed);        % should add epsi filter
        Pvel1  = P1 ./ htotal(:,1);                                      % vel spec as function of k
        Pvel2  = P2 ./ htotal(:,2);                                      % vel spec as function of k
        P1a(j,:)=P1a(j,:)./ sqrt([htotal(:,1).' flipud(htotal(:,1)).' ]); 
        P2a(j,:)=P2a(j,:)./ sqrt([htotal(:,2).' flipud(htotal(:,2)).' ]);
        P3a(j,:)=P3a(j,:)./ [h_freq.shear.' h_freq.shear.'];
        
        Sh1(j,:)=(2*pi*k).^2 .* Pvel1;
        Sh2(j,:)=(2*pi*k).^2 .* Pvel2;
        P1a(j,:)=(2*pi*k1) .* P1a(j,:);
        P2a(j,:)=(2*pi*k1) .* P2a(j,:);
        P13a(j,:)=conj(P1a(j,:)).*P3a(j,:)./(Lscan)^2/df*wc2;
        P23a(j,:)=conj(P2a(j,:)).*P3a(j,:)./(Lscan)^2/df*wc2;
        P13=2*P13a(j,indf1);
        P23=2*P23a(j,indf1);
        H1(j,:)= P13./P3(j,:); %transfer function between measured (mm) shear and acceleration (aa)
        H2(j,:)= P23./P3(j,:); %transfer function between measured (mm) shear and acceleration (aa)
    end
end
       
mH1=smoothdata(H1,1,'movmean',20);
mH2=smoothdata(H2,1,'movmean',20);

for j=1:length(ind_cast)
    speed=Epsilon.w(j); % convert to m/s
    if speed > .2

        Sh1a(j,:) = Sh1(j,:) - abs(mH1(j,:)).^2.*P3(j,:);  % Css(f) = Cmm(f) - |H(f)|^2*Caa(f)
        Sh2a(j,:) = Sh2(j,:) - abs(mH2(j,:)).^2.*P3(j,:);  % Css(f) = Cmm(f) - |H(f)|^2*Caa(f)
        
        Epsilon.Pshear1{j} = Sh1a*speed;                          % shear spec  as function of k
        Epsilon.Pshear2{j} = Sh2a*speed;                          % shear spec  as function of k

        Epsilon.k{j}      = k;                                             % k
        Epsilon.kvis{j}   = kvis;
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



