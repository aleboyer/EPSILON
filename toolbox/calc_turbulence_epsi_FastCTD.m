function [MS]=calc_turbulence_epsi_FastCTD(Profile,tscan,f,Sv)
% calc_eps_mmp.m
%   
%	REQUIRES the following, set up in batchprocess4G_mmp.m, setup_epschi3_mmp.m:
%		drop,procdata,cruise,mmpid, hfperscan,dt_hf,eps_step,FS_hf, f,df;
%		neps,cntr_scan,pr_eps,w_eps, t,s,kvis;  save_eps_spec,displ_shear_spec;
%	CREATES file eps<drop>.mat, saving:  
%       MS structure for Micro Structure. Inside MS you ll find
%       temperature spectra in degC Hz^-1
%       Horizontal  velocity spectra in m^2/s^-2 Hz^-1
%       Acceleration/speed spectra in s^-1 Hz^-1 
%
%
%
% REVISED  June-2017 by Arnaud Le Boyer
% REVISED  June-2000 by Dave Winkel,
%		from 02sep96 M.Gregg version calc_eps3_mmp.m

%% Gravity  ... of the situation :)
G       = 9.81;

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
Profile.w=Profile.w/100;

% on line calibration of the FPO7s
CALFPO7_1=polyfit(Profile.Sensor1,Profile.T,3);
CALFPO7_2=polyfit(Profile.Sensor2,Profile.T,3);
%CALFPO7_2=polyfit(Profile.Sensor2,Profile.T,1);


%% define the index in the profile for each scan
total_indscan = arrayfun(@(x) (1+floor(Lscan/2)*(x-1):1+floor(Lscan/2)*(x-1)+Lscan-1),1:nbscan,'un',0);
total_w       = cellfun(@(x) nanmean(Profile.w(x)),total_indscan); 

ind_downcast = find((total_w)>.20);

MS.indscan   = total_indscan(ind_downcast);
MS.nbscan    = nbscan;
MS.fmax      = 45; % arbitrary cut off frequency usually extract from coherence spectra shear/accel 
MS.nbchannel = nbscan;

MS.w       = cellfun(@(x) nanmean(Profile.w(x)),total_indscan(ind_downcast)); 
MS.t       = cellfun(@(x) nanmean(Profile.T(x)),total_indscan(ind_downcast)); % needed to compute Kvis 
MS.s       = cellfun(@(x) nanmean(Profile.S(x)),total_indscan(ind_downcast)); % needed to compute Kvis 
MS.pr      = cellfun(@(x) nanmean(Profile.P(x)),total_indscan(ind_downcast)); % needed to compute Kvis 


% loop to estimate eps
% get and split the data
data(1,:,:) = cell2mat(cellfun(@(x) polyval(CALFPO7_1,Profile.Sensor1(x)),MS.indscan,'un',0).');
data(2,:,:) = cell2mat(cellfun(@(x) polyval(CALFPO7_2,Profile.Sensor2(x)),MS.indscan,'un',0).');
data(3,:,:) = cell2mat(cellfun(@(x) Profile.Sensor3(x),MS.indscan,'un',0).');
data(4,:,:) = cell2mat(cellfun(@(x) Profile.Sensor4(x),MS.indscan,'un',0).');
% acceleration divided by speed to get the same units (s-1)
data(5,:,:) = cell2mat(cellfun(@(x,y) Profile.Sensor6(x)./y,MS.indscan,num2cell(MS.w),'un',0).');
data(6,:,:) = cell2mat(cellfun(@(x,y) Profile.Sensor7(x)./y,MS.indscan,num2cell(MS.w),'un',0).');
data(7,:,:) = cell2mat(cellfun(@(x,y) Profile.Sensor8(x)./y,MS.indscan,num2cell(MS.w),'un',0).');

% compute kinematic viscosity
MS.kvis=nu(MS.s,MS.t,MS.pr);
MS.ktemp=kt(MS.s,MS.t,MS.pr).';

% remove the mean to compute the fft

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
    shiftdim(repmat(ones(nbscan,1)*h_freq.electAccel,[1,1,3]),2);

%% correct transfert functions for shear spectra

TF1 =@(x) (Sv.'.*x/(2*G)).^2 .* h_freq.shear .* haf_oakey(f1,x);     % should add epsi filter
TFshear=cell2mat(cellfun(@(x) TF1(x),num2cell(MS.w),'un',0).');
TFshear=reshape(TFshear,[2,nbscan,Lf1]);
P11(3:4,:,:) = P11(3:4,:,:) ./ TFshear;      % vel frequency spectra m^2/s^-2 Hz^-1

%% correct transfert functions for temperature spectra
Emp_Corr_fac=1;
TFtemp=cell2mat(cellfun(@(x) h_freq.FPO7(x),num2cell(MS.w),'un',0).');
TFtemp=shiftdim(repmat(TFtemp,[1,1,2]),2);
P11(1:2,:,:) = Emp_Corr_fac * P11(1:2,:,:)./TFtemp;


% convert frequency to wavenumber
k=cell2mat(cellfun(@(x) f1/x, num2cell(MS.w),'un',0).');
dk=cell2mat(cellfun(@(x) df/x, num2cell(MS.w),'un',0));

dk_all=nanmin(nanmean(diff(k,1,2),2));
k_all=nanmin(k(:)):dk_all:nanmax(k(:));

% temperature, vel and accell spec as function of k
nbchannel=7;
P11k  = P11.* shiftdim(repmat(ones(nbchannel,1)*MS.w,[1,1,Lf1]),3);   

MS.f  =  f1;
MS.k  = k_all;
MS.Pf =  P11;
% Set kmax for integration to highest bin below pump spike,
% which is between 49 and 52 Hz in a 1024-pt spectrum
MS.kmax=MS.fmax./MS.w; % Lowest estimate below pump spike in 1024-pt record



% calc epsilon by integrating to k with 90% variance of Panchev spec
% unless spectrum is noisy at lower k.
% Check that data window > 0.5 m, as needed for initial estimate
for j=1:nbscan
    MS.PphiT_k(j,:,1)  = (2*pi*k_all).^2 .* interp1(k(j,:),squeeze(P11k(1,j,:)),k_all);        % T1_k spec  as function of k
    MS.PphiT_k(j,:,2)  = (2*pi*k_all).^2 .* interp1(k(j,:),squeeze(P11k(2,j,:)),k_all);        % T2_k spec  as function of k
    MS.Pshear_k(j,:,1) = (2*pi*k_all).^2 .* interp1(k(j,:),squeeze(P11k(3,j,:)),k_all);        % shear spec  as function of k
    MS.Pshear_k(j,:,2) = (2*pi*k_all).^2 .* interp1(k(j,:),squeeze(P11k(4,j,:)),k_all);        % shear spec  as function of k
    MS.Paccel_k(j,:,1)     =  interp1(k(j,:),squeeze(P11k(5,j,:)),k_all);        % Accel 1 spec  as function of k
    MS.Paccel_k(j,:,2)     =  interp1(k(j,:),squeeze(P11k(6,j,:)),k_all);        % Accel 2 spec  as function of k
    MS.Paccel_k(j,:,3)     =  interp1(k(j,:),squeeze(P11k(7,j,:)),k_all);        % Accel 3 spec  as function of k
    
    % compute epsilon in eps1_mmp
    [MS.epsilon(j,1),MS.kc(j,1)]=eps1_mmp(k_all,MS.Pshear_k(j,:,1),MS.kvis(j),dk(j),MS.kmax(j)); 
    [MS.epsilon(j,2),MS.kc(j,2)]=eps1_mmp(k_all,MS.Pshear_k(j,:,2),MS.kvis(j),dk(j),MS.kmax(j));
    [kpan,Ppan] = panchev(MS.epsilon(j,1),MS.kvis(j));
    MS.Ppan(j,:,1)=interp1(kpan,Ppan,k_all);
    [kpan,Ppan] = panchev(MS.epsilon(j,2),MS.kvis(j));
    MS.Ppan(j,:,2)=interp1(kpan,Ppan,k_all);

    %% TODO check to see if we need P1 in V^2 Hz^{-1} or if we can change it to degC Hz^{-1}
    MS.fc_index(j,1)=FPO7_cutoff(f1,squeeze(P11(1,j,:).*TFtemp(1,j,:)));
    MS.fc_index(j,2)=FPO7_cutoff(f1,squeeze(P11(2,j,:).*TFtemp(2,j,:)));
    MS.chi(j,1)=6*MS.ktemp(j)*dk(j).*nansum(MS.PphiT_k(j,1:MS.fc_index(j,1),1));
    MS.chi(j,2)=6*MS.ktemp(j)*dk(j).*nansum(MS.PphiT_k(j,1:MS.fc_index(j,2),2));
end



